"""CARSPA: Carpinteri & Spagnoli (2001), Int J Fatigue 23:135-145, doi:10.1016/S0142-1123(00)00075-X,
with the Euler-angle averaging of Carpinteri et al. (1999), Int J Fatigue 21:83-88,
doi:10.1016/S0142-1123(98)00046-2, and the notation of Carpinteri, Brighenti & Spagnoli (2000),
Fatigue Fract Engng Mater Struct 23:355-364, doi:10.1046/j.1460-2695.2000.00265.x. All three
checked against the PDFs supplied by the maintainer.

Tables 2-4 (p.142) print, for 30 Nishihara-Kawamoto bending-torsion limits, the angle eta_cal
between the specimen axis and 1^ (c = 0.5) and the error index I of Eq. (30).

Findings:
- the criterion step: with the printed eta_cal, Eqs. (24)-(27) and w = cos(delta) 1^ +
  sin(delta) 3^ give all 30 printed I within their rounding (and the rounding of eta_cal). The
  opposite rotation misses by up to 15 points;
- the averaging step: 27 of 30 eta_cal within 1.5 deg (the table tends to lie 0.5-1.4 deg above).
  For the three out-of-phase (beta = 90) rows with tau_a/sigma_a = 0.2, the table prints 0 where
  Eqs. (22)-(23) with the reduction of the 1999 paper give 2.3-3.5 deg: off-peak instants with
  sigma_1 >= c sigma_af have principal directions on both sides of the axis, which the 2nd
  reduction folds to one side. No time step, threshold or exponent gives 0 here and also the
  printed 12 deg of hard steel test 11 (strict xfails, unexplained);
- Findley, Matake and McDiarmid columns are transcribed for reference only.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.cp_methods.car_spa import (
    CarSpa,
    carspa_critical_plane,
    carspa_delta,
    euler_matrix,
    plane_stresses,
    reduced_euler_angles,
)
from multifat.models.errors import LoadingScopeError, MaterialScopeError, MissingParamsError
from multifat.models.registry import get_models

TABLE = load_table("carpinteri2001_tables2_4")
ETA_UNEXPLAINED = {("hard", 10), ("mild", 8), ("cast", 6)}
ETA_TOL = 1.5  # deg


def harmonics(sigma_a, tau_a, beta_deg):
    """sigma_l = sigma_a sin(wt) on yy (axis Y), tau = tau_a sin(wt - beta) on xy (Eq. 1)."""
    b = np.radians(beta_deg)
    S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3))
    S_sin[1, 1] = sigma_a
    S_sin[0, 1] = S_sin[1, 0] = tau_a * np.cos(b)
    S_cos[0, 1] = S_cos[1, 0] = -tau_a * np.sin(b)
    return np.zeros((3, 3)), S_sin, S_cos


def params(r):
    return {"Sigm1": r.sigma_af, "Taum1": r.tau_af, "m_sigma": r.m_sigma}


def error_index_at_eta(r, eta_deg):
    """I of Eq. (30) on the plane built from 1^ at eta from the axis, 2^ = -Z (2000, case c)."""
    e = np.radians(eta_deg)
    e1, e2 = np.array([np.sin(e), np.cos(e), 0.0]), np.array([0.0, 0.0, -1.0])
    e3 = np.cross(e1, e2)
    d = carspa_delta(r.sigma_af, r.tau_af)
    n_m, n_a, c_a = plane_stresses(
        *harmonics(r.sigma_a, r.tau_a, r.beta), np.cos(d) * e1 + np.sin(d) * e3
    )
    s_eq = np.sqrt((n_m + n_a) ** 2 + (r.sigma_af / r.tau_af * c_a) ** 2)
    return 100.0 * (s_eq - r.sigma_af) / r.sigma_af


def computed_eta(r):
    res = carspa_critical_plane(*harmonics(r.sigma_a, r.tau_a, r.beta), params(r))
    return float(np.degrees(np.arccos(abs(res["axes"][1, 0])))), res


def _ids(r):
    return f"{r.material}-{r.test}"


@pytest.mark.parametrize("r", list(TABLE.itertuples()), ids=_ids)
def test_error_index_at_printed_eta(r):
    """Eqs. (24)-(27) reproduce I from the printed eta_cal (both printed to integers)."""
    slope = abs(error_index_at_eta(r, r.eta_cal + 0.5) - error_index_at_eta(r, r.eta_cal - 0.5))
    assert abs(error_index_at_eta(r, r.eta_cal) - r.I_present) <= 0.5 + 0.5 * slope + 1e-9


def test_rotation_sense_is_fixed_by_the_tables():
    """w rotated from 1^ towards -3^ instead misses hard steel test 11 by about 15 points."""
    r = TABLE[(TABLE.material == "hard") & (TABLE.test == 11)].iloc[0]
    e = np.radians(r.eta_cal)
    e1 = np.array([np.sin(e), np.cos(e), 0.0])
    e3 = np.cross(e1, [0.0, 0.0, -1.0])
    d = carspa_delta(r.sigma_af, r.tau_af)
    n_m, n_a, c_a = plane_stresses(
        *harmonics(r.sigma_a, r.tau_a, r.beta), np.cos(d) * e1 - np.sin(d) * e3
    )
    s_eq = np.sqrt((n_m + n_a) ** 2 + (r.sigma_af / r.tau_af * c_a) ** 2)
    assert 100.0 * (s_eq - r.sigma_af) / r.sigma_af < -14.0


def _eta_rows():
    reason = (
        "Table prints eta_cal = 0 for this beta = 90 deg, tau_a/sigma_a = 0.2 row; Eqs. (22)-(23) "
        "with the 1999 reduction give 2.3-3.5 deg (unexplained, see module docstring)"
    )
    return [
        (
            pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=reason), id=_ids(r))
            if (r.material, r.test) in ETA_UNEXPLAINED
            else pytest.param(r, id=_ids(r))
        )
        for r in TABLE.itertuples()
    ]


@pytest.mark.parametrize("r", _eta_rows())
def test_fracture_plane_angle(r):
    eta, _ = computed_eta(r)
    assert abs(eta - r.eta_cal) <= ETA_TOL


@pytest.mark.parametrize(
    "r", [r for r in TABLE.itertuples() if (r.material, r.test) not in ETA_UNEXPLAINED], ids=_ids
)
def test_error_index_end_to_end(r):
    """Full method: I within the rounding of I plus the effect of the eta_cal difference."""
    eta, res = computed_eta(r)
    slope = abs(error_index_at_eta(r, eta + 0.5) - error_index_at_eta(r, eta - 0.5))
    tol = 0.5 + slope * (abs(eta - r.eta_cal) + 0.5)
    assert abs(100.0 * (res["sigma_eq"] - r.sigma_af) / r.sigma_af - r.I_present) <= tol + 1e-9


@pytest.mark.parametrize("gamma", [0.0, 12.0, 30.0, 44.0])
def test_plane_stress_cases_of_the_2000_paper(gamma):
    """Reduced Euler angles of cases (a1)-(c2), Carpinteri et al. (2000) p.358-359."""
    g = np.radians(gamma)
    X, Y, Z = np.eye(3)

    def at(angle):  # in-plane unit vector at angle from X
        return np.cos(angle) * X + np.sin(angle) * Y

    def red(e1, e2, e3):
        return np.degrees(reduced_euler_angles(np.stack([e1, e2, e3], axis=-1)))

    # (a1)/(a2): sigma_3 = 0 along Z, 1 at gamma (or 90 - gamma) from X
    np.testing.assert_allclose(red(at(g), at(g + np.pi / 2), Z), [gamma, 90, 90], atol=1e-9)
    np.testing.assert_allclose(
        red(at(np.pi / 2 - g), at(np.pi - g), Z), [90 - gamma, 90, 90], atol=1e-9
    )
    # (c1)/(c2): sigma_2 = 0 along Z, phi = gamma (or 90 - gamma), theta = 90, psi = 0
    np.testing.assert_allclose(red(at(g), Z, at(g + np.pi / 2)), [gamma, 90, 0], atol=1e-9)
    np.testing.assert_allclose(
        red(at(np.pi / 2 - g), Z, at(np.pi - g)), [90 - gamma, 90, 0], atol=1e-9
    )
    # (b1): sigma_1 = 0 along Z, theta = 0 and phi + psi = gamma
    phi, theta, psi = red(Z, at(g), at(g + np.pi / 2))
    assert theta == pytest.approx(0.0, abs=1e-9) and phi + psi == pytest.approx(gamma, abs=1e-9)


def test_reduction_round_trip():
    """Matrix (2) of the reduced angles spans the same axes, for random frames."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        f, t, p = reduced_euler_angles(Q)
        assert (
            0 <= f <= np.pi / 2 + 1e-9 and 0 <= t <= np.pi / 2 + 1e-9 and abs(p) <= np.pi / 2 + 1e-9
        )
        A = euler_matrix(f, t, p)
        # direction 1 is folded into the first octant: the same axis up to mirror images
        np.testing.assert_allclose(np.abs(A[:, 2]), np.abs(Q[:, 0]), atol=1e-9)


@pytest.mark.parametrize("beta", [0, 30, 60, 90])
def test_mirror_invariance(beta):
    """tau -> -tau (mirror X) and swapping the sign of both loads give the same sigma_eq."""
    r = TABLE[(TABLE.material == "hard") & (TABLE.test == 6)].iloc[0]
    p = params(r)
    base = carspa_critical_plane(*harmonics(255.1, 127.5, beta), p)["sigma_eq"]
    mirrored = carspa_critical_plane(*harmonics(255.1, -127.5, beta), p)["sigma_eq"]
    shifted = carspa_critical_plane(*harmonics(-255.1, -127.5, beta), p)["sigma_eq"]
    assert mirrored == pytest.approx(base, rel=1e-9) and shifted == pytest.approx(base, rel=1e-6)


def test_brittle_ductile_threshold_closed_forms():
    """tau_af/sigma_af = 1/sqrt(3): delta = 45 deg, Eq. (28); tension and torsion at their limits."""
    s_af = 300.0
    p = {"Sigm1": s_af, "Taum1": s_af / np.sqrt(3.0)}
    assert np.degrees(carspa_delta(**{"f_m1": s_af, "t_m1": s_af / np.sqrt(3.0)})) == pytest.approx(
        45.0
    )
    tension = carspa_critical_plane(*harmonics(s_af, 0.0, 0), p)
    torsion = carspa_critical_plane(*harmonics(0.0, s_af / np.sqrt(3.0), 0), p)
    assert tension["sigma_eq"] == pytest.approx(s_af) and torsion["sigma_eq"] == pytest.approx(s_af)


def test_extremely_hard_limit():
    """tau_af = sigma_af: delta = 0, the critical plane is normal to 1^ (pure tension: sigma_eq = sigma_a)."""
    p = {"Sigm1": 100.0, "Taum1": 100.0}
    res = carspa_critical_plane(*harmonics(80.0, 0.0, 0), p)
    assert res["delta"] == pytest.approx(0.0) and res["sigma_eq"] == pytest.approx(80.0)


def test_mean_normal_stress_enters_n_max():
    p = {"Sigm1": 300.0, "Taum1": 300.0 / np.sqrt(3.0)}
    S_mean, S_sin, S_cos = harmonics(200.0, 0.0, 0)
    S_mean[1, 1] = 50.0
    res = carspa_critical_plane(S_mean, S_sin, S_cos, p)
    assert res["N_m"] == pytest.approx(25.0) and res["N_a"] == pytest.approx(100.0)


def test_scope_errors():
    with pytest.raises(MaterialScopeError):
        carspa_critical_plane(*harmonics(100.0, 0.0, 0), {"Sigm1": 300.0, "Taum1": 160.0})
    S_mean, S_sin, S_cos = harmonics(100.0, 0.0, 0)
    S_mean[1, 1] = -150.0
    with pytest.raises(LoadingScopeError):
        carspa_critical_plane(S_mean, S_sin, S_cos, {"Sigm1": 300.0, "Taum1": 200.0})


def test_m_sigma_only_for_rotating_directions():
    """In-phase loading runs without m_sigma; out-of-phase needs m_sigma or b_fat."""
    p = {"Sigm1": 313.9, "Taum1": 196.2}
    carspa_critical_plane(*harmonics(255.1, 127.5, 0), p)
    with pytest.raises(MissingParamsError):
        carspa_critical_plane(*harmonics(255.1, 127.5, 60), p)
    via_b = carspa_critical_plane(*harmonics(255.1, 127.5, 60), {**p, "b_fat": -1 / 8.7})
    direct = carspa_critical_plane(*harmonics(255.1, 127.5, 60), {**p, "m_sigma": 8.7})
    assert via_b["sigma_eq"] == pytest.approx(direct["sigma_eq"])


def test_below_threshold_load_uses_the_peak_direction():
    """sigma_1 < c sigma_af everywhere: sigma_eq is homogeneous below the window (repository choice)."""
    p = {"Sigm1": 313.9, "Taum1": 196.2, "m_sigma": 8.7}
    lo = carspa_critical_plane(*harmonics(25.51, 12.75, 60), p)["sigma_eq"]
    lower = carspa_critical_plane(*harmonics(12.755, 6.375, 60), p)["sigma_eq"]
    assert lo == pytest.approx(2.0 * lower, rel=1e-9)


def test_registered_case_model():
    (model,) = get_models(["CARSPA"])
    assert isinstance(model, CarSpa)
    r = TABLE[(TABLE.material == "hard") & (TABLE.test == 12)].iloc[0]
    S_mean, S_sin, S_cos = harmonics(r.sigma_a, r.tau_a, r.beta)
    out = model.evaluate_case(
        S0=None, E0=None, S1=None, E1=None, params=params(r),
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos},
    )  # fmt: skip
    assert 100 * (out.values[""] - r.sigma_af) / r.sigma_af == pytest.approx(r.I_present, abs=0.5)
