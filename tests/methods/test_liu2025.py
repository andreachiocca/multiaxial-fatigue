"""LIU2025: Liu, Xie & Hao (2025), J Mater Eng Perform 34:24315-24322, doi:10.1007/s11665-025-11163-z,
checked against the PDF supplied by the maintainer.

Eqs. (5)-(13), (21)-(27). The life results are figures and band fractions only (Table 2) and the
experimental data are in the cited sources, so the tests check the stated identities and Fig. 3:
the shear and normal stress amplitudes on the critical plane against the phase angle, for
lambda = sigma_a/tau_a = sqrt(3) and 0.866. As Fig. 2 shows, the amplitudes of each phase are
scaled so that the largest von Mises stress of the cycle is 160 MPa (a circle of radius 160 at
90 deg for lambda = sqrt(3)). Digitised values, +-2 MPa:
- lambda = 0.866: both amplitudes at every phase angle;
- lambda = sqrt(3): the shear amplitude at every phase angle and the normal amplitude at 0 and
  90 deg; at 30 and 60 deg the equations give 69.0 and 102.7 MPa against about 75 and 108 read from
  the figure (unresolved; the same values follow on planes normal to the surface only).
"""

import numpy as np
import pytest

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.liu2025 import Liu2025, gerber_amplitude, serensen_amplitude
from multifat.planes.sampling import build_R_list

R = build_R_list()
Z = np.zeros((3, 3))
# Fig. 3 digitised: (lambda, phi) -> (shear amplitude, normal amplitude), MPa
FIG3 = {
    ("sqrt3", 0): (86.0, 56.0),
    ("sqrt3", 30): (86.0, 75.0),
    ("sqrt3", 45): (86.0, None),
    ("sqrt3", 60): (86.5, 108.0),
    ("sqrt3", 90): (92.0, 160.0),
    ("0.866", 0): (91.0, 36.0),
    ("0.866", 30): (90.5, 47.5),
    ("0.866", 45): (91.0, 57.5),
    ("0.866", 60): (91.0, 67.0),
    ("0.866", 90): (92.0, 80.0),
}
LAMBDA = {"sqrt3": np.sqrt(3.0), "0.866": 0.866}


class _Part(Liu2025):
    """Plane quantity of the critical plane: same plane rule, the requested amplitude as damage."""

    def __init__(self, part):
        self.part = part

    def evaluate_on_planes(self, **kw):
        S0r, S1r = kw["S0r"], kw["S1r"]
        dp, tau_a = super().evaluate_on_planes(**kw)
        sigma_a = 0.5 * np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
        return (tau_a if self.part == "shear" else sigma_a), tau_a + 1e-12 * dp


def _fig3_load(lam, phi):
    p = np.radians(phi)
    wt = np.linspace(0.0, 2.0 * np.pi, 20001)
    tau = 160.0 / np.max(np.sqrt((lam * np.sin(wt)) ** 2 + 3.0 * np.sin(wt - p) ** 2))

    def t(x, y):
        return np.array([[x, y, 0.0], [y, 0.0, 0.0], [0.0, 0.0, 0.0]])

    return t(lam * tau, tau * np.cos(p)), t(0.0, -tau * np.sin(p))


def _plane_value(lam, phi, part):
    S_sin, S_cos = _fig3_load(lam, phi)
    return best_cp_for_case(
        S0=S_sin, S1=-S_sin, E0=Z, E1=Z, R_list=R, method=_Part(part),
        params={"Su": 1000.0, "Tu": 700.0}, S_sin=S_sin, S_cos=S_cos, E_sin=Z, E_cos=Z,
    ).by_metric  # fmt: skip


@pytest.mark.parametrize("key", sorted(FIG3), ids=lambda k: f"{k[0]}-{k[1]}")
def test_fig3_shear_amplitude(key):
    assert _plane_value(LAMBDA[key[0]], key[1], "shear") == pytest.approx(FIG3[key][0], abs=2.0)


def _normal_rows():
    rows = []
    for key in sorted(FIG3):
        if FIG3[key][1] is None:
            continue
        if key in (("sqrt3", 30), ("sqrt3", 60)):
            reason = "equations give 69.0 / 102.7 MPa against about 75 / 108 read in Fig. 3"
            rows.append(pytest.param(key, marks=pytest.mark.xfail(strict=True, reason=reason)))
        else:
            rows.append(key)
    return rows


@pytest.mark.parametrize("key", _normal_rows(), ids=lambda k: f"{k[0]}-{k[1]}")
def test_fig3_normal_amplitude(key):
    assert _plane_value(LAMBDA[key[0]], key[1], "normal") == pytest.approx(FIG3[key][1], abs=2.0)


def test_eq26_branches_join_at_r0_and_eq24():
    """At R = 0 both branches give sqrt(2) sigma_a, i.e. S_0 = sqrt(2) S_-1 (Eq. 24)."""
    su, sa = 900.0, 200.0
    below = serensen_amplitude(sa, sa * (1 - 1e-12), su)
    above = serensen_amplitude(sa, sa * (1 + 1e-12), su)
    assert below == pytest.approx(np.sqrt(2) * sa, rel=1e-9)
    assert above == pytest.approx(np.sqrt(2) * sa, rel=1e-9)
    assert serensen_amplitude(sa, 0.0, su) == pytest.approx(sa)
    # compressive mean: linear branch, clipped at zero
    assert serensen_amplitude(sa, -100.0, su) == pytest.approx(sa - (np.sqrt(2) - 1) * 100.0)
    assert serensen_amplitude(sa, -1e4, su) == 0.0


def test_eq25_gerber():
    assert gerber_amplitude(100.0, 0.0, 600.0) == pytest.approx(100.0)
    assert gerber_amplitude(100.0, 300.0, 600.0) == pytest.approx(100.0 / 0.75)
    assert np.isinf(gerber_amplitude(100.0, 600.0, 600.0))


def test_reversed_tension_and_torsion_are_von_mises():
    """No mean stress: 45 deg plane gives sqrt(s^2/4 + 3 s^2/4) = s; torsion gives sqrt(3) tau."""
    model = Liu2025()

    def run(S):
        return best_cp_for_case(
            S0=S, S1=-S, E0=Z, E1=Z, R_list=R, method=model, params={"Su": 900.0, "Tu": 600.0},
        ).by_metric  # fmt: skip

    assert run(np.diag([300.0, 0.0, 0.0])) == pytest.approx(300.0, rel=1e-9)
    T = np.zeros((3, 3))
    T[0, 1] = T[1, 0] = 200.0
    assert run(T) == pytest.approx(np.sqrt(3) * 200.0, rel=1e-9)


def test_tu_estimate_warns():
    S = np.zeros((3, 3))
    S[0, 2] = S[2, 0] = 50.0
    with pytest.warns(RuntimeWarning, match="Tu"):
        Liu2025().evaluate_on_plane(S0r=S, S1r=-S, E0r=Z, E1r=Z, params={"Su": 900.0})


def test_mean_shear_beyond_tau_b_is_out_of_scope():
    from multifat.models.errors import LoadingScopeError

    S0, S1 = np.zeros((3, 3)), np.zeros((3, 3))
    S0[0, 2] = S0[2, 0] = 700.0
    S1[0, 2] = S1[2, 0] = 500.0  # mean shear 600 MPa on the plane
    with pytest.raises(LoadingScopeError, match="Gerber"):
        Liu2025().evaluate_on_plane(S0r=S0, S1r=S1, E0r=Z, E1r=Z, params={"Su": 900.0, "Tu": 550.0})


def test_minimum_normal_stress_beyond_sigma_b_is_out_of_scope():
    from multifat.models.errors import LoadingScopeError

    S0, S1 = np.diag([0.0, 0.0, 1300.0]), np.diag([0.0, 0.0, 1000.0])  # sigma_min = 1000 > 900
    with pytest.raises(LoadingScopeError, match="Serensen"):
        Liu2025().evaluate_on_plane(S0r=S0, S1r=S1, E0r=Z, E1r=Z, params={"Su": 900.0, "Tu": 600.0})


def test_fatigue_index_load_search_survives_trial_loads_beyond_sigma_b():
    """30NCD16 row 24 (sigma_m = 600, tau_a = 350 MPa): FI_load is finite (it stopped before)."""
    from multifat.eval.endurance import fatigue_index, reference_cases

    S_mean = np.diag([0.0, 600.0, 0.0])
    S_sin = np.zeros((3, 3))
    S_sin[0, 1] = S_sin[1, 0] = 350.0
    case = {
        "S0": S_mean + S_sin, "S1": S_mean - S_sin, "E0": Z, "E1": Z,
        "H": {"S_mean": S_mean, "S_sin": S_sin, "S_cos": Z},
    }  # fmt: skip
    params = {"Su": 1160.0, "Sigm1": 685.0, "Taum1": 430.0}
    refs = reference_cases(685.0, 430.0, None, None)
    (row,) = fatigue_index(Liu2025Adapter(), params, [case], refs, build_R_list(n1=10, n2=20))
    assert np.isfinite(row["FI_load"]) and np.isfinite(row["FI_axial"])


def Liu2025Adapter():
    from multifat.models.registry import get_models

    return get_models(["LIU2025"])[0]
