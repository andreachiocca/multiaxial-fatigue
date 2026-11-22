"""Bolchoun, Baumgartner & Kaufmann (2017), doi:10.1016/j.ijfatigue.2017.01.012: Eqs. (9)-(15).

The paper gives no table of computed values for a reproducible stress state (its results are
figures for notched welded joints), so the tests check what the equations fix in closed form:
the hemisphere integral of Eq. (12), the out-of-phase measure of Eq. (9), the end points of the
interpolation factor and of Eq. (11), and the calibration identity of Eqs. (14)-(15).
"""

from __future__ import annotations

import numpy as np
import pytest

from multifat.eval.material_fits import out_of_phase_fill
from multifat.models.errors import LoadingScopeError, MissingParamsError
from multifat.models.integral_methods.bolchoun import (
    REQUIRED,
    Bolchoun,
    WoehlerSet,
    bolchoun_quantities,
    interpolation_factor,
    out_of_phase_measure,
    substitute_amplitude,
    tension_torsion_path,
    unit_factors,
)

ZERO = np.zeros((3, 3))
# Axial and torsion Basquin curves of a steel, and an out-of-phase curve above the in-phase one
PARAMS = {
    "sigma_f": 1154.0,
    "b_fat": -0.061,
    "tau_f": 940.0,
    "b0": -0.065,
    "sigma_f_OP": 773.0,
    "b_OP": -0.036,
    "lambda_OP": 0.5773,
    "delta_OP": 90.0,
}


def _axial(amplitude):
    S = np.zeros((3, 3))
    S[1, 1] = amplitude
    return S


def _shear(amplitude):
    S = np.zeros((3, 3))
    S[0, 1] = S[1, 0] = amplitude
    return S


def test_eq_12_unit_tension_and_torsion():
    """sigma_n,a = sigma |cos^2 a| and tau |2 n_x n_y|: hemisphere integrals 2 pi / 3 and 8 / 3."""
    c_ax, c_tor = unit_factors()
    assert c_ax == pytest.approx(2.0 * np.pi / 3.0, rel=1e-9)
    # |2 n_x n_y| has kinks on the grid lines; the same grid serves every case, so the ratio is
    # consistent across the model
    assert c_tor == pytest.approx(8.0 / 3.0, rel=2e-3)


def test_eq_12_is_linear_in_the_amplitude():
    S_sin, S_cos = tension_torsion_path(300.0, 0.6, 90.0)
    assert substitute_amplitude(2 * S_sin, 2 * S_cos) == pytest.approx(
        2 * substitute_amplitude(S_sin, S_cos), rel=1e-12
    )


def test_eq_9_vanishes_for_proportional_loading():
    S = _axial(200.0) + _shear(120.0)
    assert out_of_phase_measure(S, ZERO) == 0.0
    assert out_of_phase_measure(S, 0.5 * S) == pytest.approx(0.0, abs=1e-9)


def test_eq_9_ellipse_area_matches_the_time_integral():
    """1/2 int |tau_0 x dtau_0/dt| dt over one period is pi |tau_s x tau_c| on every plane."""
    S_sin, S_cos = tension_torsion_path(250.0, 0.7, 60.0)
    rng = np.random.default_rng(3)
    t = np.linspace(0.0, 2.0 * np.pi, 20001)
    for n in rng.normal(size=(5, 3)):
        n /= np.linalg.norm(n)

        def shear(S):
            v = S @ n
            return v - (n @ v) * n

        tau_s, tau_c = shear(S_sin), shear(S_cos)
        tau = np.outer(np.sin(t), tau_s) + np.outer(np.cos(t), tau_c)
        rate = np.outer(np.cos(t), tau_s) - np.outer(np.sin(t), tau_c)
        integrand = 0.5 * np.linalg.norm(np.cross(tau, rate), axis=1)
        swept = np.trapezoid(integrand, t)
        assert swept == pytest.approx(np.pi * np.linalg.norm(np.cross(tau_s, tau_c)), rel=1e-9)


def test_eq_9_grows_with_the_square_of_the_amplitude():
    S_sin, S_cos = tension_torsion_path(100.0, 0.5, 90.0)
    assert out_of_phase_measure(3 * S_sin, 3 * S_cos) == pytest.approx(
        9 * out_of_phase_measure(S_sin, S_cos), rel=1e-12
    )


def test_interpolation_factor_end_points():
    """Section 3: tau_max,a = sigma_max,a / 2 in tension, tau_max,a = sigma_max,a in torsion."""
    assert interpolation_factor(_axial(1.0), ZERO) == pytest.approx(0.5, rel=1e-9)
    assert interpolation_factor(_shear(1.0), ZERO) == pytest.approx(1.0, rel=1e-9)


def test_interpolation_factor_proportional_closed_form():
    sigma, tau = 100.0, 60.0
    radius = np.hypot(sigma / 2.0, tau)
    f = interpolation_factor(_axial(sigma) + _shear(tau), ZERO)
    assert f == pytest.approx(radius / (sigma / 2.0 + radius), rel=1e-9)


def test_eq_11_end_points_are_the_axial_and_torsion_curves():
    curves = WoehlerSet(PARAMS)
    c_ax, c_tor = unit_factors()
    k_ax, l_ax = curves.interpolated(0.5)
    k_tor, l_tor = curves.interpolated(1.0)
    assert k_ax == pytest.approx(-1.0 / PARAMS["b_fat"])
    assert k_tor == pytest.approx(-1.0 / PARAMS["b0"])
    assert l_ax == pytest.approx(c_ax * PARAMS["sigma_f"] * 2e6 ** PARAMS["b_fat"])
    assert l_tor == pytest.approx(c_tor * PARAMS["tau_f"] * 2e6 ** PARAMS["b0"])
    # f_I outside [1/2, 1] uses the nearest curve (Section 3)
    assert curves.interpolated(0.3) == curves.interpolated(0.5)
    assert curves.interpolated(1.2) == curves.interpolated(1.0)


@pytest.mark.parametrize("amplitude", [450.0, 600.0, 800.0])
def test_uniaxial_tension_and_torsion_return_their_basquin_curves(amplitude):
    model, curves = Bolchoun(), WoehlerSet(PARAMS)
    for S, coef, exp in (
        (_axial(amplitude), PARAMS["sigma_f"], PARAMS["b_fat"]),
        (_shear(amplitude), PARAMS["tau_f"], PARAMS["b0"]),
    ):
        q = bolchoun_quantities(ZERO, S, ZERO)
        assert q["F_OP"] == 0.0
        aux = {"dp": q["sigma_as"], "F_OP": q["F_OP"], "f_I": q["f_I"]}
        life = model.native_life(aux, PARAMS)
        assert life == pytest.approx(0.5 * (amplitude / coef) ** (1.0 / exp), rel=1e-6)
        assert life == pytest.approx(curves.life(q["sigma_as"], 0.0, q["f_I"]), rel=1e-12)


@pytest.mark.parametrize(
    "op_curve, sign",
    [((773.0, -0.036), -1.0), ((900.0, -0.08), 1.0)],
    ids=["benefit (minus variant)", "extra damage"],
)
def test_eqs_14_15_reproduce_the_out_of_phase_curve_on_its_own_path(op_curve, sign):
    """On the calibration path k_I,C = k_I,OP, and Eq. (15) returns D_total of Eq. (14)."""
    params = {**PARAMS, "sigma_f_OP": op_curve[0], "b_OP": op_curve[1]}
    curves = WoehlerSet(params)
    for sigma_a in (480.0, 520.0, 600.0):
        S_sin, S_cos = tension_torsion_path(sigma_a, params["lambda_OP"], params["delta_OP"])
        q = bolchoun_quantities(ZERO, S_sin, S_cos)
        assert curves.out_of_phase_damage(q["F_OP"])[1] == sign
        expected = 0.5 * (sigma_a / op_curve[0]) ** (1.0 / op_curve[1])
        assert curves.life(q["sigma_as"], q["F_OP"], q["f_I"]) == pytest.approx(expected, rel=1e-9)


def test_mean_stress_is_out_of_scope():
    with pytest.raises(LoadingScopeError):
        bolchoun_quantities(_axial(50.0), _axial(200.0), ZERO)


def test_missing_out_of_phase_curve_skips_the_model():
    with pytest.raises(MissingParamsError):
        Bolchoun().native_life({"dp": 500.0}, {k: PARAMS[k] for k in REQUIRED[:4]})


def _case(sigma_a, lambda_OP, delta_OP, life, S_mean=None):
    S_sin, S_cos = tension_torsion_path(sigma_a, lambda_OP, delta_OP)
    mean = ZERO if S_mean is None else S_mean
    return {
        "harmonics": {"S_mean": mean, "S_sin": S_sin, "S_cos": S_cos},
        "meta": {"Nf_cycles": life},
    }


def test_out_of_phase_fill_uses_the_largest_path_and_its_true_ratio():
    lives = np.array([1e4, 3e4, 1e5, 3e5, 1e6])
    amps = 800.0 * (2 * lives) ** -0.05
    cases = {}
    for i, (a, n) in enumerate(zip(amps, lives)):
        # tau_a / sigma_a scatters around 0.5774 inside one rounded group
        cases[i] = _case(a, 0.5774 + 0.0002 * (i - 2), 90.0, n)
    cases[10] = _case(500.0, 1.0, 45.0, 2e5)  # a smaller out-of-phase path
    cases[11] = _case(500.0, 0.6, 0.0, 2e5)  # proportional: not an out-of-phase test
    cases[12] = _case(500.0, 0.6, 90.0, 2e5, S_mean=_axial(100.0))  # mean stress: excluded
    out = out_of_phase_fill({"M": {"cases": cases}}, "M")
    p = out["params"]
    assert sorted(out["design_case_ids"]) == [0, 1, 2, 3, 4]
    assert p["sigma_f_OP"] == pytest.approx(800.0) and p["b_OP"] == pytest.approx(-0.05)
    assert p["lambda_OP"] == pytest.approx(0.5774, rel=1e-12)
    assert p["delta_OP"] == pytest.approx(90.0)


def test_out_of_phase_fill_without_a_path():
    out = out_of_phase_fill({"M": {"cases": {0: _case(500.0, 0.6, 0.0, 1e5)}}}, "M")
    assert out["params"] == {} and out["design_case_ids"] == []
