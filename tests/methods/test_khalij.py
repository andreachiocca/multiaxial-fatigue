"""KHALIJ (N19): Khalij, Pagnacco & Troian (2017), Int J Fatigue 99:137-150,
doi:10.1016/j.ijfatigue.2017.02.016, checked against the PDF supplied by the maintainer.

Eqs. (12), (22), (23), (26), pp.138-140. Table 6 (p.148) prints the "Damage" of the zero-mean
criterion for 53 in-phase fatigue limits of eight materials (Table 1, p.141).

Findings:
- the OP column is the left-hand side of Eq. (23) with alpha = 1/(3 r^2) and r = t_-1/s_-1 from the
  limits of Table 1: 51 of 53 rows to print precision. With the two-digit r printed in Table 6,
  7 rows miss, so the table was computed with the unrounded ratio;
- two 3%C cast-iron rows (x = 0.99, t = 0.22 and x = 0.59, t = 0.75) give 1.035 and 1.063 against
  the printed 1.06 and 1.01, with either r (strict xfails, unexplained);
- the mean-stress extension, Eqs. (29)-(31), is not implemented: the printed gamma_2 term has
  units of 1/MPa and no reading reproduces Table 7 (see the module docstring).
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.errors import LoadingScopeError
from multifat.models.invariant_methods.khalij import Khalij, khalij_alpha, khalij_damage

TABLE = load_table("khalij2017_table6")
MISPRINTED = {("3%C cast iron", 0.99, 0.22), ("3%C cast iron", 0.59, 0.75)}


def _plane(sxx, syy, sxy):
    return np.array([[sxx, sxy, 0.0], [sxy, syy, 0.0], [0.0, 0.0, 0.0]])


def eq23(x, y, t, s_m1, t_m1, alpha, cos_phi=1.0):
    """Left-hand side of Eq. (23) as printed, in normalised amplitudes."""
    i1 = np.sqrt(max(0.0, x * x + y * y + 2 * x * y * cos_phi))
    return alpha * (x * x + y * y - x * y * cos_phi) + t * t + (1 - alpha) * i1


def _rows():
    reason = (
        "Table 6 prints a Damage that Eq. (23) does not give for this row (0.025 and 0.053 off)"
    )
    return [
        (
            pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=reason))
            if (r.material, r.x, r.t) in MISPRINTED
            else r
        )
        for r in TABLE.itertuples()
    ]


@pytest.mark.parametrize("r", _rows(), ids=lambda r: f"{r.material}-{r.x}-{r.y}-{r.t}")
def test_table6_op_column(r):
    S = _plane(r.x * r.s_m1, r.y * r.s_m1, r.t * r.t_m1)
    # the inputs are printed to two digits, so the printed Damage is matched to one unit of its last digit
    assert abs(round(khalij_damage(S, np.zeros((3, 3)), r.s_m1, r.t_m1), 2) - r.D_OP) <= 0.01 + 1e-9


def test_table6_needs_the_unrounded_ratio():
    """With the printed two-digit r, seven OP rows miss (the table used t_-1/s_-1 of Table 1)."""
    misses = sum(
        abs(round(eq23(r.x, r.y, r.t, r.s_m1, r.t_m1, 1 / (3 * r.r_printed**2)), 2) - r.D_OP)
        > 0.01 + 1e-9
        for r in TABLE.itertuples()
    )
    assert misses == 7


@pytest.mark.parametrize("phi", [0.0, 60.0, 90.0, 180.0])
def test_invariant_form_equals_eq23(phi):
    """J_2,a/t_-1^2 + (1 - alpha) I_1,a/s_-1 is Eq. (23) with phase phi_yy and a phase-shifted shear."""
    s_m1, t_m1, sxx, syy, sxy = 343.0, 204.0, 150.0, 120.0, 80.0
    p, q = np.radians(phi), np.radians(37.0)  # the shear phase must not matter (Eq. 12)
    S_sin = _plane(sxx, syy * np.cos(p), sxy * np.cos(q))
    S_cos = _plane(0.0, -syy * np.sin(p), -sxy * np.sin(q))
    expected = eq23(
        sxx / s_m1, syy / s_m1, sxy / t_m1, s_m1, t_m1, khalij_alpha(s_m1, t_m1), np.cos(p)
    )
    assert khalij_damage(S_sin, S_cos, s_m1, t_m1) == pytest.approx(expected, rel=1e-12)


def test_limits_and_frame_invariance():
    s_m1, t_m1 = 361.0, 228.0
    zero = np.zeros((3, 3))
    assert khalij_damage(_plane(s_m1, 0, 0), zero, s_m1, t_m1) == pytest.approx(1.0)
    assert khalij_damage(_plane(0, 0, t_m1), zero, s_m1, t_m1) == pytest.approx(1.0)
    # pure shear written in its principal frame (Eq. 24-26: OP is the isotropic setting)
    assert khalij_damage(_plane(t_m1, -t_m1, 0), zero, s_m1, t_m1) == pytest.approx(1.0)
    c, s = np.cos(0.4), np.sin(0.4)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    S = _plane(200.0, 90.0, 60.0)
    assert khalij_damage(R @ S @ R.T, zero, s_m1, t_m1) == pytest.approx(
        khalij_damage(S, zero, s_m1, t_m1), rel=1e-12
    )


def test_scope_mean_and_out_of_plane():
    params = {"Sigm1": 343.0, "Taum1": 204.0}
    zero = np.zeros((3, 3))
    harmonics = {"S_mean": _plane(50, 0, 0), "S_sin": _plane(200, 0, 0), "S_cos": zero}
    with pytest.raises(LoadingScopeError, match="mean-stress"):
        Khalij().evaluate_case(
            S0=None, E0=None, S1=None, E1=None, params=params, harmonics=harmonics
        )
    S_sin = _plane(200, 0, 0)
    S_sin[2, 2] = 30.0
    harmonics = {"S_mean": zero, "S_sin": S_sin, "S_cos": zero}
    with pytest.raises(LoadingScopeError, match="plane stress"):
        Khalij().evaluate_case(
            S0=None, E0=None, S1=None, E1=None, params=params, harmonics=harmonics
        )
