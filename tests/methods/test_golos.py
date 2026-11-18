"""GOLOS (N13): Golos, Debski & Debski (2014), Theor Appl Fract Mech 73:3-8,
doi:10.1016/j.tafmec.2014.07.005, checked against the PDF supplied by the maintainer.

Eqs. (26)-(32), pp.5-7. The paper reports its validation as figures (Figs. 2-7) and as the scatter
index of Table 2, computed from normal distances to the curves (Fig. 8) for data points that are
not tabulated; Table 2 is therefore not reproduced. The tests check the stated identities:
- A + B = 1 (bending limit on the curve) for Eqs. (30) and (31);
- Eq. (30) gives A = 1, B = 0 at kappa = sqrt(3) (Gough ellipse, Eq. 26) and A = 1/3, B = 2/3 at
  kappa = 4/3, where Eq. (29) coincides with Debski's equivalent stress, Eqs. (27)-(28);
- Eq. (31) approximates Eq. (30) within 0.01 over the paper's range 1.5 <= kappa <= 1.75;
- finding: the text also states A = 1, B = 0 for kappa = 2, which Eq. (30) does not give (1.45).
"""

import numpy as np
import pytest
from _paper import bending_torsion, frame_with_normal, rotate

from multifat.models.errors import LoadingScopeError
from multifat.models.invariant_methods.golos import (
    Golos,
    bending_torsion_amplitudes,
    golos_coefficients,
    golos_index,
)

Z = np.zeros((3, 3))


@pytest.mark.parametrize("kappa", [4 / 3, 1.53, 1.6, np.sqrt(3), 1.75, 2.0])
@pytest.mark.parametrize("exact", [True, False])
def test_bending_limit_lies_on_the_curve(kappa, exact):
    A, B = golos_coefficients(kappa * 200.0, 200.0, exact=exact)
    assert A + B == pytest.approx(1.0, abs=1e-12)


def test_eq30_special_cases():
    A, B = golos_coefficients(np.sqrt(3) * 100, 100, exact=True)
    assert (A, B) == pytest.approx((1.0, 0.0), abs=1e-12)
    A, B = golos_coefficients(4 / 3 * 100, 100, exact=True)
    assert (A, B) == pytest.approx((1 / 3, 2 / 3), abs=1e-12)
    A, _ = golos_coefficients(200, 100, exact=True)
    assert A == pytest.approx(1.4480, abs=1e-4)  # not the A = 1 stated in the text for kappa = 2


def test_eq31_approximates_eq30_on_the_paper_range():
    for kappa in np.linspace(1.5, 1.75, 26):
        exact = np.array(golos_coefficients(kappa, 1.0, exact=True))
        approx = np.array(golos_coefficients(kappa, 1.0))
        assert np.max(np.abs(exact - approx)) < 0.01


def test_kappa_4_3_is_debski_equivalent_stress_eq27():
    """On Eq. (28), sigma/3 + 2/3 sqrt(sigma^2 + 4 tau^2) = f_-1 (Eq. 27 with f/t = 4/3)."""
    f, t = 400.0, 300.0
    for s in np.linspace(0.0, f, 9):
        tau = t * np.sqrt(1 - (s / f) ** 2 / 3 - 2 / 3 * s / f)  # Eq. (28) solved for tau
        assert s / 3 + 2 / 3 * np.sqrt(s**2 + 4 * tau**2) == pytest.approx(f, rel=1e-12)


def test_limits_and_gough_ellipse():
    f, t = 313.9, 196.2
    assert golos_index(f, 0.0, f, t) == pytest.approx(1.0, abs=1e-12)
    assert golos_index(0.0, t, f, t) == pytest.approx(1.0, abs=1e-12)
    # kappa = sqrt(3) with Eq. (31): close to the Gough ellipse
    t = f / np.sqrt(3)
    assert golos_index(f / np.sqrt(2), t / np.sqrt(2), f, t) == pytest.approx(1.0, abs=0.01)


def test_rotated_bending_torsion_state_is_recognised():
    S = bending_torsion(250.0, 120.0)
    R = frame_with_normal([0.3, -0.5, 0.8])
    assert bending_torsion_amplitudes(Z, rotate(S, R), Z) == pytest.approx((250.0, 120.0))
    # the cosine harmonic alone, and a 180 deg (anti-phase) sign change, are still in phase
    assert bending_torsion_amplitudes(Z, Z, -S) == pytest.approx((250.0, 120.0))


@pytest.mark.parametrize(
    "parts,match",
    [
        ((bending_torsion(50.0, 0.0), bending_torsion(200.0, 80.0), Z), "fully reversed"),
        ((Z, bending_torsion(200.0, 0.0), bending_torsion(0.0, 80.0)), "in-phase"),
        ((Z, np.diag([200.0, 100.0, 0.0]), Z), "equal sign"),
        ((Z, np.diag([200.0, -100.0, 50.0]), Z), "plane stress"),
    ],
)
def test_out_of_scope_cases_are_rejected(parts, match):
    with pytest.raises(LoadingScopeError, match=match):
        bending_torsion_amplitudes(*parts)


def test_registered_model_returns_the_index():
    S = bending_torsion(200.0, 100.0)
    out = Golos().evaluate_case(
        S0=S,
        S1=-S,
        E0=None,
        E1=None,
        params={"Sigm1": 313.9, "Taum1": 196.2},
        harmonics={"S_mean": Z, "S_sin": S, "S_cos": Z},
    )
    assert out.values[""] == pytest.approx(golos_index(200.0, 100.0, 313.9, 196.2))


def test_golos_needs_no_strain():
    """The stress-only probe (a 90 deg case) is out of scope for GOLOS, not a strain failure."""
    from multifat.models.registry import stress_only_models

    assert "GOLOS" in stress_only_models()


def test_load_scaling_treats_out_of_scope_trial_loads_as_infinite_damage():
    """Searches over scaled loads (FI_load, load error) step back instead of stopping."""
    from multifat.eval.endurance import scaled_damage

    class Bounded:
        name = "BOUNDED"

        def evaluate_case(self, *, S0, E0, S1, E1, params, R_list=None, harmonics=None):
            from multifat.models.base import CaseResult

            amp = float(np.abs(harmonics["S_sin"]).max())
            if amp > 150.0:
                raise LoadingScopeError("BOUNDED: beyond the domain")
            return CaseResult(values={"": amp})

    S = bending_torsion(100.0, 0.0)
    case = {"S0": S, "S1": -S, "E0": Z, "E1": Z, "H": {"S_mean": Z, "S_sin": S, "S_cos": Z}}
    assert scaled_damage(Bounded(), {}, case, 1.2, None) == pytest.approx(120.0)
    assert scaled_damage(Bounded(), {}, case, 2.0, None) == float("inf")


def test_load_factor_brackets_against_an_infinite_upper_end():
    from multifat.eval.endurance import load_factor

    class Bounded:
        name = "BOUNDED"

        def evaluate_case(self, *, S0, E0, S1, E1, params, R_list=None, harmonics=None):
            from multifat.models.base import CaseResult

            amp = float(np.abs(harmonics["S_sin"]).max())
            if amp > 150.0:
                raise LoadingScopeError("BOUNDED: beyond the domain")
            return CaseResult(values={"": amp})

    S = bending_torsion(100.0, 0.0)
    case = {"S0": S, "S1": -S, "E0": Z, "E1": Z, "H": {"S_mean": Z, "S_sin": S, "S_cos": Z}}
    assert load_factor(Bounded(), {}, case, 140.0, None) == pytest.approx(1.4, rel=1e-9)
