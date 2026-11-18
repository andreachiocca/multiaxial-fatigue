"""GAM (N08): Goncalves, Araujo & Mamiya (2005), Int J Fatigue 27:177-187,
doi:10.1016/j.ijfatigue.2004.05.006, checked against the PDF supplied by the maintainer.

Eqs. (11)-(13), p.181, and Appendix A, pp.185-186. Tables 1-4 (pp.182-183) give the error index
I = 100 (LHS - lambda)/lambda, Eq. (14), of the current model (column I^d) for 41 fatigue limits.

Findings:
- every in-phase row (beta = 0) is reproduced to print precision by Eq. (12) as written;
- the out-of-phase rows were computed with the in-phase peak principal stress
  s/2 + sqrt(s^2/4 + x^2) (s = sigma_a + sigma_m, x = tau_a + tau_m), not with Eq. (11), the maximum
  over the cycle of the largest principal stress. With that substitution all 41 rows match to
  0.02 points; with Eq. (11) the 27 out-of-phase rows are lower by 0.3 to 13.4 points. The model
  implements Eq. (11); those rows are strict xfails.
"""

import numpy as np
import pytest
from _paper import bending_torsion, load_table

from multifat.models.invariant_methods.goncalves_araujo_mamiya import (
    gam_constants,
    gam_index,
    prism_shear_amplitude,
)

TABLE = load_table("goncalves2005_tables1_4")


def _parts(r):
    b = np.radians(r.beta)
    return (
        bending_torsion(r.sigma_m, r.tau_m),
        bending_torsion(r.sigma_a, r.tau_a * np.cos(b)),
        bending_torsion(0.0, -r.tau_a * np.sin(b)),
    )


def _rows():
    reason = "Tables 1-4 use the in-phase peak principal stress instead of Eq. (11)"
    return [
        r if r.beta == 0 else pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=reason))
        for r in TABLE.itertuples()
    ]


@pytest.mark.parametrize("r", _rows(), ids=lambda r: r.row)
def test_tables_1_4_eq12(r):
    I = 100 * (gam_index(*_parts(r), r.f_m1, r.t_m1) - 1)
    assert I == pytest.approx(r.I_gam, abs=0.02)


@pytest.mark.parametrize("r", list(TABLE.itertuples()), ids=lambda r: r.row)
def test_tables_1_4_used_the_in_phase_peak(r):
    """Every printed I^d follows from Eq. (12) with the phase-free peak principal stress."""
    weight, limit = gam_constants(r.f_m1, r.t_m1)
    s, x = r.sigma_a + r.sigma_m, r.tau_a + r.tau_m
    peak = s / 2 + np.hypot(s / 2, x)
    _, S_sin, S_cos = _parts(r)
    I = 100 * ((prism_shear_amplitude(S_sin, S_cos) + weight * peak) / limit - 1)
    assert I == pytest.approx(r.I_gam, abs=0.02)


def test_eq13_constants_put_both_limits_on_the_threshold():
    """Eq. (13) makes the index 1 at f_-1 (bending) and t_-1 (torsion)."""
    f, t = 313.9, 196.2
    zero = np.zeros((3, 3))
    assert gam_index(zero, bending_torsion(f, 0.0), zero, f, t) == pytest.approx(1.0, abs=1e-12)
    assert gam_index(zero, bending_torsion(0.0, t), zero, f, t) == pytest.approx(1.0, abs=1e-12)


def test_prism_amplitude_is_basis_invariant_sum_of_a_i():
    """Eq. (A.17) in the basis of Eq. (A.3) equals the closed form, for an out-of-phase path."""
    rng = np.random.default_rng(3)
    A, B = (0.5 * (M + M.T) for M in rng.normal(size=(2, 3, 3)) * 100)
    times = np.linspace(0.0, 2 * np.pi, 20001)
    S = np.sin(times)[:, None, None] * A + np.cos(times)[:, None, None] * B
    dev = S - np.trace(S, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
    s = np.stack(
        [
            np.sqrt(1.5) * dev[:, 0, 0],
            (dev[:, 1, 1] - dev[:, 2, 2]) / np.sqrt(2),
            np.sqrt(2) * dev[:, 0, 1],
            np.sqrt(2) * dev[:, 0, 2],
            np.sqrt(2) * dev[:, 1, 2],
        ]
    )
    a = 0.5 * (s.max(axis=1) - s.min(axis=1))
    assert np.sqrt(np.sum(a**2)) == pytest.approx(prism_shear_amplitude(A, B), rel=1e-6)


def test_invalid_limits_are_rejected_and_soft_metals_warned():
    zero = np.zeros((3, 3))
    with pytest.raises(ValueError):
        gam_constants(200.0, 250.0)
    with pytest.warns(RuntimeWarning, match="hard-metal"):
        gam_index(zero, bending_torsion(100.0, 0.0), zero, 400.0, 200.0)
