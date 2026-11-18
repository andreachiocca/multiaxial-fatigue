"""Goncalves, Araujo & Mamiya (2005) multiaxial HCF endurance criterion for hard metals (N08).

C. A. Goncalves, J. A. Araujo, E. N. Mamiya, Int. J. Fatigue 27 (2005) 177-187,
doi:10.1016/j.ijfatigue.2004.05.006, Eqs. (11)-(13) and Appendix A::

    sqrt(sum_i a_i^2) + kappa sigma_p,max <= lambda                                   (12)
    a_i = (max_t s_i(t) - min_t s_i(t)) / 2,  i = 1..5                               (A.17)
    sigma_p,max = max_t max_i sigma_p,i(t)                                           (11)
    kappa  = sqrt(2) / (f_-1 - t_-1) (t_-1 - f_-1/sqrt(3))
    lambda = sqrt(2) t_-1 f_-1 / (f_-1 - t_-1) (1 - 1/sqrt(3))                       (13)

s_i(t) are the components of the deviatoric stress in an orthonormal basis of the deviatoric space
(Eq. A.3); a_i are the half-sides of the rectangular prism that encloses the deviatoric path. The
returned fatigue index is LHS / lambda (1 on the fatigue limit); the paper's error index is
I = 100 (index - 1), Eq. (14). The paper restricts the criterion to hard metals,
1.3 <= f_-1/t_-1 < sqrt(3).

For one harmonic s_i(t) = m_i + p_i sin t + q_i cos t, so a_i^2 = p_i^2 + q_i^2 and, the basis being
orthonormal, sum_i a_i^2 = ||dev S_sin||_F^2 + ||dev S_cos||_F^2 in every basis (Proposition 1:
the sum equals the squared semi-axes of the enclosing ellipse).

Repository choices:
  - sigma_p,max is the maximum over the cycle of the largest principal stress, as Eq. (11) states.
    Tables 1-4 of the paper were instead computed with the in-phase peak
    (s/2 + sqrt(s^2/4 + x^2), s = sigma_a + sigma_m, x = tau_a + tau_m) for every phase angle;
    their out-of-phase rows are therefore not reproduced (strict xfails in the paper test);
  - f_-1 <= t_-1 is rejected (kappa, lambda undefined); outside the hard-metal range the criterion
    is evaluated with a warning (f_-1/t_-1 >= sqrt(3) makes kappa <= 0);
  - the fatigue-limit index is used as the damage parameter of the shared master curve
    DP = A N^b for finite-life comparisons (not part of the criterion).

User guide: docs/method-guides/GAM.md.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from ..base import CaseResult
from .stress_invariants import harmonic_parts, peak_principal_stress

SQRT2, SQRT3 = np.sqrt(2.0), np.sqrt(3.0)


def _deviatoric(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    return T - np.trace(T) / 3.0 * np.eye(3)


def prism_shear_amplitude(S_sin, S_cos) -> float:
    """sqrt(sum a_i^2) of Eq. (A.18) for one harmonic: Frobenius norms of the deviatoric parts."""
    return float(np.sqrt(np.sum(_deviatoric(S_sin) ** 2) + np.sum(_deviatoric(S_cos) ** 2)))


def gam_constants(f_m1: float, t_m1: float) -> tuple[float, float]:
    """(kappa, lambda) of Eq. (13) from the fully reversed bending and torsion limits."""
    f_m1, t_m1 = float(f_m1), float(t_m1)
    if not (f_m1 > t_m1 > 0.0):
        raise ValueError(f"GAM needs f_-1 > t_-1 > 0, got f_-1 = {f_m1}, t_m1 = {t_m1}")
    weight = SQRT2 / (f_m1 - t_m1) * (t_m1 - f_m1 / SQRT3)
    limit = SQRT2 * t_m1 * f_m1 / (f_m1 - t_m1) * (1.0 - 1.0 / SQRT3)
    return float(weight), float(limit)


def gam_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    """Fatigue index (sqrt(sum a_i^2) + kappa sigma_p,max) / lambda of Eq. (12)."""
    weight, limit = gam_constants(f_m1, t_m1)
    kappa = float(f_m1) / float(t_m1)  # repository notation: f_-1/t_-1, not Eq. (12)'s kappa
    if not (1.3 <= kappa < SQRT3):
        warnings.warn(
            f"GAM: f_-1/t_-1 = {kappa:.3f} lies outside the hard-metal range [1.3, sqrt(3)) "
            "of the source",
            RuntimeWarning,
            stacklevel=2,
        )
    sigma_p_max = peak_principal_stress(
        np.asarray(S_mean, float), np.asarray(S_sin, float), np.asarray(S_cos, float)
    )
    return (prism_shear_amplitude(S_sin, S_cos) + weight * sigma_p_max) / limit


class GoncalvesAraujoMamiya:
    name = "GAM"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self,
        *,
        S0,
        E0,
        S1,
        E1,
        params: Mapping[str, float],
        R_list=None,
        harmonics=None,
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        index = gam_index(S_mean, S_sin, S_cos, float(params["Sigm1"]), float(params["Taum1"]))
        return CaseResult(values={"": index})
