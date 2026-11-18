"""Probabilistic volume fatigue model of Flacelière and Morel (2004).

Published in *Fatigue & Fracture of Engineering Materials & Structures* 27 (2004),
pp. 1123-1135, doi:10.1111/j.1460-2695.2004.00823.x (roadmap item G24), Eqs. (3), (12)-(14)::

    P_F = 1 - exp(-(V / (V0 * sigma_u^m)) * sigma_eq,max^m * H_m)
    H_m = (1 / (V * sigma_eq,max^m)) * int_V sigma_eq^m dV
    sigma_eq = M_sigma + p * H_max

where ``M_sigma`` is the RMS resolved shear amplitude (``t_-1`` torsion, ``s_-1/sqrt(3)``
tension) and ``H_m`` takes closed-form specimen factors (Eq. 14: ``1`` tension,
``1/(m+1)`` 4-pt plane, ``2/(m+2)`` 4-pt rotating, ``1/(m+1)^2`` 3-pt plane,
``2/((m+1)(m+2))`` 3-pt rotating).

Repository choices:
- Supplied as point kernels (``heterogeneity_factor``, ``failure_probability``,
  ``median_strength``) rather than a registered model because workbooks lack stressed
  specimen volumes ``V`` and the paper omits ``V`` and ``V0 * sigma_u^m`` for its tests.
- Uses exact analytical heterogeneity factors from Eq. (14) directly instead of integrating
  a stress field numerically.
"""

from __future__ import annotations

import numpy as np

_H_FACTORS = {
    "tension": lambda m: 1.0,
    "four_point_plane": lambda m: 1.0 / (m + 1.0),
    "four_point_rotating": lambda m: 2.0 / (m + 2.0),
    "three_point_plane": lambda m: 1.0 / (m + 1.0) ** 2,
    "three_point_rotating": lambda m: 2.0 / ((m + 1.0) * (m + 2.0)),
}


def heterogeneity_factor(case: str, m: float) -> float:
    """H_m of Eq. (14)."""
    return float(_H_FACTORS[case](m))


def failure_probability(sigma_eq_max: float, m: float, volume: float, v0_su_m: float, h_m: float):
    """P_F of Eq. (13)."""
    return float(1.0 - np.exp(-(volume / v0_su_m) * sigma_eq_max**m * h_m))


def median_strength(m: float, volume: float, v0_su_m: float, h_m: float) -> float:
    """sigma_eq,max evaluated at P_F = 50 %."""
    return float((np.log(2.0) * v0_su_m / (volume * h_m)) ** (1.0 / m))
