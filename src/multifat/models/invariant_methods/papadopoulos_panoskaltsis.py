"""Gradient-sensitive multiaxial fatigue criterion of Papadopoulos & Panoskaltsis (1996).

Published in *Engineering Fracture Mechanics* 55 (1996), pp. 513-528,
doi:10.1016/S0013-7944(96)00047-1 (roadmap item G29), Eqs. (13)-(23)::

    sqrt(J2,a) + alpha * P_max * (1 - beta * <G / P_max>^n) <= gamma
    gamma = t_-1,   alpha = 3 * t_-1 / s_-1 - sqrt(3),   n = 1 / 2
    beta = R^n * (1 - s_-1 / f) / (1 - s_-1 / (sqrt(3) * t_-1))

where ``sqrt(J2,a)`` is the deviatoric amplitude (Eq. 11), ``P_max = (sigma_kk,a + sigma_kk,m)/3``
(Eq. 12), ``G = |grad P_max|`` is the hydrostatic gradient norm, ``<G / P_max> = 0`` when
``P_max <= 0``, and ``f`` is the constant-moment bending limit at radius ``R``.

Repository choices:
- Implemented as pointwise kernel functions (``pp_beta``, ``pp_index``) rather than a
  registered case model because experimental workbooks do not carry stress gradients ``G``.
- Valid for proportional (in-phase) loading (Eq. 8).
- Uses ``sqrt(J2,a) = s / sqrt(3)`` in reversed tension rather than the misprinted
  ``s * sqrt(3)`` after Eq. (13), which is required by Eq. (14).
"""

from __future__ import annotations

import numpy as np

from .stress_invariants import sqrt_j2_amplitude


def pp_beta(
    f_m1: float,
    t_m1: float,
    f: float,
    R: float,
    n: float = 0.5,
) -> float:
    """beta of Eq. (23) from the constant-moment bending limit f at radius R."""
    return R**n * (1.0 - f_m1 / f) / (1.0 - f_m1 / (np.sqrt(3.0) * t_m1))


def pp_index(
    S_mean,
    S_sin,
    G: float,
    f_m1: float,
    t_m1: float,
    beta: float,
    n: float = 0.5,
) -> float:
    """Ratio of the left-hand side of Eq. (17) to gamma for proportional loading
    S(t) = S_mean + S_sin f(t) with |f| <= 1."""
    S_mean = np.asarray(S_mean, dtype=float)
    S_sin = np.asarray(S_sin, dtype=float)
    alpha = 3.0 * t_m1 / f_m1 - np.sqrt(3.0)
    P_max = (abs(np.trace(S_sin)) + np.trace(S_mean)) / 3.0
    G_over_P_max = G / P_max if P_max > 0.0 else 0.0
    attenuation = 1.0 - beta * max(G_over_P_max, 0.0) ** n
    LHS = sqrt_j2_amplitude(S_sin, np.zeros((3, 3))) + alpha * P_max * attenuation
    return float(LHS / t_m1)
