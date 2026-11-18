"""Vu, Halm & Nadot (2010) stress-invariant fatigue criterion.

Published in *International Journal of Fatigue* 32 (2010), pp. 1004-1014,
doi:10.1016/j.ijfatigue.2009.11.006 (roadmap N20), Eqs. (6)-(21)::

    f(t) = sqrt(gamma1 * J'_2(t)^2 + gamma2 * J_{2,mean}^2 + gamma3 * I_f(I_{1,a}, I_{1,m})) <= beta
    beta = t_-1,   gamma3 = (t_-1^2 - f_-1^2 / 3) / f_-1,   gamma1 + 4 * gamma2 / pi^2 = 1

where ``J'_2(t) = sqrt(S^a(t) : S^a(t) / 2)`` acts on the alternating deviatoric stress
``Sigma(t) - Sigma_m``, ``J_{2,mean}`` is its cycle mean, and ``I_{1,a}, I_{1,m}`` are
the amplitude and mean of ``I_1`` (with ``gamma1, gamma2`` and ``I_f`` branched at
``R_m = 750`` MPa, Eqs. 17-20).

Repository choices:
- For single-harmonic cycles the alternating deviator traces an ellipse with J2 semi-axes
  ``p >= q``, so ``max J'_2 = p`` and ``J_{2,mean} = (2 / pi) * p * E(1 - q^2 / p^2)``
  are evaluated in closed form via the arithmetic-geometric mean (``_elliptic_second_kind``).
- The ultimate-strength branch switches at ``Su < 750`` vs ``Su >= 750`` MPa, with finite
  positive ``Sigm1``, ``Taum1``, and ``Su`` required. Multi-frequency histories (Tables
  13-14) are outside the single-harmonic scope.
- Negative radicands under compressive mean stress clamp to ``0`` before taking the root,
  and the normalized index ``f_max / t_-1`` is returned (Eq. 21 error index is
  ``100 * (index - 1)``).
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from .stress_invariants import _deviatoric, harmonic_parts

R_M_SPLIT = 750.0


def _elliptic_second_kind(m: float) -> float:
    """Complete elliptic integral of the second kind E(m) for 0 <= m <= 1, via the AGM."""
    if m >= 1.0:
        return 1.0
    cur, comp, first = 1.0, np.sqrt(1.0 - m), m
    acc, weight = 0.5 * first, 0.5
    while abs(cur - comp) > 1e-15 * cur:
        cur, comp, diff = 0.5 * (cur + comp), np.sqrt(cur * comp), 0.5 * (cur - comp)
        weight *= 2.0
        acc += weight * diff * diff
    return float(np.pi / (2.0 * cur) * (1.0 - acc))


def j2_alternating(S_sin, S_cos) -> tuple[float, float]:
    """Largest J'2(t) over one period together with its period mean J2,mean (Eq. 7)."""
    S_sin_dev, S_cos_dev = _deviatoric(S_sin), _deviatoric(S_cos)
    g_ss, g_cc, g_sc = (
        0.5 * np.sum(x * y)
        for x, y in ((S_sin_dev, S_sin_dev), (S_cos_dev, S_cos_dev), (S_sin_dev, S_cos_dev))
    )
    mid, spread = 0.5 * (g_ss + g_cc), np.hypot(0.5 * (g_ss - g_cc), g_sc)
    hi2, lo2 = mid + spread, max(mid - spread, 0.0)
    if hi2 <= 0.0:
        return 0.0, 0.0
    hi = float(np.sqrt(hi2))
    return hi, 2.0 / np.pi * hi * _elliptic_second_kind(1.0 - lo2 / hi2)


def vhn_parameters(f_m1: float, t_m1: float, r_m: float) -> tuple[float, float, float, float]:
    """gamma1, gamma2, gamma3 and the weight carried by I1,m inside I_f."""
    gamma3 = (t_m1**2 - f_m1**2 / 3.0) / f_m1
    if r_m < R_M_SPLIT:
        return 0.65, 0.8636, gamma3, 1.0
    return 0.3, 1.7272, gamma3, f_m1 / t_m1


def vhn_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float, r_m: float) -> float:
    gamma1, gamma2, gamma3, mean_weight = vhn_parameters(f_m1, t_m1, r_m)
    J2_max, J2_mean = j2_alternating(S_sin, S_cos)
    I1_a = float(np.hypot(np.trace(np.asarray(S_sin, float)), np.trace(np.asarray(S_cos, float))))
    I1_m = float(np.trace(np.asarray(S_mean, float)))
    f_sq = gamma1 * J2_max**2 + gamma2 * J2_mean**2 + gamma3 * (I1_a + mean_weight * I1_m)
    return float(np.sqrt(max(f_sq, 0.0)) / t_m1)


class VuHalmNadot:
    """Vu, Halm and Nadot (2010) invariant criterion, returning f_max / t_-1 (Eq. 8)."""

    name = "VHN"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1", "Su"}

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
        supplied = [params.get(k) for k in ("Sigm1", "Taum1", "Su")]
        if any(val is None or not np.isfinite(float(val)) or float(val) <= 0 for val in supplied):
            raise ValueError("VHN needs finite, strictly positive Sigm1, Taum1 and Su (R_m branch)")
        f_m1, t_m1, r_m = (float(val) for val in supplied)
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": vhn_index(S_mean, S_sin, S_cos, f_m1, t_m1, r_m)})
