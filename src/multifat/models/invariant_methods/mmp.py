"""Manson-McKnight-Papuga (MMP) fatigue limit criterion of Papuga & Fojtik (2017) (G05).

J. Papuga, F. Fojtik, Int. J. Fatigue 104 (2017) 27-42, doi:10.1016/j.ijfatigue.2017.07.001,
Eqs. (19)-(25), variant MMP_WP (Walker mean-stress correction with the factor of Eq. 22)::

    sigma_a = sqrt(1/2 [(s_x,a - s_y,a)^2 + (s_y,a - s_z,a)^2 + (s_z,a - s_x,a)^2
                        + 2 kappa^2 (t_xy,a^2 + t_yz,a^2 + t_zx,a^2)])                     (19)
    sigma_m = sqrt(1/2 [(same with mean values)] + 2 X_m^2 (t_xy,m^2 + t_yz,m^2 + t_zx,m^2))  (24)
    sigma_a^w (sigma_a + beta_P sigma_m)^(1-w) <= f_-1                                    (20)
    beta_P = s_1,max / (s_1,max - s_3,min)  if |s_1,max| >= |s_3,min|, else
             s_3,min / (s_1,max - s_3,min)                                                 (22)
    w = log(f_0 / f_-1) / log 2,   X_m = 2 kappa [(2 t_-1 / t_0)^(1/(1-w)) - 1]          (21), (25)

kappa = f_-1/t_-1; f_0 and t_0 are the repeated (R = 0) limits as maximum stress of the cycle;
s_1,max and s_3,min are the largest and smallest principal stresses over the whole cycle. The
amplitudes and means are taken component by component, so the criterion ignores the phase shift
between load channels (the paper's stated limitation). The error index is
Delta FI = 100 (D_p - f_-1)/f_-1 (Eq. 27). The fully reversed axial and torsion limits, repeated
tension (beta_P = 1, Eq. 21) and repeated torsion (beta_P = 1/2, Eq. 25) are met exactly.

Repository choices:
  - DP = D_p / f_-1, the fatigue index (1 on the fatigue limit, Delta FI = 100 (DP - 1));
  - maintainer decision: missing f_0 (Sig0) is estimated from Su and f_-1 by Boehme et al. (2026)
    Eq. (1), as for PCR and MMMH; missing t_0 (Tau0) from 4 t_-1/t_0 - 2 f_-1/f_0 = 1 (Zenner,
    quoted as the paper's Eq. 23);
  - maintainer decision: with harmonics, each component amplitude is sqrt(sin^2 + cos^2) with the
    sign of its in-phase (sin) part, so proportional loading is exact and the criterion stays
    phase-blind for non-proportional loading;
  - a compressive mean making sigma_a + beta_P sigma_m negative gives D_p = 0 (Eq. 20 undefined);
  - life from the shared master curve DP = A N^b (an extension of the fatigue-limit criterion).

User guide: docs/method-guides/MMP.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from .stress_invariants import harmonic_parts, peak_principal_stress

_SHEAR = ((0, 1), (1, 2), (2, 0))


def mmp_constants(f_m1: float, t_m1: float, f_0: float, t_0: float) -> tuple[float, float, float]:
    """kappa, the Walker exponent w (Eq. 21) and the mean shear weight X_m (Eq. 25)."""
    kappa = f_m1 / t_m1
    w = np.log(f_0 / f_m1) / np.log(2.0)
    if not 0.0 < w < 1.0:
        raise ValueError(
            f"MMP: Walker exponent w = {w:.4g} outside (0, 1) needs f_-1 < f_0 < 2 f_-1"
        )
    X_m = 2.0 * kappa * ((2.0 * t_m1 / t_0) ** (1.0 / (1.0 - w)) - 1.0)
    return float(kappa), float(w), float(X_m)


def mmp_equivalent(S, shear_weight: float) -> float:
    """sqrt of Eq. (19) / (24) for a tensor of component amplitudes or means."""
    S = np.asarray(S, dtype=float)
    sx, sy, sz = S[0, 0], S[1, 1], S[2, 2]
    shear = sum(S[i, j] ** 2 for i, j in _SHEAR)
    return float(
        np.sqrt(
            0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2 + 2.0 * shear_weight**2 * shear)
        )
    )


def mmp_beta(s_1max: float, s_3min: float) -> float:
    """Mean stress factor beta_P of Eq. (22)."""
    span = s_1max - s_3min
    if span <= 0.0:
        return 0.0
    return (s_1max if abs(s_1max) >= abs(s_3min) else s_3min) / span


def mmp_damage(sigma_a: float, sigma_m: float, beta: float, w: float) -> float:
    """D_p, the left-hand side of Eq. (20) (MPa); 0 when sigma_a + beta sigma_m <= 0."""
    base = sigma_a + beta * sigma_m
    if base <= 0.0 or sigma_a <= 0.0:
        return 0.0
    return float(sigma_a**w * base ** (1.0 - w))


def signed_amplitudes(S_sin, S_cos) -> np.ndarray:
    """Component amplitudes sqrt(sin^2 + cos^2) with the sign of the in-phase part."""
    S_sin, S_cos = np.asarray(S_sin, dtype=float), np.asarray(S_cos, dtype=float)
    sign = np.where(S_sin != 0.0, np.sign(S_sin), np.sign(S_cos))
    return sign * np.hypot(S_sin, S_cos)


def mmp_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float, f_0: float, t_0: float) -> float:
    """D_p / f_-1 for one single-harmonic cycle."""
    kappa, w, X_m = mmp_constants(f_m1, t_m1, f_0, t_0)
    S_mean, S_sin, S_cos = (np.asarray(x, dtype=float) for x in (S_mean, S_sin, S_cos))
    sigma_a = mmp_equivalent(signed_amplitudes(S_sin, S_cos), kappa)
    sigma_m = mmp_equivalent(S_mean, X_m)
    s_1max = peak_principal_stress(S_mean, S_sin, S_cos)
    s_3min = -peak_principal_stress(-S_mean, -S_sin, -S_cos)
    return mmp_damage(sigma_a, sigma_m, mmp_beta(s_1max, s_3min), w) / f_m1


class MMP:
    """Manson-McKnight-Papuga criterion MMP_WP (see the module docstring)."""

    name = "MMP"
    optional_params = {
        "Sig0": "Boehme et al. 2026 Eq. (1) from Su and Sigm1",
        "Tau0": "4 Taum1/Tau0 - 2 Sigm1/Sig0 = 1 (Zenner; Papuga & Fojtik 2017 Eq. 23)",
        "Su": "only to estimate Sig0",
    }

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        from ..integral_methods.bohme_papuga import _strengths

        f_m1, t_m1, f_0, t_0 = _strengths(params, self.name)
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": mmp_index(S_mean, S_sin, S_cos, f_m1, t_m1, f_0, t_0)})
