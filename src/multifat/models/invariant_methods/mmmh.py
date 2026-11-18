"""Modified Mohr-Mises hypothesis (MMMH) of Mertens & Liebich (2017).

Published in *International Journal of Fatigue* (2017), doi:10.1016/j.ijfatigue.2017.06.003
(roadmap G07), Eqs. (1), (4), (9), (22)-(25) for smooth surfaces (``m_a = n_el = m_o1 = m_a1 = 1``)::

    sigma_eq,MMM(eta) = V(eta) * sqrt((4 - k_a^2) * M_M^2 + k_a^2 * R_M^2),   k_a = s_-1 / t_-1
    M_M = (sigma_1 + sigma_2) / 2 - sigma_3,   |R_M| = |sigma_1 - sigma_2| / 2
    V   = sign(sqrt(4 - k_a^2) * M_M + k_a * R_M)
    sigma_eq,a = (max sigma_eq,MMM - min sigma_eq,MMM) / 2
    sigma_eq,m = sigma_o1 - sigma_a1

For one harmonic, the mean-free stress flips sign every half period, so
``sigma_eq,a = max_t sqrt(I_1^2 - k_a^2 I_2)`` reduces to ``COIN`` with ``k_phi = 0.5``,
and ``sigma_eq,m`` equals the difference between the cycle peaks of maximum principal
total stress and alternating principal stress (Eqs. 20-24).

Repository choices:
- Because the paper's load-cycle-dependent Haigh diagrams rely on unobtainable internal
  references ([16, 23]), the registered ``MMMH`` model evaluates the endurance-limit Haigh
  line ``(sigma_eq,a + M * sigma_eq,m) / s_-1`` for universal master-curve life prediction.
- Mean-stress sensitivity ``M`` uses optional ``M_MMMH`` when supplied, or
  ``M = 2 * Sigm1 / Sig0 - 1`` (Dubbel 2004, E1) using ``Sig0`` or the Boehme et al. (2026)
  Eq. (1) estimate from ``Su`` and ``Sigm1``. Requires ``0 < Sigm1 / Taum1 <= 2``.
- Time-domain sign selection for ``R_M`` drops the constant global quotient of Eq. (20),
  which cancels in the single-harmonic half-range; Eq. (8)'s text misprint
  ``|s1a - s2a|`` for ``k_a^2 = 0`` (actually ``|s1a + s2a|``) is not followed.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from .coin_liwi import coin_amplitude
from .stress_invariants import harmonic_parts, peak_principal_stress


def mmmh_amplitude(S_sin, S_cos, k_a: float) -> float:
    """sigma_eq,a of Eq. (25) for one harmonic (= COIN with k_phi = 0.5)."""
    return coin_amplitude(S_sin, S_cos, k_a, k_phi=0.5)


def mmmh_equivalent_history(sigma_x, sigma_y, tau_xy, k_a: float, sigma_z=0.0) -> np.ndarray:
    """sigma_eq,MMM(eta) of Eq. (1) for sampled mean-free surface stresses (smooth, m_a = n_el = 1)."""
    sigma_x, sigma_y, tau_xy = (
        np.asarray(comp, dtype=float) for comp in (sigma_x, sigma_y, tau_xy)
    )
    sigma_z = np.broadcast_to(np.asarray(sigma_z, dtype=float), sigma_x.shape)
    M_M = 0.5 * (sigma_x + sigma_y) - sigma_z
    half_diff = 0.5 * (sigma_x - sigma_y)
    shear_mag = np.hypot(half_diff, tau_xy)
    idx_min = int(np.argmin(shear_mag))
    angle = 0.5 * np.arctan2(tau_xy[idx_min], half_diff[idx_min])
    shear_rot = -half_diff * np.sin(2.0 * angle) + tau_xy * np.cos(2.0 * angle)
    R_M = np.sign(shear_rot) * shear_mag
    magnitude = np.sqrt(np.maximum((4.0 - k_a**2) * M_M**2 + k_a**2 * R_M**2, 0.0))
    return np.sign(np.sqrt(4.0 - k_a**2) * M_M + k_a * R_M) * magnitude


def mmmh_mean_proportional(S_mean, S_amp) -> float:
    """sigma_eq,m of Eq. (24) for proportional loading S(t) = S_mean + S_amp f(t), |f| <= 1."""
    S_mean, S_amp = np.asarray(S_mean, dtype=float), np.asarray(S_amp, dtype=float)
    top_mean = max(
        np.linalg.eigvalsh(S_mean + S_amp)[-1],
        np.linalg.eigvalsh(S_mean - S_amp)[-1],
    )
    half_range = float(np.max(np.abs(np.linalg.eigvalsh(S_amp))))
    return float(top_mean - half_range)


def mmmh_mean(S_mean, S_sin, S_cos) -> float:
    """sigma_eq,m of Eq. (24) (smooth: m_o1 = m_a1 = n_el = 1) for a single-harmonic cycle."""
    S_mean, S_sin, S_cos = (np.asarray(comp, dtype=float) for comp in (S_mean, S_sin, S_cos))
    peak_total = peak_principal_stress(S_mean, S_sin, S_cos)  # Eq. (22)
    peak_oscillating = peak_principal_stress(np.zeros((3, 3)), S_sin, S_cos)  # Eq. (23), see above
    return float(peak_total - peak_oscillating)


def mmmh_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float, M: float) -> float:
    """(sigma_eq,a + M sigma_eq,m)/s_-1 on the endurance-limit Haigh line, k_a = s_-1/t_-1."""
    k_a = f_m1 / t_m1
    if not 0.0 < k_a <= 2.0:
        raise ValueError("MMMH needs 0 < k_a = Sigm1/Taum1 <= 2 (Eq. 1)")
    return (mmmh_amplitude(S_sin, S_cos, k_a) + M * mmmh_mean(S_mean, S_sin, S_cos)) / f_m1


class MMMH:
    """Modified Mohr-Mises hypothesis at the endurance limit (see the module docstring)."""

    name = "MMMH"
    optional_params = {
        "M_MMMH": "2 Sigm1/Sig0 - 1 (Dubbel 2004, E1: mean-stress sensitivity between R = -1 and 0)",
        "Sig0": "Boehme et al. 2026 Eq. (1) from Su and Sigm1",
        "Su": "only to estimate Sig0",
    }

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        M = params.get("M_MMMH")
        if M is None:
            from ..integral_methods.bohme_papuga import _strengths

            f_m1, _, f_0, _ = _strengths(params, self.name)
            M = 2.0 * f_m1 / f_0 - 1.0
        index = mmmh_index(
            S_mean,
            S_sin,
            S_cos,
            float(params["Sigm1"]),
            float(params["Taum1"]),
            float(M),
        )
        return CaseResult(values={"": index})
