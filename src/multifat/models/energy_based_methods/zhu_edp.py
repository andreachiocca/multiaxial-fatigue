"""Energy-based equivalent damage parameter (EDP) of Zhu et al. (2019).

Published in *International Journal of Fatigue* 121 (2019), pp. 1-8 (accepted
manuscript p. 14), Eqs. (11)-(26)::

    psi = (xi * sqrt(12 * I_XX / D^2) + sqrt(12 * I_YY / D^2)) * (1 + alpha_w * F_NP) * dW_0
    F_NP = (sigma_eq / sigma_0 - 1) / (sigma_op / sigma_0 - 1)
    dW_0 = sigma_a^2 / (2 * E) + 2 * sigma_a * (2 * (sigma_a / H_c)^(1 / h_c)) * (1 - h_c) / (1 + h_c)

where ``I_XX, I_YY`` are centroidal second moments of the scaled strain path
``(eps_x * (1 + nu_bar), sqrt(3) * gamma_xy / 2)`` inside minimal enclosing circle
diameter ``D``, ``F_NP`` measures non-proportional cyclic hardening, and ``dW_0`` is
uniaxial elastoplastic work per cycle.

Repository choices:
- The case-level model ``ZhuEDP`` stays disabled via ``ensure_enabled("ZHU_EDP")``
  because workbooks supply elastic tensors only and lack stabilized non-proportional
  stress responses (``F_NP``) and cyclic constants (``H_c``, ``h_c``).
- Kernel functions (``strain_path``, ``moi_factor``, ``f_np``, ``delta_w0``, ``zhu_edp``)
  are exposed with explicit inputs, discretizing single-harmonic paths at 720 angles
  and defaulting ``nu_bar = 0.5`` in ``strain_path``.
"""

from __future__ import annotations

import math
from typing import Mapping

import numpy as np

from ..base import CaseResult

_PHASE = np.linspace(0.0, 2.0 * math.pi, 720, endpoint=False)


def strain_path(E_mean, E_sin, E_cos, nu_bar: float = 0.5, axial: int = 1, shear=(0, 1)):
    """Eq. (26) coordinates (X, Y) for a tension-torsion strain history with one harmonic."""
    E_mean, E_sin, E_cos = (np.asarray(a, dtype=float) for a in (E_mean, E_sin, E_cos))
    row, col = shear
    sin_t, cos_t = np.sin(_PHASE), np.cos(_PHASE)
    eps_axial = E_mean[axial, axial] + E_sin[axial, axial] * sin_t + E_cos[axial, axial] * cos_t
    gamma_eng = 2.0 * (E_mean[row, col] + E_sin[row, col] * sin_t + E_cos[row, col] * cos_t)
    return eps_axial * (1.0 + nu_bar), gamma_eng * math.sqrt(3.0) / 2.0


def moi_factor(X: np.ndarray, Y: np.ndarray, xi: float) -> tuple[float, float]:
    """Return F_1 from Eq. (21) and D for a closed path sampled at uniform times (Eqs. 19-21).

    A single-harmonic path is a centred ellipse (or a segment), so its smallest enclosing
    circle shares the ellipse centre and D equals twice the greatest distance from that centre.
    """
    dx, dy = np.roll(X, -1) - X, np.roll(Y, -1) - Y
    seg_len = np.hypot(dx, dy)
    perimeter = seg_len.sum()
    if perimeter <= 0.0:
        return 0.0, 0.0
    x_mid, y_mid = X + dx / 2, Y + dy / 2  # midpoint quadrature on each chord
    x_c, y_c = (x_mid * seg_len).sum() / perimeter, (y_mid * seg_len).sum() / perimeter
    # second moments of every chord about the centroid (exact for a polyline)
    I_x = ((((y_mid - y_c) ** 2) + dy**2 / 12.0) * seg_len).sum() / perimeter
    I_y = ((((x_mid - x_c) ** 2) + dx**2 / 12.0) * seg_len).sum() / perimeter
    diam = 2.0 * float(np.max(np.hypot(X - X.mean(), Y - Y.mean())))
    if diam <= 0.0:
        return 0.0, 0.0
    return float(xi * math.sqrt(12.0 * I_x) / diam + math.sqrt(12.0 * I_y) / diam), diam


def f_np(sigma_eq: float, sigma_0: float, sigma_op: float) -> float:
    """Evaluate Eq. (22) from measured von Mises stress amplitudes."""
    return (sigma_eq / sigma_0 - 1.0) / (sigma_op / sigma_0 - 1.0)


def delta_w0(eps_a: float, E: float, H_c: float, h_c: float) -> float:
    """Left-hand side of Eq. (12) for uniaxial loading:
    sigma_a^2/(2E) + dsigma deps_p (1 - h_c)/(1 + h_c).

    The stress amplitude sigma_a satisfies eps_a = sigma_a/E + (sigma_a/H_c)^(1/h_c)
    (Eq. 11 with F_NP = 0).
    """
    lower, upper = 0.0, E * eps_a
    for _ in range(100):
        trial = 0.5 * (lower + upper)
        lower, upper = (
            (trial, upper) if trial / E + (trial / H_c) ** (1.0 / h_c) < eps_a else (lower, trial)
        )
    sigma_a = 0.5 * (lower + upper)
    deps_p = 2.0 * (sigma_a / H_c) ** (1.0 / h_c)
    return sigma_a * sigma_a / (2.0 * E) + 2.0 * sigma_a * deps_p * (1.0 - h_c) / (1.0 + h_c)


def zhu_edp(F1: float, F_NP: float, dW0: float, alpha_w: float) -> float:
    """Evaluate Eq. (23)."""
    return F1 * (1.0 + alpha_w * F_NP) * dW0


class ZhuEDP:
    """Energy-based equivalent damage parameter model of Zhu et al."""

    name = "ZHU_EDP"

    def required_params(self) -> set[str]:
        return {"xi_ZHU", "alpha_w"}

    def evaluate_case(
        self,
        *,
        S0: np.ndarray,
        E0: np.ndarray,
        S1: np.ndarray,
        E1: np.ndarray,
        params: Mapping[str, float],
        R_list=None,
        harmonics=None,
    ) -> CaseResult:
        from multifat.models.references import ensure_enabled

        ensure_enabled("ZHU_EDP")
        raise NotImplementedError("ZHU_EDP requires measured F_NP and dW_0 values")
