# fatigue/models/energy_based_methods/zhu_edp.py
from __future__ import annotations

"""Zhu et al. (2019) energy-based EDP inspired model (MOI-based).

The original paper (H. Zhu et al., Int. J. Fatigue 121 (2019) 1–8) proposes
an Equivalent Damage Parameter (EDP) based on:
  - Moment Of Inertia (MOI) features of the strain path
  - A non-proportional factor F_NP
  - Material constants (α_w, ξ)

The full formulation also uses the uniaxial elastoplastic work per cycle W0.
In this repository we typically **calibrate** a scalar damage parameter against
uniaxial data; therefore we implement a practical damage scalar that preserves:
  - Loading path effects via MOI
  - Phase effects via an F_NP proxy
  - Load level via the path diameter D

Implementation summary (based on Eq.(19)–(23) of the paper):

  1) Build the 2D strain path in the (ε_axial, γ_xy/√3) diagram using the
     single-harmonic representation stored in `harmonics`.
  2) Compute the perimeter centroid and MOI about the centroid numerically.
  3) Compute F1 = Ixx/D^2 + ξ * Iyy/D^2.
  4) Compute F_NP as |sin(Δφ)| from the relative phase between axial and shear.
  5) Return damage = F1 * (1 + α_w * F_NP) * D.

Defaults for (α_w, ξ) are taken from the 316L example in the paper.
"""

from typing import Mapping
import math
import numpy as np

from ..base import CaseResult
from ..utils.param_utils import get_param


class ZhuEDP:
    name = "ZHU_EDP"

    def required_params(self) -> set[str]:
        # α_w and ξ are optional (defaults provided). No mandatory params.
        return set()

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
        if harmonics is None:
            raise ValueError(f"{self.name}: harmonics are required to build the strain path")

        alpha_w = float(get_param(params, "alpha_w", default=0.38, context=self.name))
        xi = float(get_param(params, "xi_ZHU", default=0.15, context=self.name))

        E_mean = np.asarray(harmonics.get("E_mean"), dtype=float)
        E_sin = np.asarray(harmonics.get("E_sin"), dtype=float)
        E_cos = np.asarray(harmonics.get("E_cos"), dtype=float)

        # Axial strain component (yy) and engineering shear γ_xy = 2*ε_xy
        ex_m, ex_s, ex_c = float(E_mean[1, 1]), float(E_sin[1, 1]), float(E_cos[1, 1])
        g_m, g_s, g_c = 2.0 * float(E_mean[0, 1]), 2.0 * float(E_sin[0, 1]), 2.0 * float(E_cos[0, 1])
        gy_m, gy_s, gy_c = g_m / math.sqrt(3.0), g_s / math.sqrt(3.0), g_c / math.sqrt(3.0)

        # Relative phase proxy for F_NP
        def _phase(sin_coeff: float, cos_coeff: float) -> float:
            if abs(sin_coeff) < 1e-16 and abs(cos_coeff) < 1e-16:
                return 0.0
            return math.atan2(cos_coeff, sin_coeff)

        phi_x = _phase(ex_s, ex_c)
        phi_y = _phase(gy_s, gy_c)
        phi_rel = phi_y - phi_x
        F_NP = abs(math.sin(phi_rel))

        # Sample one cycle for numerical MOI
        n = 240
        t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        sin_t = np.sin(t)
        cos_t = np.cos(t)

        x = ex_m + ex_s * sin_t + ex_c * cos_t
        y = gy_m + gy_s * sin_t + gy_c * cos_t

        dxdt = ex_s * cos_t - ex_c * sin_t
        dydt = gy_s * cos_t - gy_c * sin_t
        ds = np.sqrt(dxdt * dxdt + dydt * dydt)
        # Avoid zero-length paths
        p = float(np.sum(ds))
        if p <= 1e-20:
            return CaseResult(values={"": 0.0})

        # Perimeter centroid
        Xc = float(np.sum(x * ds) / p)
        Yc = float(np.sum(y * ds) / p)

        x0 = x - Xc
        y0 = y - Yc

        Ixx = float(np.sum((y0 * y0) * ds) / p)
        Iyy = float(np.sum((x0 * x0) * ds) / p)

        # Convex enclosure diameter (approx. as 2*max radius from centroid)
        r = np.sqrt(x0 * x0 + y0 * y0)
        D = float(2.0 * np.max(r))
        if D <= 1e-20:
            return CaseResult(values={"": 0.0})

        F1 = (Ixx / (D * D)) + xi * (Iyy / (D * D))

        damage = float(F1 * (1.0 + alpha_w * F_NP) * D)
        return CaseResult(values={"": damage})
