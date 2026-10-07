# fatigue/models/cp_methods/mgse_zhu.py
from __future__ import annotations

"""MGSE (Zhu et al.) critical-plane damage parameter with material constant.

Zhu et al. (Int. J. Fatigue 2018) propose a modified generalized strain
energy (MGSE) criterion with an additional material constant k_MGSE:

    MGSE = τ_max * (Δγ/2) + k_MGSE * σ_{n,max} * (Δε_n/2)

We implement the left-hand damage parameter for calibration/assessment.

If k_MGSE is not provided in the material file, we follow the paper's
practical recommendation and fall back to k_MGSE = 1 with a warning.
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult
from ..utils.param_utils import get_param


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


class MGSE_Zhu:
    name = "MGSE_ZHU"

    def required_params(self) -> set[str]:
        # k_MGSE is optional (default 1.0)
        return set()

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        k_mgse = get_param(params, "k_MGSE", default=1.0, context=self.name)

        # Normal stress max (tensile only)
        sig0 = max(float(S0r[2, 2]), 0.0)
        sig1 = max(float(S1r[2, 2]), 0.0)
        sig_n_max = max(sig0, sig1)

        # Shear stress magnitude at max/min load step
        tau0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = float(params.get("_tau_max", max(tau0, tau1)))

        # Normal strain range on plane
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range on plane
        d_gam_a = float(E0r[0, 2] - E1r[0, 2])
        d_gam_b = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_a, d_gam_b))

        damage = float(tau_max * (0.5 * d_gam) + float(k_mgse) * sig_n_max * (0.5 * d_eps))
        metric = damage
        return CPPlaneResult(damage=damage, metric=metric)
