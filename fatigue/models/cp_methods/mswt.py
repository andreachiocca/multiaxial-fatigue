# fatigue/models/cp_methods/mswt.py
from __future__ import annotations

"""MSWT (Modified Smith-Watson-Topper) critical-plane damage parameter.

Zhu et al. (Int. J. Fatigue 2018) summarize Jiang's MSWT damage parameter:

    DP_MSWT = a * <σ_max> * (Δε/2) + (1-a) * (Δγ/2) * (Δτ/2)

where <x> is the MacCauley bracket (i.e. max(x,0)), Δε is the normal strain
range, σ_max is the maximum normal stress, Δγ is the shear strain range and
Δτ is the shear stress range on the candidate plane.

We implement the left-hand DP for calibration. The critical plane is selected
as the plane that maximizes DP.
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


class MSWT:
    name = "MSWT"

    def required_params(self) -> set[str]:
        # a_MSWT is optional (default 0.5)
        return set()

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        a = float(get_param(params, "a_MSWT", default=0.5, context=self.name))
        if not np.isfinite(a) or not 0.0 <= a <= 1.0:
            raise ValueError("a_MSWT must be finite and in [0, 1]")

        # sigma_max on plane (MacCauley bracket)
        sig0 = float(S0r[2, 2])
        sig1 = float(S1r[2, 2])
        sig_max = max(sig0, sig1)
        sig_pos = max(sig_max, 0.0)

        # normal strain range
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # shear strain range
        d_gam_a = float(E0r[0, 2] - E1r[0, 2])
        d_gam_b = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_a, d_gam_b))

        # shear stress range
        d_tau_a = float(S0r[0, 2] - S1r[0, 2])
        d_tau_b = float(S0r[1, 2] - S1r[1, 2])
        d_tau = float(np.hypot(d_tau_a, d_tau_b))

        term1 = a * sig_pos * (0.5 * d_eps)
        term2 = (1.0 - a) * (0.5 * d_gam) * (0.5 * d_tau)
        dp = float(term1 + term2)

        return CPPlaneResult(damage=dp, metric=dp)
