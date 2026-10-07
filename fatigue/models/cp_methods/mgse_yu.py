# fatigue/models/cp_methods/mgse_yu.py
from __future__ import annotations

"""MGSE (Yu et al.) critical-plane damage parameter.

From Yu et al. (Materials 2017) a modified generalized strain energy (MGSE)
parameter on a plane close to the maximum shear strain plane:

    MGSE = τ_max * (Δγ/2) + σ_{n,max} * (Δε/2)

where Δγ is the shear strain range on the candidate plane, Δε is the normal
strain range on that plane, τ_max is the maximum shear stress magnitude on
that plane, and σ_{n,max} is the maximum normal stress on that plane.

We use the same conventions as other CP methods in this repo:
  - Candidate plane normal = local z
  - Normal stress/strain uses zz component
  - Shear components use (xz, yz)
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


class MGSE_Yu:
    name = "MGSE_YU"

    def required_params(self) -> set[str]:
        # No extra material constants required.
        return set()

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        # Normal stress max (tensile only)
        sig0 = max(float(S0r[2, 2]), 0.0)
        sig1 = max(float(S1r[2, 2]), 0.0)
        sig_n_max = max(sig0, sig1)

        # Shear stress magnitude at max/min load step
        tau0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = max(abs(tau0), abs(tau1))

        # Normal strain range on plane
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range on plane
        d_gam_a = float(E0r[0, 2] - E1r[0, 2])
        d_gam_b = float(E0r[1, 2] - E1r[1, 2])
        d_gam = float(np.hypot(d_gam_a, d_gam_b))

        damage = float(tau_max * (0.5 * d_gam) + sig_n_max * (0.5 * d_eps))
        # In Yu's work the plane is chosen near maximum shear; in practice we
        # select the plane maximizing the MGSE value.
        metric = damage
        return CPPlaneResult(damage=damage, metric=metric)
