# fatigue/models/cp_methods/liu2021.py
from __future__ import annotations

"""Liu et al. (2021) critical-plane model (TNP family) – damage parameter.

From Liu et al. (Int. J. Pressure Vessels and Piping 194 (2021) 104532),
Eq. (27) defines a modified shear strain amplitude on the critical plane:

    γ*_max / 2 = (Δγ_max / 2) * ( 1 + σ_{n,max} * sqrt(σ_{n,max} E Δε)
                                      / (2 G Δγ_max τ'_f) )

We implement γ*_max/2 as a damage parameter.

Notes on implementation:
  - Δγ_max is the shear strain range on the candidate plane.
  - Δε is taken here as the normal strain range on the same plane.
  - τ'_f is the shear fatigue strength coefficient. If not provided, it is
    derived from axial fatigue properties using von Mises relations.
  - Plane selection follows the paper's definition: maximum shear strain plane.
"""

from typing import Mapping
import math
import numpy as np

from ..base import CPPlaneResult
from ..utils.fatigue_estimation import ensure_shear_fatigue_props


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


class Liu2021:
    name = "LIU2021"

    def required_params(self) -> set[str]:
        # Need E and nu to compute G; Su is used for estimating fatigue props if needed.
        return {"E", "nu", "Su"}

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        E = float(params["E"])  # MPa
        nu = float(params["nu"])
        if E <= 0:
            raise ValueError(f"E must be > 0. Got E={E}")
        if not (-1.0 < nu < 0.5):
            raise ValueError(f"nu not plausible. Got nu={nu}")
        G = E / (2.0 * (1.0 + nu))

        tau_f = params.get("tau_f")
        if tau_f is None:
            tau_f = ensure_shear_fatigue_props(params).tau_f
        tau_f = float(tau_f)
        if tau_f <= 0:
            raise ValueError(f"tau_f must be > 0. Got tau_f={tau_f}")

        # Normal stress max (tensile only)
        sig0 = max(float(S0r[2, 2]), 0.0)
        sig1 = max(float(S1r[2, 2]), 0.0)
        sig_n_max = max(sig0, sig1)

        # Normal strain range (used as Δε)
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range Δγ
        d_gam_a = float(E0r[0, 2] - E1r[0, 2])
        d_gam_b = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_a, d_gam_b))

        # Algebraically expanded Eq. (27), including its zero-shear limit.
        rad = max(sig_n_max * E * d_eps, 0.0)
        gamma_star_over2 = 0.5*d_gam + sig_n_max*math.sqrt(rad)/(4.0*G*tau_f)

        metric = d_gam
        return CPPlaneResult(damage=float(gamma_star_over2), metric=float(metric))
