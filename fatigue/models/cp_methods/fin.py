# fatigue_cp/methods/fin.py
from __future__ import annotations
from typing import Mapping, Tuple
import numpy as np
from ..base import CPPlaneResult


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    """Ensure A is a single (3,3) tensor."""
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


def _findley_plane_value(
    *,
    S0r: np.ndarray,
    S1r: np.ndarray,
    k_FI: float,
) -> Tuple[float, float]:
    """Findley value on one candidate plane.

    Same conventions as your previous cp_fun.cp_FIN:
      - Plane normal = local z-axis
      - Normal stress = Szz (signed maximum, including compression)
      - Shear stress range uses Sxz and Syz

    Returns: (damage, metric) where metric = damage.
    """
    S0r = _as33("S0r", S0r)
    S1r = _as33("S1r", S1r)
    k_FI = float(k_FI)

    if k_FI < 0.0:
        raise ValueError(f"k_FI must be >= 0. Got k_FI={k_FI}")

    # Signed maximum normal stress, as used in the Findley expression.
    sig0 = float(S0r[2, 2])
    sig1 = float(S1r[2, 2])
    smax = max(sig0, sig1)

    # tau_a = Sxz, tau_b = Syz
    dtau_a = float(S0r[0, 2] - S1r[0, 2])
    dtau_b = float(S0r[1, 2] - S1r[1, 2])
    delta_tau = float(np.sqrt(dtau_a**2 + dtau_b**2))

    fin_value = float(delta_tau / 2.0 + k_FI * smax)
    metric = fin_value
    return fin_value, metric


class Findley:
    name = "FIN"

    def required_params(self) -> set[str]:
        return {"k_FI"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        # Findley uses stresses only (E0r/E1r are ignored here, but kept for a uniform interface).
        val, metric = _findley_plane_value(S0r=S0r, S1r=S1r, k_FI=params["k_FI"])
        return CPPlaneResult(damage=float(val), metric=float(metric))