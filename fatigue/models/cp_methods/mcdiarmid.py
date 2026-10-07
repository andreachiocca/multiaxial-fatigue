# fatigue/models/cp_methods/mcdiarmid.py
from __future__ import annotations

from typing import Mapping, Tuple
import numpy as np

from ..base import CPPlaneResult


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


def _mcdiarmid_plane_value(
    *,
    S0r: np.ndarray,
    S1r: np.ndarray,
    Su: float,
    Taum1: float,
) -> Tuple[float, float]:
    """
    McDiarmid criterion (common HCF form) on one candidate plane.

    Critical plane: maximum shear stress amplitude/range (MSSR-style).

    A typical normalized statement is:
      tau_a / Taum1 + sigma_n,max / (2*Su) = 1
    Rearranged to a "damage-like" scalar (monotone with severity):
      tau_a + k * sigma_n,max
      where k = Taum1 / (2*Su)

    Conventions:
      - candidate plane normal is local z-axis
      - shear stress components on plane: Sxz and Syz
      - normal stress: Szz (clipped at 0 when taking maximum)
      - two load steps: (0=max), (1=min)
    """
    S0r = _as33("S0r", S0r)
    S1r = _as33("S1r", S1r)

    Su = float(Su)
    Taum1 = float(Taum1)

    if Su <= 0.0:
        raise ValueError(f"Su (ultimate tensile strength) must be > 0. Got {Su}")
    if Taum1 <= 0.0:
        raise ValueError(f"Taum1 (torsional fatigue limit) must be > 0. Got {Taum1}")

    # k from existing material properties (no new Excel fields)
    k = Taum1 / (2.0 * Su)

    # Shear stress RANGE vector components on plane (xz, yz)
    dtau_xz = float(S0r[0, 2] - S1r[0, 2])
    dtau_yz = float(S0r[1, 2] - S1r[1, 2])
    delta_tau = float(np.sqrt(dtau_xz**2 + dtau_yz**2))

    # shear stress amplitude
    tau_a = 0.5 * delta_tau

    # max normal stress on plane (clip compression)
    sig0 = max(float(S0r[2, 2]), 0.0)
    sig1 = max(float(S1r[2, 2]), 0.0)
    sigma_n_max = max(sig0, sig1)

    damage = float(tau_a + k * sigma_n_max)

    # MSSR metric: maximize shear stress range (or amplitude)
    metric = delta_tau
    return damage, metric


class McDiarmid:
    """
    McDiarmid critical plane method.
    Computes k internally from:
      k = Taum1 / (2*Su)
    """
    name = "MCD"

    def required_params(self) -> set[str]:
        return {"Su", "Taum1"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,  # unused
        S1r,
        E1r,  # unused
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        if params.get("Su", None) is None:
            raise KeyError("MCD requires params['Su'] (Tensile_strength_MPa)")
        if params.get("Taum1", None) is None:
            raise KeyError("MCD requires params['Taum1'] (Torsional_fatigue_limit_MPa)")

        damage, metric = _mcdiarmid_plane_value(
            S0r=S0r,
            S1r=S1r,
            Su=float(params["Su"]),
            Taum1=float(params["Taum1"]),
        )
        return CPPlaneResult(damage=damage, metric=metric)