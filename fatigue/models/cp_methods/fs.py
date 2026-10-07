# fatigue_cp/methods/fs.py
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


def _fatemi_socie_plane_value(
    *,
    S0r: np.ndarray,
    E0r: np.ndarray,
    S1r: np.ndarray,
    E1r: np.ndarray,
    k_FS: float,
    Sy: float,
) -> Tuple[float, float]:
    """Compute the Fatemi–Socie parameter on a single candidate plane.

    Conventions (same as your previous cp_fun.cp_FS implementation):
      - Normal direction of the candidate plane is the local z-axis.
      - Normal stress is Szz (clipped at 0 when taking the maximum).
      - Shear strain range is based on Exz and Eyz components.

    Returns
    -------
    damage : float
        FS value for this plane.
    metric : float
        Plane-selection metric (DeltaGamma).
    """

    S0r = _as33("S0r", S0r)
    E0r = _as33("E0r", E0r)
    S1r = _as33("S1r", S1r)
    E1r = _as33("E1r", E1r)

    k_FS = float(k_FS)
    Sy = float(Sy)

    if Sy <= 0.0:
        raise ValueError(f"Sy must be > 0. Got Sy={Sy}")
    if k_FS < 0.0:
        raise ValueError(f"k_FS must be >= 0. Got k_FS={k_FS}")

    # sigma_nn = Szz (clipped at 0)
    sig0 = max(float(S0r[2, 2]), 0.0)
    sig1 = max(float(S1r[2, 2]), 0.0)
    smax = max(sig0, sig1)

    # Shear strain range in the plane: gamma_a = Exz, gamma_b = Eyz
    dgam_a = float(E0r[0, 2] - E1r[0, 2])
    dgam_b = float(E0r[1, 2] - E1r[1, 2])
    delta_gamma = float(np.sqrt(dgam_a**2 + dgam_b**2))

    fs_value = float(delta_gamma * (1.0 + k_FS * (smax / Sy)))

    # Keep the same plane-selection metric you were using before.
    metric = delta_gamma

    return fs_value, metric


class FatemiSocie:
    name = "FS"

    def required_params(self) -> set[str]:
        return {"Sy", "k_FS"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        val, metric = _fatemi_socie_plane_value(
            S0r=S0r,
            E0r=E0r,
            S1r=S1r,
            E1r=E1r,
            k_FS=params["k_FS"],
            Sy=params["Sy"],
        )
        return CPPlaneResult(damage=float(val), metric=float(metric))