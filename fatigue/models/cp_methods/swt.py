# fatigue_cp/methods/swt.py
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

class SWT:
    name = "SWT"

    def required_params(self) -> set[str]:
        return set()

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        """Smith-Watson-Topper value for one candidate plane.

        Same convention as cp_fun.cp_SWT:
          - plane normal is local z
          - normal stress uses Szz (clipped at 0) and Smax = max(Szz@step1, Szz@step2)
          - normal strain range uses |Ezz(step1) - Ezz(step2)|

        Returns (damage, metric) where metric = DeltaEps.
        """
        S0r = _as33("S0r", S0r)
        E0r = _as33("E0r", E0r)
        S1r = _as33("S1r", S1r)
        E1r = _as33("E1r", E1r)

        # sigma_nn = Szz (clipped at 0)
        sig0 = max(float(S0r[2, 2]), 0.0)
        sig1 = max(float(S1r[2, 2]), 0.0)
        smax = max(sig0, sig1)

        # eps_nn = Ezz
        delta_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        swt_value = float(delta_eps * smax / 2.0)
        metric = float(delta_eps)
        return CPPlaneResult(damage=float(swt_value), metric=float(metric))