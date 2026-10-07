# fatigue/models/cp_methods/mkbm.py
from __future__ import annotations

"""MKBM (Modified KBM) critical-plane strain parameter.

Li et al. (Int. J. Fatigue 2011) propose a modification of the
Kandil–Brown–Miller parameter to include additional cyclic hardening
through the maximum normal stress on the maximum shear strain plane:

    Δε*_eq / 2 = Δγ_max / 2 + (1 + σ_{n,max} / σ_y) * (Δε_n / 2)

We implement the left-hand strain-like parameter. The plane-selection
follows the original MKBM definition: the plane that maximizes Δγ (maximum
shear strain range plane).
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


class MKBM:
    name = "MKBM"

    def required_params(self) -> set[str]:
        # Requires yield strength (σ_y)
        return {"Sy"}

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        Sy = float(params["Sy"])
        if Sy <= 0:
            raise ValueError(f"Sy must be > 0. Got Sy={Sy}")

        # Normal stress max (tensile only)
        sig0 = max(float(S0r[2, 2]), 0.0)
        sig1 = max(float(S1r[2, 2]), 0.0)
        sig_n_max = max(sig0, sig1)

        # Normal strain range
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range
        d_gam_a = float(E0r[0, 2] - E1r[0, 2])
        d_gam_b = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_a, d_gam_b))

        dmg = float(0.5 * d_gam + (1.0 + sig_n_max / Sy) * (0.5 * d_eps))
        metric = float(d_gam)  # maximum shear strain plane
        return CPPlaneResult(damage=dmg, metric=metric)
