# fatigue/models/cp_methods/gse.py
from __future__ import annotations

"""GSE (Generalized Strain Energy) critical-plane damage parameter.

Ince & Glinka (Int. J. Fatigue 62 (2014) 34–41) define the generalized
strain energy (GSE) damage parameter on the *maximum damage plane*:

  GSE = [
          τ_max * (Δγ_e/2)  +  (Δτ/2) * (Δγ_p/2)
        + σ_{n,max} * (Δε_{n,e}/2) + (Δσ_n/2) * (Δε_{n,p}/2)
        ]_max

where elastic/plastic parts are computed from the stress/strain response
on each candidate plane.

Repository convention:
  - candidate plane normal is the local z-axis
  - normal components use zz
  - in-plane shear uses (xz, yz) as a 2D shear vector

Because experimental files store total strains, we approximate elastic
strains using isotropic Hooke's law with (E, nu) and set plastic strain as:
  ε_p = ε_total - ε_e

This implementation returns the left-hand GSE scalar; the case-level
pipeline calibrates this scalar against uniaxial fatigue data.
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


def _elastic_strain_from_stress_iso(sig: np.ndarray, *, E: float, nu: float) -> np.ndarray:
    """Small-strain isotropic compliance: ε = (1/E)[(1+ν)σ - ν tr(σ) I]."""
    sig = np.asarray(sig, dtype=float)
    tr = float(np.trace(sig))
    I = np.eye(3)
    return ((1.0 + float(nu)) / float(E)) * sig - (float(nu) / float(E)) * tr * I


class GSE:
    name = "GSE"

    def required_params(self) -> set[str]:
        # Needed for elastic/plastic strain decomposition.
        return {"E", "nu"}

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        E = float(params["E"])
        nu = float(params["nu"])
        if E <= 0.0:
            raise ValueError(f"{self.name}: E must be > 0. Got E={E}")
        if not (-1.0 < nu < 0.5):
            raise ValueError(f"{self.name}: nu must be in (-1,0.5). Got nu={nu}")

        # --- elastic/plastic strain tensors (two load steps) ---
        Ee0 = _elastic_strain_from_stress_iso(S0r, E=E, nu=nu)
        Ee1 = _elastic_strain_from_stress_iso(S1r, E=E, nu=nu)
        Ep0 = E0r - Ee0
        Ep1 = E1r - Ee1

        # --- normal stress max (tension only) ---
        sig_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)

        # --- shear stress max magnitude ---
        tau0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = max(tau0, tau1)

        # stress ranges
        d_tau = float(np.hypot(S0r[0, 2] - S1r[0, 2], S0r[1, 2] - S1r[1, 2]))
        d_sig_n = abs(float(S0r[2, 2] - S1r[2, 2]))

        # strain ranges (elastic/plastic)
        d_gam_e = float(np.hypot(Ee0[0, 2] - Ee1[0, 2], Ee0[1, 2] - Ee1[1, 2]))
        d_gam_p = float(np.hypot(Ep0[0, 2] - Ep1[0, 2], Ep0[1, 2] - Ep1[1, 2]))

        d_eps_e = abs(float(Ee0[2, 2] - Ee1[2, 2]))
        d_eps_p = abs(float(Ep0[2, 2] - Ep1[2, 2]))

        # GSE (energy-like) damage
        gse = (
            tau_max * (0.5 * d_gam_e)
            + (0.5 * d_tau) * (0.5 * d_gam_p)
            + sig_n_max * (0.5 * d_eps_e)
            + (0.5 * d_sig_n) * (0.5 * d_eps_p)
        )
        gse = float(gse)

        # Critical plane: maximum damage parameter plane -> metric = damage
        return CPPlaneResult(damage=gse, metric=gse)
