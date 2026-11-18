"""GSE (Generalized Strain Energy) critical-plane damage parameter.

Ince & Glinka (Int. J. Fatigue 62 (2014) 34–41) define the generalized
strain energy (GSE) damage parameter on the *maximum damage plane*:

  GSE = [
          τ_max * (Δγ_e/2)  +  (Δτ/2) * (Δγ_p/2)
        + σ_{n,max} * (Δε_{n,e}/2) + (Δσ_n/2) * (Δε_{n,p}/2)
        ]_max

where elastic/plastic parts are computed from the stress/strain response
on each candidate plane.

Repository choices:
  - candidate plane normal is the local z-axis
  - normal components use zz
  - in-plane shear uses (xz, yz) as a 2D shear vector

The workbooks store total strains; elastic strains follow from isotropic
Hooke's law with (E, nu) and plastic strain is the difference:
  ε_p = ε_total - ε_e

This implementation returns the left-hand GSE scalar; the case-level
pipeline calibrates this scalar against uniaxial fatigue data.

User guide: docs/method-guides/GSE.md.
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, elastic_strain_planes, tau_max_planes

from ..base import CPPlaneResult


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

    def evaluate_on_plane(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]
    ) -> CPPlaneResult:
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
        sigma_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)

        # --- shear stress max magnitude ---
        tau_0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau_1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = float(params.get("_tau_max", max(tau_0, tau_1)))

        # stress ranges
        d_tau = float(np.hypot(S0r[0, 2] - S1r[0, 2], S0r[1, 2] - S1r[1, 2]))
        d_sigma_n = abs(float(S0r[2, 2] - S1r[2, 2]))

        # strain ranges (elastic/plastic)
        d_gam_e = 2.0 * float(np.hypot(Ee0[0, 2] - Ee1[0, 2], Ee0[1, 2] - Ee1[1, 2]))
        d_gam_p = 2.0 * float(np.hypot(Ep0[0, 2] - Ep1[0, 2], Ep0[1, 2] - Ep1[1, 2]))

        d_eps_e = abs(float(Ee0[2, 2] - Ee1[2, 2]))
        d_eps_p = abs(float(Ep0[2, 2] - Ep1[2, 2]))
        if "_strain_ranges" in params:
            d_gam_e, d_gam_p, d_eps_e, d_eps_p = params["_strain_ranges"]

        # GSE (energy-like) damage
        gse = (
            tau_max * (0.5 * d_gam_e)
            + (0.5 * d_tau) * (0.5 * d_gam_p)
            + sigma_n_max * (0.5 * d_eps_e)
            + (0.5 * d_sigma_n) * (0.5 * d_eps_p)
        )
        gse = float(gse)

        # Critical plane: maximum damage parameter plane -> metric = damage
        return CPPlaneResult(damage=gse, metric=gse)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        E, nu = float(params["E"]), float(params["nu"])
        if E <= 0.0:
            raise ValueError(f"{self.name}: E must be > 0. Got E={E}")
        if not (-1.0 < nu < 0.5):
            raise ValueError(f"{self.name}: nu must be in (-1,0.5). Got nu={nu}")
        Ee0 = elastic_strain_planes(S0r, E=E, nu=nu)
        Ee1 = elastic_strain_planes(S1r, E=E, nu=nu)
        Ep0, Ep1 = E0r - Ee0, E1r - Ee1
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], S1r[:, 2, 2]), 0.0)
        tau_max_plane = tau_max_planes(S0r, S1r, params, tau_max)
        d_gam_e = 2.0 * np.hypot(Ee0[:, 0, 2] - Ee1[:, 0, 2], Ee0[:, 1, 2] - Ee1[:, 1, 2])
        d_gam_p = 2.0 * np.hypot(Ep0[:, 0, 2] - Ep1[:, 0, 2], Ep0[:, 1, 2] - Ep1[:, 1, 2])
        d_eps_e = np.abs(Ee0[:, 2, 2] - Ee1[:, 2, 2])
        d_eps_p = np.abs(Ep0[:, 2, 2] - Ep1[:, 2, 2])
        if strain_ranges is not None:
            d_gam_e, d_gam_p, d_eps_e, d_eps_p = (np.asarray(x, dtype=float) for x in strain_ranges)
        elif "_strain_ranges" in params:
            d_gam_e, d_gam_p, d_eps_e, d_eps_p = (
                np.full(S0r.shape[0], float(x)) for x in params["_strain_ranges"]
            )
        d_tau = np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        d_sigma_n = np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
        gse = (
            tau_max_plane * (0.5 * d_gam_e)
            + (0.5 * d_tau) * (0.5 * d_gam_p)
            + sigma_n_max * (0.5 * d_eps_e)
            + (0.5 * d_sigma_n) * (0.5 * d_eps_p)
        )
        return gse, gse
