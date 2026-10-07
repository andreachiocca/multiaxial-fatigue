# fatigue/models/cp_methods/gsa.py
from __future__ import annotations

"""GSA (Generalized Strain Amplitude) critical-plane damage parameter.

Ince & Glinka (Int. J. Fatigue 62 (2014) 34–41) convert the generalized
strain energy (GSE) form into a generalized strain amplitude (GSA) form by
normalizing the *elastic* strain amplitudes with stress correction factors.

In a practical form (see Eq.(5) in the paper):

  GSA = [
          (τ_max/τ'_f) * (Δγ_e/2) + (Δγ_p/2)
        + (σ_{n,max}/σ'_f) * (Δε_{n,e}/2) + (Δε_{n,p}/2)
        ]_max

where σ'_f and τ'_f are uniaxial fatigue strength coefficients.

When fatigue properties are missing from the material Excel file we estimate
them using the *non-hardness* empirical relations summarized in Zhao et al.
(2020) (Uniform Material Law). Shear properties are then derived from axial
properties using common von-Mises relations.
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult
from ..utils.fatigue_estimation import ensure_axial_fatigue_props, ensure_shear_fatigue_props


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


def _elastic_strain_from_stress_iso(sig: np.ndarray, *, E: float, nu: float) -> np.ndarray:
    sig = np.asarray(sig, dtype=float)
    tr = float(np.trace(sig))
    I = np.eye(3)
    return ((1.0 + float(nu)) / float(E)) * sig - (float(nu) / float(E)) * tr * I


class GSA:
    name = "GSA"

    def required_params(self) -> set[str]:
        # E, nu required for elastic/plastic split; Su is used if fatigue properties are missing.
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

        # Ensure fatigue properties
        axial = ensure_axial_fatigue_props(params)
        shear = ensure_shear_fatigue_props(params)
        sigma_f = float(axial.sigma_f)
        tau_f = float(shear.tau_f)
        if sigma_f <= 0.0 or tau_f <= 0.0:
            raise ValueError(f"{self.name}: invalid fatigue strength coefficients sigma_f={sigma_f}, tau_f={tau_f}")

        # elastic/plastic strain tensors
        Ee0 = _elastic_strain_from_stress_iso(S0r, E=E, nu=nu)
        Ee1 = _elastic_strain_from_stress_iso(S1r, E=E, nu=nu)
        Ep0 = E0r - Ee0
        Ep1 = E1r - Ee1

        # maxima
        sig_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)
        tau0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = max(tau0, tau1)

        # elastic/plastic strain ranges
        d_gam_e = float(np.hypot(Ee0[0, 2] - Ee1[0, 2], Ee0[1, 2] - Ee1[1, 2]))
        d_gam_p = float(np.hypot(Ep0[0, 2] - Ep1[0, 2], Ep0[1, 2] - Ep1[1, 2]))
        d_eps_e = abs(float(Ee0[2, 2] - Ee1[2, 2]))
        d_eps_p = abs(float(Ep0[2, 2] - Ep1[2, 2]))

        gsa = (
            (tau_max / tau_f) * (0.5 * d_gam_e)
            + (0.5 * d_gam_p)
            + (sig_n_max / sigma_f) * (0.5 * d_eps_e)
            + (0.5 * d_eps_p)
        )
        gsa = float(gsa)

        return CPPlaneResult(damage=gsa, metric=gsa)
