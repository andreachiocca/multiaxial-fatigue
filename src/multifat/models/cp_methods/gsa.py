"""GSA: Ince & Glinka (2014) generalized strain amplitude critical-plane parameter.

A. Ince, G. Glinka, Int. J. Fatigue 62 (2014) 34-41, doi:10.1016/j.ijfatigue.2013.10.007,
Eq. (5), the strain-amplitude form of the generalized strain energy (GSE, Eq. 1)::

    de*/2 = [ (tau_max/tau'_f) dgamma_e/2 + dgamma_p/2
              + (sigma_n,max/sigma'_f) de_e/2 + de_p/2 ]_max

maximised over the planes (the plane rule is the parameter itself).

Repository choices:
  - sigma_n,max is clipped at zero;
  - plastic strain = total strain - isotropic elastic strain of the stresses (E, nu); workbook
    strains are elastic, so the plastic terms vanish for imported histories;
  - missing sigma'_f is estimated as 1.5 Su (uniform material law, Zhao et al. 2020; a steel
    relation, with a warning); missing tau'_f = sigma'_f / sqrt(3) (the paper gives no value);
  - life from the shared master curve DP = A N^b instead of the strain-life curve.

User guide: docs/method-guides/GSA.md.
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, elastic_strain_planes, tau_max_planes

from ..base import CPPlaneResult
from ..utils.fatigue_estimation import ensure_axial_fatigue_props


def _elastic_strain_from_stress_iso(sig: np.ndarray, *, E: float, nu: float) -> np.ndarray:
    sig = np.asarray(sig, dtype=float)
    tr = float(np.trace(sig))
    I = np.eye(3)
    return ((1.0 + float(nu)) / float(E)) * sig - (float(nu) / float(E)) * tr * I


class GSA:
    name = "GSA"
    optional_params = {
        "sigma_f": "uniform material law from Su and E",
        "tau_f": "sigma_f / sqrt(3)",
        "Su": "only to estimate sigma_f",
    }

    def required_params(self) -> set[str]:
        # E, nu required for elastic/plastic split; Su is used if fatigue properties are missing.
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

        # Ensure fatigue properties
        sigma_f = params.get("sigma_f")
        tau_f = params.get("tau_f")
        if sigma_f is None:
            sigma_f = ensure_axial_fatigue_props(params).sigma_f
        if tau_f is None:
            tau_f = float(sigma_f) / np.sqrt(3.0)
        sigma_f, tau_f = float(sigma_f), float(tau_f)
        if sigma_f <= 0.0 or tau_f <= 0.0:
            raise ValueError(
                f"{self.name}: invalid fatigue strength coefficients sigma_f={sigma_f}, tau_f={tau_f}"
            )

        # elastic/plastic strain tensors
        Ee0 = _elastic_strain_from_stress_iso(S0r, E=E, nu=nu)
        Ee1 = _elastic_strain_from_stress_iso(S1r, E=E, nu=nu)
        Ep0 = E0r - Ee0
        Ep1 = E1r - Ee1

        # maxima
        sigma_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)
        tau_0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau_1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = float(params.get("_tau_max", max(tau_0, tau_1)))

        # elastic/plastic strain ranges
        d_gam_e = 2.0 * float(np.hypot(Ee0[0, 2] - Ee1[0, 2], Ee0[1, 2] - Ee1[1, 2]))
        d_gam_p = 2.0 * float(np.hypot(Ep0[0, 2] - Ep1[0, 2], Ep0[1, 2] - Ep1[1, 2]))
        d_eps_e = abs(float(Ee0[2, 2] - Ee1[2, 2]))
        d_eps_p = abs(float(Ep0[2, 2] - Ep1[2, 2]))
        if "_strain_ranges" in params:
            d_gam_e, d_gam_p, d_eps_e, d_eps_p = params["_strain_ranges"]

        gsa = (
            (tau_max / tau_f) * (0.5 * d_gam_e)
            + (0.5 * d_gam_p)
            + (sigma_n_max / sigma_f) * (0.5 * d_eps_e)
            + (0.5 * d_eps_p)
        )
        gsa = float(gsa)

        return CPPlaneResult(damage=gsa, metric=gsa)

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
        sigma_f = params.get("sigma_f")
        tau_f = params.get("tau_f")
        if sigma_f is None:
            sigma_f = ensure_axial_fatigue_props(params).sigma_f
        if tau_f is None:
            tau_f = float(sigma_f) / np.sqrt(3.0)
        sigma_f, tau_f = float(sigma_f), float(tau_f)
        if sigma_f <= 0.0 or tau_f <= 0.0:
            raise ValueError(
                f"{self.name}: invalid fatigue strength coefficients sigma_f={sigma_f}, tau_f={tau_f}"
            )
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
        gsa = (
            (tau_max_plane / tau_f) * (0.5 * d_gam_e)
            + (0.5 * d_gam_p)
            + (sigma_n_max / sigma_f) * (0.5 * d_eps_e)
            + (0.5 * d_eps_p)
        )
        return gsa, gsa
