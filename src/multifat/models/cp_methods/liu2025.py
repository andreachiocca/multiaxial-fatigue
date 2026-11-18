"""LIU2025: Liu, Xie & Hao (2025) multiaxial HCF life model on the maximum shear stress plane.

X.-Y. Liu, J. Xie, W. Hao, J. Mater. Eng. Perform. 34 (2025) 24315-24322,
doi:10.1007/s11665-025-11163-z, Eqs. (21)-(27), on the plane of maximum shear stress amplitude::

    sqrt( (dsigma'_n/2)^2 + 3 (dtau'_max/2)^2 ) = sigma'_f (2 N_f)^b                      (27)
    dtau'_max/2 = (dtau_max/2) / (1 - (tau_max,m / tau_b)^2)          (Gerber, Eq. 25)
    dsigma'_n/2 = dsigma_n/2 + 0.414 sigma_n,m        for -1 <= R_n <= 0
                = sqrt(2) sigma_b (dsigma_n/2) / (sigma_b + dsigma_n/2 - sigma_n,m)   for R_n > 0
                                                       (Serensen, Eqs. 23, 24, 26)

R_n = (sigma_n,m - dsigma_n/2)/(sigma_n,m + dsigma_n/2); tau_max,m and sigma_n,m are the mean shear
and normal stress on the critical plane; sigma_b and tau_b the tensile and torsional ultimate
strengths; sigma'_f and b the uniaxial Basquin constants. Of the two planes of maximum shear
amplitude (Eq. 9), the one with the larger normal stress peak is the critical plane. The left-hand
side (MPa) is returned as damage parameter.

Repository choices:
  - shear amplitude and mean on the plane from the longest chord of the shear path (the paper's
    half range of the scalar shear on surface planes, Eqs. 10-11); ties resolved towards the
    larger damage, which on the source's planes is the larger normal stress peak;
  - 0.414 is sqrt(2) - 1, the value implied by Eqs. (23)-(24) (the paper prints it rounded);
  - a compressive mean normal stress (sigma_n,m < 0) uses the linear branch, continuous at
    sigma_n,m = 0, with the corrected amplitude clipped at zero; Eq. (26) does not cover
    -dsigma_n/2 < sigma_n,m < 0, and its hyperbolic branch would make a fully compressive cycle
    damaging;
  - tau_b: workbook column ``Tu``, else estimated as Su / sqrt(3) with a warning (the paper's
    Table 1 ratios tau_b/sigma_b are 0.78, 0.68 and 0.94, so the estimate is coarse); a plane
    with |tau_m| >= tau_b, where Eq. (25) is undefined, raises LoadingScopeError and the model is
    skipped for that material; so is a plane whose minimum normal stress reaches sigma_b, where
    the denominator of Eq. (26) vanishes (in load-scaling searches such trial loads count as
    infinite damage, see multifat.eval.endurance.scaled_damage);
  - LIU2025_ext is the largest damage parameter over all planes (repository variant); life from
    the shared master curve DP = A N^b instead of the Basquin constants of Eq. (27).

User guide: docs/method-guides/LIU2025.md.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult
from ..errors import LoadingScopeError

SERENSEN = np.sqrt(2.0) - 1.0  # (2 S_-1 - S_0)/S_0 with S_0 = sqrt(2) S_-1, Eqs. (23)-(24)


def serensen_amplitude(sigma_a, sigma_m, su):
    """Equivalent fully reversed normal amplitude of Eq. (26) (arrays allowed)."""
    sigma_a, sigma_m = np.asarray(sigma_a, float), np.asarray(sigma_m, float)
    hyperbolic = np.sqrt(2.0) * su * sigma_a / np.maximum(su + sigma_a - sigma_m, 1e-300)
    linear = np.maximum(sigma_a + SERENSEN * sigma_m, 0.0)
    return np.where(sigma_m > sigma_a, hyperbolic, linear)


def gerber_amplitude(tau_a, tau_m, tu):
    """Equivalent fully reversed shear amplitude of Eq. (25); infinite for |tau_m| >= tau_b."""
    tau_a, tau_m = np.asarray(tau_a, float), np.asarray(tau_m, float)
    denom = 1.0 - (tau_m / tu) ** 2
    with np.errstate(divide="ignore"):
        return np.where(denom > 0.0, tau_a / np.where(denom > 0.0, denom, 1.0), np.inf)


class Liu2025:
    """Liu-Xie-Hao (2025) equivalent stress on the maximum shear stress amplitude plane."""

    name = "LIU2025"
    np_shear_rule = "longest_chord"
    optional_params = {"Tu": "Su / sqrt(3) (warning)"}

    def required_params(self) -> set[str]:
        return {"Su"}

    def _strengths(self, params: Mapping[str, float]) -> tuple[float, float]:
        su = float(params["Su"])
        tu = params.get("Tu")
        if tu is None:
            warnings.warn(
                "LIU2025: torsional ultimate strength Tu not supplied; estimated as Su/sqrt(3)",
                RuntimeWarning,
                stacklevel=3,
            )
            tu = su / np.sqrt(3.0)
        tu = float(tu)
        if su <= 0.0 or tu <= 0.0:
            raise ValueError(f"LIU2025 needs Su > 0 and Tu > 0, got Su = {su}, Tu = {tu}")
        return su, tu

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]):
        damage, metric = self.evaluate_on_planes(
            S0r=np.asarray(S0r, float)[None],
            E0r=np.asarray(E0r, float)[None],
            S1r=np.asarray(S1r, float)[None],
            E1r=np.asarray(E1r, float)[None],
            params=params,
        )
        return CPPlaneResult(damage=float(damage[0]), metric=float(metric[0]))

    def evaluate_on_planes(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
        tau_max=None,
        strain_ranges=None,
    ):
        """All planes at once: Eq. (27) left-hand side and the shear stress amplitude."""
        check_planes(S0r, E0r, S1r, E1r)
        su, tu = self._strengths(params)
        tau_a = 0.5 * np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        tau_m = 0.5 * np.hypot(S0r[:, 0, 2] + S1r[:, 0, 2], S0r[:, 1, 2] + S1r[:, 1, 2])
        sigma_a = 0.5 * np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
        sigma_m = 0.5 * (S0r[:, 2, 2] + S1r[:, 2, 2])
        if np.any((sigma_m > sigma_a) & (su + sigma_a - sigma_m <= 0.0)):
            raise LoadingScopeError(
                f"LIU2025: a plane has a minimum normal stress >= sigma_b = {su:g} MPa; the "
                "Serensen correction, Eq. (26), is undefined"
            )
        sigma_eq = serensen_amplitude(sigma_a, sigma_m, su)
        if np.any(np.abs(tau_m) >= tu):
            raise LoadingScopeError(
                f"LIU2025: a plane carries a mean shear stress >= tau_b = {tu:g} MPa; the Gerber "
                "correction, Eq. (25), is undefined (supply Tu if it was estimated from Su)"
            )
        tau_eq = gerber_amplitude(tau_a, tau_m, tu)
        return np.sqrt(sigma_eq**2 + 3.0 * tau_eq**2), tau_a
