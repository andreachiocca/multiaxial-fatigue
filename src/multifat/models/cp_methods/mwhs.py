"""MWHS: Wu, Li, Fang & Song (2017) modified Wu-Hu-Song critical-plane strain parameter (N28).

Z. R. Wu, X. Li, L. Fang, Y. D. Song, Proc. IMechE Part C: J. Mech. Eng. Sci. (2017),
doi:10.1177/0954406217708214, Eq. (19), on the maximum shear strain plane::

    dgamma_max/2 + 1/2 (1 + sigma_n,max/sigma_y) (sigma_n,max de_n / E)^0.5
        = tau'_f/G (2 N_f)^b0 + gamma'_f (2 N_f)^c0                                      (19)

dgamma_max/2 is the maximum (engineering) shear strain amplitude, sigma_n,max and de_n the maximum
normal stress and the normal strain range on that plane, sigma_y the yield strength. Eq. (19)
replaces the life-dependent weight k of the WHS parameter, Eqs. (13) and (18), by the stress-
dependent factor (1 + sigma_n,max/sigma_y)/2, so that the left-hand side is related to life by the
torsional strain-life curve alone.

Repository choices:
  - out-of-phase shear range: longest chord of the shear-strain path (``longest_chord``); with it
    the maximum shear strain of Fig. 14(a) is reproduced at every phase angle (rss is not);
  - sigma_n,max is the maximum normal stress over the cycle, clipped at zero (the square root
    needs a non-negative product; a compressive plane carries no normal-strain term);
  - ``MWHS`` takes the plane of largest shear strain range (ties: larger DP); ``MWHS_ext`` the
    largest DP over all planes (repository variant);
  - workbook (monotonic) Sy for sigma_y; elastic workbook strains (the source fits measured
    elastic-plastic strains);
  - life: the shared power-law master curve DP = A N^b fitted, like the right-hand side of
    Eq. (19), on the fully reversed pure torsion failures (``design_loading = "torsion"``); the
    uniaxial tests are predictions. The plastic term of the torsional Coffin-Manson curve is
    not represented by the single power law.

User guide: docs/method-guides/MWHS.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult


class MWHS:
    """Modified Wu-Hu-Song parameter on the maximum shear strain plane."""

    name = "MWHS"
    np_shear_rule = "longest_chord"
    # Eq. (19) relates the parameter to the torsional strain-life curve, whose constants come from
    # torsion tests only: the master curve is fitted on the pure torsion R = -1 failures.
    design_loading = "torsion"

    def required_params(self) -> set[str]:
        return {"Sy", "E"}

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
        """All planes at once: Eq. (19) left-hand side and the engineering shear strain range."""
        check_planes(S0r, E0r, S1r, E1r)
        sigma_y, E = float(params["Sy"]), float(params["E"])
        if sigma_y <= 0.0 or E <= 0.0:
            raise ValueError(f"MWHS needs Sy > 0 and E > 0, got Sy = {sigma_y}, E = {E}")
        d_gam = 2.0 * np.hypot(E0r[:, 0, 2] - E1r[:, 0, 2], E0r[:, 1, 2] - E1r[:, 1, 2])
        d_eps = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], S1r[:, 2, 2]), 0.0)
        normal_term = 0.5 * (1.0 + sigma_n_max / sigma_y) * np.sqrt(sigma_n_max * d_eps / E)
        return 0.5 * d_gam + normal_term, d_gam
