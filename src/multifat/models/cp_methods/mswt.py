"""Modified Smith-Watson-Topper (MSWT) critical-plane damage parameter (Zhu et al. 2018).

Range form of Jiang's plastic-work increment dY = a sigma d eps^p + (1-a)/2 tau d gamma^p
(Jiang 2000, Eq. 6; Jiang et al. 2007, Eq. 4) as written in Zhu et al., Int. J. Fatigue
113 (2018), Eq. (13), p. 8::

    DP_MSWT = 2 * a * Delta_eps * <sigma_max> + (1 - a) / 2 * Delta_gamma * Delta_tau

Assuming a plastic Poisson ratio of 0.5 recovers Jiang's crack-mode bounds (shear for
a <= 0.375, mixed for 0.375-0.5, tensile for a >= 0.5). Up to 2026-10-18 the
implementation evaluated sigma_max d_eps/2 + (1-a) d_gam d_tau/4, which gave the shear
term twice its intended weight (finding U2); the former a carries over as a' = a/(2 - a)
and yields identical lives.

Repository choices:
  - Implements the scalar range form of Zhu et al. (2018) Eq. (13) on elastic workbook
    strains, not Jiang's incremental cyclic plasticity model.
  - Candidate plane normal is the local z-axis; shear uses (xz, yz) and normal uses zz.
  - Peak normal stress sigma_max is Macaulay-bracketed at zero (<x> = max(x, 0)).
  - Non-proportional plane shear range uses the default root-sum-square (rss) rule.
  - Critical plane maximises the full DP, so standard and ``_ext`` variants coincide.
  - ``a_MSWT`` is optional, validated in [0, 1], and defaults to 0.5 (with a warning)
    when omitted; ``multifat-pipeline`` calibrates it by scatter minimisation.
  - Life prediction uses the shared power-law master curve.
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult
from ..utils.param_utils import get_param


class MSWT:
    name = "MSWT"
    optional_params = {"a_MSWT": "0.5"}

    def required_params(self) -> set[str]:
        # a_MSWT is not mandatory (falls back to 0.5)
        return set()

    def evaluate_on_plane(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]
    ) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        a = float(get_param(params, "a_MSWT", default=0.5, context=self.name))
        if not np.isfinite(a) or not 0.0 <= a <= 1.0:
            raise ValueError("a_MSWT must be finite and lie within [0, 1]")

        # peak normal stress on the plane, MacCauley-bracketed
        sigma_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)

        # range of the normal strain
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # range of the shear strain
        d_gam_xz = float(E0r[0, 2] - E1r[0, 2])
        d_gam_yz = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_xz, d_gam_yz))

        # range of the shear stress
        d_tau_xz = float(S0r[0, 2] - S1r[0, 2])
        d_tau_yz = float(S0r[1, 2] - S1r[1, 2])
        d_tau = float(np.hypot(d_tau_xz, d_tau_yz))

        dY_normal = 2.0 * a * d_eps * sigma_max
        dY_shear = 0.5 * (1.0 - a) * d_gam * d_tau
        DP = float(dY_normal + dY_shear)

        return CPPlaneResult(damage=DP, metric=DP)

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
        """Handle every plane in one shot; matches ``evaluate_on_plane`` exactly per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        a = float(get_param(params, "a_MSWT", default=0.5, context=self.name))
        if not np.isfinite(a) or not 0.0 <= a <= 1.0:
            raise ValueError("a_MSWT must be finite and lie within [0, 1]")
        sigma_max = np.maximum(np.maximum(S0r[:, 2, 2], S1r[:, 2, 2]), 0.0)
        d_eps = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        d_gam = 2.0 * np.hypot(E0r[:, 0, 2] - E1r[:, 0, 2], E0r[:, 1, 2] - E1r[:, 1, 2])
        d_tau = np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        DP = 2.0 * a * d_eps * sigma_max + 0.5 * (1.0 - a) * d_gam * d_tau
        return DP, DP
