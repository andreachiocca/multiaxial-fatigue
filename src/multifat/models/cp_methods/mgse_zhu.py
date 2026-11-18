"""Zhu et al. (2018) modified generalized strain energy (MGSE_ZHU) criterion.

Zhu et al., Int. J. Fatigue 113 (2018), Eqs. (7)-(12), pp. 7-8, define::

    MGSE = tau_max * (Delta_gamma / 2) + k_MGSE * sigma_{n,max} * (Delta_eps_n / 2)

where tau_max is the maximum shear stress on the critical plane, Delta_gamma / 2
and Delta_eps_n / 2 are the engineering shear and normal strain amplitudes,
sigma_{n,max} is the peak normal stress, and k_MGSE weights the normal term.

Repository choices:
  - Candidate plane normal is the local z-axis; shear uses (xz, yz) and normal uses zz.
  - Peak normal stress sigma_{n,max} is clipped at zero (tension-only).
  - Non-proportional plane shear range uses the default root-sum-square (rss) rule;
    tau_max comes from the 720-sample harmonic shear maximum via ``params['_tau_max']``
    (falling back to endpoint max when called outside ``prepare_planes``).
  - Critical plane maximises the full MGSE parameter, so standard and ``_ext`` coincide.
  - ``k_MGSE`` is optional and defaults to 1.0 (with a warning) when omitted; in
    ``multifat-pipeline`` it is derived from Zhu (2018) Eq. (12) over N = 1e4..2e6
    when Basquin fits are identifiable, else fitted by scatter minimisation.
  - Life prediction uses the shared elastic master curve instead of a strain-life curve.
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, tau_max_planes as _tau_max_planes

from ..base import CPPlaneResult
from ..utils.param_utils import get_param


class MGSE_Zhu:
    name = "MGSE_ZHU"
    optional_params = {"k_MGSE": "1.0"}

    def required_params(self) -> set[str]:
        # k_MGSE is optional (default 1.0)
        return set()

    def evaluate_on_plane(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]
    ) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        k_mgse = get_param(params, "k_MGSE", default=1.0, context=self.name)

        # Normal stress max (tensile only)
        sigma_z0 = max(float(S0r[2, 2]), 0.0)
        sigma_z1 = max(float(S1r[2, 2]), 0.0)
        sigma_n_max = max(sigma_z0, sigma_z1)

        # Shear stress magnitude at max/min load step
        tau_0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau_1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))
        tau_max = float(params.get("_tau_max", max(tau_0, tau_1)))

        # Normal strain range on plane
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range on plane
        d_gam_xz = float(E0r[0, 2] - E1r[0, 2])
        d_gam_yz = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_xz, d_gam_yz))

        damage = float(tau_max * (0.5 * d_gam) + float(k_mgse) * sigma_n_max * (0.5 * d_eps))
        metric = damage
        return CPPlaneResult(damage=damage, metric=metric)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        k_mgse = get_param(params, "k_MGSE", default=1.0, context=self.name)
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], 0.0), np.maximum(S1r[:, 2, 2], 0.0))
        tau_max_plane = _tau_max_planes(S0r, S1r, params, tau_max)
        d_eps = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        d_gam = 2.0 * np.hypot(E0r[:, 0, 2] - E1r[:, 0, 2], E0r[:, 1, 2] - E1r[:, 1, 2])
        damage = tau_max_plane * (0.5 * d_gam) + float(k_mgse) * sigma_n_max * (0.5 * d_eps)
        return damage, damage
