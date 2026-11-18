"""Carpinteri & Spagnoli (2001) critical-plane equivalent stress (disabled).

Published in *International Journal of Fatigue* 23 (2001), pp. 135-145,
doi:10.1016/S0142-1123(00)00075-X (roadmap item I_CARSPA)::

    sigma_eq = f_-1 * sqrt((sigma_n,max / f_-1)^2 + (d_tau / (2 * t_-1))^2)

where ``sigma_n,max`` and ``d_tau`` are the maximum normal stress and shear stress
range on the paper's averaged principal-stress critical plane.

Repository choices:
- Disabled at evaluation time via ``ensure_enabled("CARSPA")`` because the published
  principal-direction averaging rule to orient the critical plane is not implemented;
  scanning planes by maximum shear range or maximum damage does not match the criterion.
"""

from __future__ import annotations
from typing import Mapping
from ..base import CPPlaneResult
import numpy as np


class CarSpa:
    name = "CARSPA"

    def required_params(self) -> set[str]:
        # Sigm1: fully reversed tensile fatigue limit (sigma_-1)
        # Taum1: fully reversed torsional fatigue limit (tau_-1)
        return {"Sigm1", "Taum1"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,  # S0r: max state (loadstep 1) in rotated frame
        S1r,
        E1r,  # S1r: min state (loadstep 2) in rotated frame
        params: Mapping[str, float],
    ) -> CPPlaneResult:

        from multifat.models.references import ensure_enabled

        ensure_enabled("CARSPA")
        f_m1 = float(params["Sigm1"])
        t_m1 = float(params["Taum1"])

        # Car-Spa equivalent stress:
        sigma_n_max = max(S0r[2, 2], S1r[2, 2])  # Maximum normal stress on the CP
        #        tau_a = Sxz, tau_b = Syz
        d_tau_xz = float(S0r[0, 2] - S1r[0, 2])
        d_tau_yz = float(S0r[1, 2] - S1r[1, 2])
        d_tau = np.sqrt(d_tau_xz**2 + d_tau_yz**2)  # Shear stress range on the CP
        # Evaluate at the two extrema (max/min). For proportional sinusoidal loading,
        # max over time is captured by these endpoints.
        DP = f_m1 * np.sqrt((sigma_n_max / f_m1) ** 2 + (d_tau / 2 / t_m1) ** 2)

        # Stress-equivalent parameter in MPa (option A); the plane metric is the shear range.
        damage = DP
        metric = d_tau

        return CPPlaneResult(damage=float(damage), metric=float(metric))
