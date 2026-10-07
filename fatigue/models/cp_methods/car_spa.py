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
        self, *,
        S0r, E0r,   # S0r: max state (loadstep 1) in rotated frame
        S1r, E1r,   # S1r: min state (loadstep 2) in rotated frame
        params: Mapping[str, float]
    ) -> CPPlaneResult:
 
        sigm1 = float(params["Sigm1"])
        taum1 = float(params["Taum1"])
 
        # Car-Spa equivalent stress:
        sn_max = max(S0r[2,2],S1r[2,2]) # Maximum normal stress on the CP
#        tau_a = Sxz, tau_b = Syz
        dtau_a = float(S0r[0, 2] - S1r[0, 2])
        dtau_b = float(S0r[1, 2] - S1r[1, 2])
        DeltaTau = np.sqrt(dtau_a**2 + dtau_b**2)   # Shear stress range on the CP
        # Evaluate at the two extrema (max/min). For proportional sinusoidal loading,
        # max over time is captured by these endpoints.
        cs_max = sigm1*np.sqrt((sn_max/sigm1)**2+(DeltaTau/2/taum1)**2)
 
        # Optional clamp to avoid negative values causing your fitter to drop points
        #cs_max = max(0.0, dv_max)
 
        # You have two reasonable choices:
        #   (A) return dv_max (units: MPa)  -> like a stress-equivalent parameter
        #   (B) return dv_max/b (dimensionless) -> "utilization" vs fatigue limit
        damage = cs_max          # choose (A)
        metric = DeltaTau          # no metric in Car-Spa model

        return CPPlaneResult(damage=float(damage), metric=float(metric))