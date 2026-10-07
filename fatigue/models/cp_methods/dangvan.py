from __future__ import annotations
from typing import Mapping
from ..base import CPPlaneResult
import numpy as np
 
def _sigma_h(S: np.ndarray) -> float:
    """Hydrostatic stress = tr(S)/3 for a single 3x3 tensor."""
    S = np.asarray(S, dtype=float)
    if S.shape != (3, 3):
        raise ValueError(f"Stress must be (3,3), got {S.shape}")
    return float((S[0, 0] + S[1, 1] + S[2, 2]) / 3.0)
 
class DangVan:
    name = "DANGVAN"
 
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
 
        from fatigue.models.references import ensure_enabled
        ensure_enabled("DANGVAN")
        sigm1 = float(params["Sigm1"])
        taum1 = float(params["Taum1"])
 
        if sigm1 <= 0.0 or taum1 <= 0.0:
            raise ValueError(f"Sigm1 and Taum1 must be > 0. Got Sigm1={sigm1}, Taum1={taum1}")
 
        # Dang Van parameters:
        #   b = tau_-1
        #   1/2*sigma_-1 + a*(1/3*sigma_-1) = b  =>  a = 3*(b/sigma_-1 - 1/2)
        b = taum1
        a = 3.0 * (b / sigm1 - 0.5)
 
#        tau_a = Sxz, tau_b = Syz
        dtau_a = float(S0r[0, 2] - S1r[0, 2])
        dtau_b = float(S0r[1, 2] - S1r[1, 2])
        DeltaTau = np.sqrt(dtau_a**2 + dtau_b**2)
        # Evaluate at the two extrema (max/min). For proportional sinusoidal loading,
        # max over time is captured by these endpoints.
        dv0 = DeltaTau/2 + a * _sigma_h(S0r)
        dv1 = DeltaTau/2 + a * _sigma_h(S1r)
 
        dv_max = float(max(dv0, dv1))
 
        # Optional clamp to avoid negative values causing your fitter to drop points
        dv_max = max(0.0, dv_max)
 
        # You have two reasonable choices:
        #   (A) return dv_max (units: MPa)  -> like a stress-equivalent parameter
        #   (B) return dv_max/b (dimensionless) -> "utilization" vs fatigue limit
        damage = dv_max          # choose (A)
        metric = DeltaTau          # Dang Van plane selection is the same maximization
 
        return CPPlaneResult(damage=float(damage), metric=float(metric))