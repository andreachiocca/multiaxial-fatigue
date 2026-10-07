# fatigue/models/cp_methods/li.py
from __future__ import annotations

"""Li et al. (2021) generalized strain-energy-density critical-plane parameter.

Li et al. propose (Eq.(3) / Eq.(25)):

  (Δτ_t/2 + |τ_m|) * (Δγ_max/2) + L * σ_{n,max} * (Δε_n/2) = f(N_f)

where the critical plane is defined as the *maximum shear strain energy density*
plane with the larger normal strain energy density.

In this repository we implement the left-hand scalar as a damage parameter and
let the main pipeline calibrate it to fatigue life using the uniaxial subset.

Implementation notes (repository conventions):
  - candidate plane normal is local z-axis
  - shear traction vector is (xz, yz) and we use its 2-norm
  - shear strain vector is (xz, yz) and we use its 2-norm
  - σ_{n,max} uses Szz, clipped to tension-only (>=0)
  - material constant L is read from params key "L_LI" (default 1.0)
"""

from typing import Mapping
import numpy as np

from ..base import CPPlaneResult
from ..utils.param_utils import get_param


def _as33(name: str, A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
    return A


class LI:
    name = "LI"

    def required_params(self) -> set[str]:
        # L_LI is optional (default 1.0)
        return set()

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        L = float(get_param(params, "L_LI", default=1.0, context=self.name))

        # Normal stress max (tension only)
        sig_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)

        # Shear stress magnitudes at the two load steps
        tau0 = float(np.hypot(S0r[0, 2], S0r[1, 2]))
        tau1 = float(np.hypot(S1r[0, 2], S1r[1, 2]))

        # Mean shear stress magnitude (Li uses |τ_m| and notes both signs are detrimental)
        tau_m = 0.5 * (tau0 + tau1)

        # Shear stress range magnitude
        d_tau = float(np.hypot(S0r[0, 2] - S1r[0, 2], S0r[1, 2] - S1r[1, 2]))

        # Shear strain range magnitude ("Δγ_max" in the paper)
        d_gam = float(np.hypot(E0r[0, 2] - E1r[0, 2], E0r[1, 2] - E1r[1, 2]))

        # Normal strain range
        d_eps_n = abs(float(E0r[2, 2] - E1r[2, 2]))

        shear_term = (0.5 * d_tau + abs(tau_m)) * (0.5 * d_gam)
        normal_term = L * sig_n_max * (0.5 * d_eps_n)
        damage = float(shear_term + normal_term)

        # Plane selection: prioritize the maximum shear strain energy density plane;
        # use a tiny tie-break towards larger normal energy.
        metric = float(shear_term + 1e-12 * normal_term)
        return CPPlaneResult(damage=damage, metric=metric)
