"""Li et al. (2021) generalized strain-energy-density critical-plane parameter.

Li et al. propose (Eq.(3) / Eq.(25)):

  (Δτ_t/2 + |τ_m|) * (Δγ_max/2) + L * σ_{n,max} * (Δε_n/2) = f(N_f)

where the critical plane is defined as the *maximum shear strain energy density*
plane with the larger normal strain energy density.

The left-hand scalar is returned as the damage parameter; the main pipeline
calibrates it to fatigue life on the uniaxial subset.

Repository choices:
  - candidate plane normal is local z-axis
  - shear traction vector is (xz, yz) and its 2-norm is used
  - shear strain vector is (xz, yz) and its 2-norm is used
  - σ_{n,max} uses Szz, clipped to tension-only (>=0)
  - material constant L is read from params key "L_LI" (default 1.0)
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult
from ..utils.param_utils import get_param


class LI:
    name = "LI"
    fitted_constants = ("L_LI",)
    input_classes = ("elastic constants",)
    optional_params = {"L_LI": "1.0"}
    # Source defines the plane shear range over the cycle (longest chord), see SHEAR_RULES.
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        # L_LI is optional (default 1.0)
        return set()

    def evaluate_on_plane(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]
    ) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        L = float(get_param(params, "L_LI", default=1.0, context=self.name))

        # Normal stress max (tension only)
        sigma_n_max = max(float(S0r[2, 2]), float(S1r[2, 2]), 0.0)

        # Mean shear stress magnitude (Li uses |τ_m| and notes both signs are detrimental)
        tau_m = float(np.hypot(0.5 * (S0r[0, 2] + S1r[0, 2]), 0.5 * (S0r[1, 2] + S1r[1, 2])))

        # Shear stress range magnitude
        d_tau = float(np.hypot(S0r[0, 2] - S1r[0, 2], S0r[1, 2] - S1r[1, 2]))

        # Shear strain range magnitude ("Δγ_max" in the paper)
        d_gam = 2.0 * float(np.hypot(E0r[0, 2] - E1r[0, 2], E0r[1, 2] - E1r[1, 2]))

        # Normal strain range
        d_eps_n = abs(float(E0r[2, 2] - E1r[2, 2]))

        shear_term = (0.5 * d_tau + abs(tau_m)) * (0.5 * d_gam)
        normal_term = L * sigma_n_max * (0.5 * d_eps_n)
        damage = float(shear_term + normal_term)

        # Plane selection: prioritize the maximum shear strain energy density plane;
        # ties are resolved by damage in the case-level evaluator.
        metric = float(shear_term)
        return CPPlaneResult(damage=damage, metric=metric)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        L = float(get_param(params, "L_LI", default=1.0, context=self.name))
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], S1r[:, 2, 2]), 0.0)
        tau_m = np.hypot(0.5 * (S0r[:, 0, 2] + S1r[:, 0, 2]), 0.5 * (S0r[:, 1, 2] + S1r[:, 1, 2]))
        d_tau = np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        d_gam = 2.0 * np.hypot(E0r[:, 0, 2] - E1r[:, 0, 2], E0r[:, 1, 2] - E1r[:, 1, 2])
        d_eps_n = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        shear_term = (0.5 * d_tau + np.abs(tau_m)) * (0.5 * d_gam)
        return shear_term + L * sigma_n_max * (0.5 * d_eps_n), shear_term
