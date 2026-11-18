"""MKBM: Li, Zhang, Sun & Li (2011) modified Kandil-Brown-Miller critical-plane strain parameter.

J. Li et al., Int. J. Fatigue 33 (2011) 90-101, doi:10.1016/j.ijfatigue.2010.07.003, Eq. (27),
on the plane of maximum shear strain range (Eq. 19 for general paths)::

    de*_eq/2 = dgamma_max/2 + (1 + sigma_n,max / sigma_y) de_n/2

The left-hand equivalent strain amplitude is returned (dgamma: engineering shear strain range).
``MKBM`` evaluates it on the plane of largest dgamma; ``MKBM_ext`` is the maximum over all planes.

Repository choices:
  - sigma_n,max is clipped at zero (compression never lowers the factor below 1);
  - the workbook (monotonic) yield strength Sy replaces sigma_y;
  - out-of-phase shear range: longest chord of the shear-strain path (``longest_chord``, Eq. 19);
  - strains are linear-elastic; life from the shared master curve DP = A N^b instead of the
    strain-life curve.

User guide: docs/method-guides/MKBM.md.
"""

from __future__ import annotations

from typing import Mapping
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult


class MKBM:
    name = "MKBM"
    # Source defines the plane shear range over the cycle (longest chord), see SHEAR_RULES.
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        # Requires yield strength (σ_y)
        return {"Sy"}

    def evaluate_on_plane(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]
    ) -> CPPlaneResult:
        S0r = _as33("S0r", S0r)
        S1r = _as33("S1r", S1r)
        E0r = _as33("E0r", E0r)
        E1r = _as33("E1r", E1r)

        sigma_y = float(params["Sy"])
        if sigma_y <= 0:
            raise ValueError(f"Sy must be > 0. Got Sy={sigma_y}")

        # Normal stress max (tensile only)
        sigma_z0 = max(float(S0r[2, 2]), 0.0)
        sigma_z1 = max(float(S1r[2, 2]), 0.0)
        sigma_n_max = max(sigma_z0, sigma_z1)

        # Normal strain range
        d_eps = abs(float(E0r[2, 2] - E1r[2, 2]))

        # Shear strain range
        d_gam_xz = float(E0r[0, 2] - E1r[0, 2])
        d_gam_yz = float(E0r[1, 2] - E1r[1, 2])
        d_gam = 2.0 * float(np.hypot(d_gam_xz, d_gam_yz))

        damage = float(0.5 * d_gam + (1.0 + sigma_n_max / sigma_y) * (0.5 * d_eps))
        metric = float(d_gam)  # maximum shear strain plane
        return CPPlaneResult(damage=damage, metric=metric)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        sigma_y = float(params["Sy"])
        if sigma_y <= 0:
            raise ValueError(f"Sy must be > 0. Got Sy={sigma_y}")
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], 0.0), np.maximum(S1r[:, 2, 2], 0.0))
        d_eps = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        d_gam = 2.0 * np.hypot(E0r[:, 0, 2] - E1r[:, 0, 2], E0r[:, 1, 2] - E1r[:, 1, 2])
        return 0.5 * d_gam + (1.0 + sigma_n_max / sigma_y) * (0.5 * d_eps), d_gam
