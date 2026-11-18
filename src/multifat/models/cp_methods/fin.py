"""Findley (1959), J. Eng. Ind., doi:10.1115/1.4008327: stress-based critical-plane criterion.

Eqs. (2), (8)-(10), p. 302, on the plane of orientation theta::

    tau_thetaA + k sigma_thetaM = f        maximised over the planes

with tau_thetaA the shear stress amplitude, sigma_thetaM the maximum normal stress on the plane
and k a material constant. Fully reversed torsion and bending give f = t_-1 sqrt(1 + k^2)
(Eq. 19) and f = f_-1 (k + sqrt(1 + k^2)) / 2 (Eq. 24).

The plane damage parameter is DP = d_tau / 2 + k_FI sigma_max (MPa), with the selection metric
equal to DP, so FIN and FIN_ext coincide.

Repository choices:
  - sigma_max keeps its sign (compression lowers DP); no clipping;
  - out-of-phase shear amplitude: root-sum-square of the shear-ellipse semi-axes (the default
    ``rss`` rule of ``multifat.eval.evaluator.SHEAR_RULES``); the paper treats in-phase loading;
  - life from the shared master curve DP = A N^b, and k_FI calibrated by minimising the life
    scatter of tension and torsion failures (``multifat.cli.calibrate``).

User guide: docs/method-guides/FIN.md.
"""

from __future__ import annotations
from typing import Mapping, Tuple
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, pypow2
from ..base import CPPlaneResult


def _findley_plane_value(
    *,
    S0r: np.ndarray,
    S1r: np.ndarray,
    k_FI: float,
) -> Tuple[float, float]:
    """Findley DP on one candidate plane.

    Plane normal = local z-axis; normal stress = Szz (signed maximum, including
    compression); shear stress range from Sxz and Syz.

    Returns: (DP, metric) where metric = DP.
    """
    S0r = _as33("S0r", S0r)
    S1r = _as33("S1r", S1r)
    k_FI = float(k_FI)

    if k_FI < 0.0:
        raise ValueError(f"k_FI must be >= 0. Got k_FI={k_FI}")

    # Signed maximum normal stress, as used in the Findley expression.
    sigma_z0 = float(S0r[2, 2])
    sigma_z1 = float(S1r[2, 2])
    sigma_max = max(sigma_z0, sigma_z1)

    # tau_a = Sxz, tau_b = Syz
    d_tau_xz = float(S0r[0, 2] - S1r[0, 2])
    d_tau_yz = float(S0r[1, 2] - S1r[1, 2])
    d_tau = float(np.sqrt(d_tau_xz**2 + d_tau_yz**2))

    DP = float(d_tau / 2.0 + k_FI * sigma_max)
    metric = DP
    return DP, metric


class Findley:
    name = "FIN"
    fitted_constants = ("k_FI",)
    input_classes = ()

    def required_params(self) -> set[str]:
        return {"k_FI"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        # Findley uses stresses only (E0r/E1r are ignored here, but kept for a uniform interface).
        damage, metric = _findley_plane_value(S0r=S0r, S1r=S1r, k_FI=params["k_FI"])
        return CPPlaneResult(damage=float(damage), metric=float(metric))

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        k_FI = float(params["k_FI"])
        if k_FI < 0.0:
            raise ValueError(f"k_FI must be >= 0. Got k_FI={k_FI}")
        sigma_max = np.maximum(S0r[:, 2, 2], S1r[:, 2, 2])
        d_tau_xz = S0r[:, 0, 2] - S1r[:, 0, 2]
        d_tau_yz = S0r[:, 1, 2] - S1r[:, 1, 2]
        d_tau = np.sqrt(pypow2(d_tau_xz) + pypow2(d_tau_yz))
        DP = d_tau / 2.0 + k_FI * sigma_max
        return DP, DP
