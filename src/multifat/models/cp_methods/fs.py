"""Fatemi & Socie (1988), FFEMS 11, doi:10.1111/j.1460-2695.1988.tb01169.x: critical-plane strain
criterion.

Eq. (8), on the plane of maximum shear strain amplitude::

    (Delta gamma_max / 2) (1 + n sigma_n^max / sigma_y) = constant

with Delta gamma_max / 2 the engineering shear strain amplitude, sigma_n^max the maximum normal
stress on that plane and n (``k_FS``) a material constant.

The plane damage parameter is DP = d_gam (1 + k_FS sigma_max / Sy): the range of the tensorial
shear strain on the plane equals the engineering shear amplitude, so no factor two is applied.
The selection metric is d_gam; FS_ext takes the largest DP over the planes.

Repository choices:
  - sigma_max is clipped at zero (compression never lowers DP);
  - the workbook (monotonic) yield strength Sy replaces the cyclic yield strength;
  - out-of-phase shear range: longest chord of the shear-strain path (``longest_chord``);
  - strains are linear-elastic; life from the shared master curve DP = A N^b instead of Eq. (9).

User guide: docs/method-guides/FS.md.
"""

from __future__ import annotations
from typing import Mapping, Tuple
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, pypow2
from ..base import CPPlaneResult


def _fatemi_socie_plane_value(
    *,
    S0r: np.ndarray,
    E0r: np.ndarray,
    S1r: np.ndarray,
    E1r: np.ndarray,
    k_FS: float,
    sigma_y: float,
) -> Tuple[float, float]:
    """Fatemi-Socie DP on one candidate plane.

    Plane normal = local z-axis; normal stress = Szz (clipped at 0 when taking the
    maximum); shear strain range from Exz and Eyz.

    Returns
    -------
    DP : float
        FS value for this plane.
    metric : float
        Plane-selection metric (d_gam).
    """

    S0r = _as33("S0r", S0r)
    E0r = _as33("E0r", E0r)
    S1r = _as33("S1r", S1r)
    E1r = _as33("E1r", E1r)

    k_FS = float(k_FS)
    sigma_y = float(sigma_y)

    if sigma_y <= 0.0:
        raise ValueError(f"Sy must be > 0. Got Sy={sigma_y}")
    if k_FS < 0.0:
        raise ValueError(f"k_FS must be >= 0. Got k_FS={k_FS}")

    # sigma_nn = Szz (clipped at 0)
    sigma_z0 = max(float(S0r[2, 2]), 0.0)
    sigma_z1 = max(float(S1r[2, 2]), 0.0)
    sigma_max = max(sigma_z0, sigma_z1)

    # Tensor shear RANGE equals engineering shear AMPLITUDE (gamma = 2*epsilon).
    d_gam_xz = float(E0r[0, 2] - E1r[0, 2])
    d_gam_yz = float(E0r[1, 2] - E1r[1, 2])
    d_gam = float(np.sqrt(d_gam_xz**2 + d_gam_yz**2))

    DP = float(d_gam * (1.0 + k_FS * (sigma_max / sigma_y)))

    # Plane-selection metric (d_gam).
    metric = d_gam

    return DP, metric


class FatemiSocie:
    name = "FS"
    fitted_constants = ("k_FS",)
    input_classes = ("elastic constants", "static strength")
    # Source defines the plane shear range over the cycle (longest chord), see SHEAR_RULES.
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        return {"Sy", "k_FS"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        damage, metric = _fatemi_socie_plane_value(
            S0r=S0r,
            E0r=E0r,
            S1r=S1r,
            E1r=E1r,
            k_FS=params["k_FS"],
            sigma_y=params["Sy"],
        )
        return CPPlaneResult(damage=float(damage), metric=float(metric))

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        k_FS, sigma_y = float(params["k_FS"]), float(params["Sy"])
        if sigma_y <= 0.0:
            raise ValueError(f"Sy must be > 0. Got Sy={sigma_y}")
        if k_FS < 0.0:
            raise ValueError(f"k_FS must be >= 0. Got k_FS={k_FS}")
        sigma_max = np.maximum(np.maximum(S0r[:, 2, 2], 0.0), np.maximum(S1r[:, 2, 2], 0.0))
        d_gam_xz = E0r[:, 0, 2] - E1r[:, 0, 2]
        d_gam_yz = E0r[:, 1, 2] - E1r[:, 1, 2]
        d_gam = np.sqrt(pypow2(d_gam_xz) + pypow2(d_gam_yz))
        return d_gam * (1.0 + k_FS * (sigma_max / sigma_y)), d_gam
