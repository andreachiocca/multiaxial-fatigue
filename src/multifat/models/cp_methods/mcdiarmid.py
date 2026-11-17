"""McDiarmid (1991), FFEMS 14, doi:10.1111/j.1460-2695.1991.tb00673.x: critical-plane criterion.

Eqs. (7)-(8), p. 437, on the plane of maximum shear stress amplitude::

    tau_a / t_A,B + sigma_n,max / (2 sigma_u) = 1,        0 <= sigma_n,max <= sigma_u

with sigma_u the ultimate tensile strength and t_A,B the shear fatigue strength of case-A or
case-B cracks. The plane damage parameter is DP = d_tau / 2 + t_-1 / (2 Su) sigma_max (MPa),
equal to t_-1 on the limit; the selection metric is d_tau and MCD_ext takes the largest DP.

Repository choices:
  - case A only, with t_A = t_-1 (Taum1);
  - sigma_max is clipped at zero and the upper validity bound is not enforced;
  - the out-of-phase plane association of the source (plane of maximum principal-stress range)
    is not implemented; the out-of-phase shear amplitude is half the longest chord of the shear
    path (``longest_chord``, the circumscribed circle of Papadopoulos et al. 1997, whose Table 1
    it reproduces; maintainer decision, 17 Nov 2026; rss before);
  - life from the shared master curve DP = A N^b.

User guide: docs/method-guides/MCD.md.
"""

from __future__ import annotations

from typing import Mapping, Tuple
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, pypow2

from ..base import CPPlaneResult


def _mcdiarmid_plane_value(
    *,
    S0r: np.ndarray,
    S1r: np.ndarray,
    su: float,
    t_m1: float,
) -> Tuple[float, float]:
    """
    McDiarmid criterion (common HCF form) on one candidate plane.

    Critical plane: maximum shear stress amplitude/range (MSSR-style).

    A typical normalized statement is:
      tau_a / Taum1 + sigma_n,max / (2*Su) = 1
    Rearranged to a "damage-like" scalar (monotone with severity):
      tau_a + k * sigma_n,max
      where k = Taum1 / (2*Su)

    Conventions:
      - candidate plane normal is local z-axis
      - shear stress components on plane: Sxz and Syz
      - normal stress: Szz (clipped at 0 when taking maximum)
      - two load steps: (0=max), (1=min)
    """
    S0r = _as33("S0r", S0r)
    S1r = _as33("S1r", S1r)

    su = float(su)
    t_m1 = float(t_m1)

    if su <= 0.0:
        raise ValueError(f"Su (ultimate tensile strength) must be > 0. Got {su}")
    if t_m1 <= 0.0:
        raise ValueError(f"Taum1 (torsional fatigue limit) must be > 0. Got {t_m1}")

    # k from existing material properties (no new Excel fields)
    k = t_m1 / (2.0 * su)

    # Shear stress RANGE vector components on plane (xz, yz)
    d_tau_xz = float(S0r[0, 2] - S1r[0, 2])
    d_tau_yz = float(S0r[1, 2] - S1r[1, 2])
    d_tau = float(np.sqrt(d_tau_xz**2 + d_tau_yz**2))

    # shear stress amplitude
    tau_a = 0.5 * d_tau

    # max normal stress on plane (clip compression)
    sigma_z0 = max(float(S0r[2, 2]), 0.0)
    sigma_z1 = max(float(S1r[2, 2]), 0.0)
    sigma_n_max = max(sigma_z0, sigma_z1)

    damage = float(tau_a + k * sigma_n_max)

    # MSSR metric: maximize shear stress range (or amplitude)
    metric = d_tau
    return damage, metric


class McDiarmid:
    """
    McDiarmid critical plane method.
    Computes k internally from:
      k = Taum1 / (2*Su)
    """

    name = "MCD"
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        return {"Su", "Taum1"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,  # unused
        S1r,
        E1r,  # unused
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        if params.get("Su", None) is None:
            raise KeyError("MCD requires params['Su'] (Tensile_strength_MPa)")
        if params.get("Taum1", None) is None:
            raise KeyError("MCD requires params['Taum1'] (Torsional_fatigue_limit_MPa)")

        damage, metric = _mcdiarmid_plane_value(
            S0r=S0r,
            S1r=S1r,
            su=float(params["Su"]),
            t_m1=float(params["Taum1"]),
        )
        return CPPlaneResult(damage=damage, metric=metric)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        if params.get("Su", None) is None:
            raise KeyError("MCD requires params['Su'] (Tensile_strength_MPa)")
        if params.get("Taum1", None) is None:
            raise KeyError("MCD requires params['Taum1'] (Torsional_fatigue_limit_MPa)")
        su, t_m1 = float(params["Su"]), float(params["Taum1"])
        if su <= 0.0:
            raise ValueError(f"Su (ultimate tensile strength) must be > 0. Got {su}")
        if t_m1 <= 0.0:
            raise ValueError(f"Taum1 (torsional fatigue limit) must be > 0. Got {t_m1}")
        k = t_m1 / (2.0 * su)
        d_tau_xz = S0r[:, 0, 2] - S1r[:, 0, 2]
        d_tau_yz = S0r[:, 1, 2] - S1r[:, 1, 2]
        d_tau = np.sqrt(pypow2(d_tau_xz) + pypow2(d_tau_yz))
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], 0.0), np.maximum(S1r[:, 2, 2], 0.0))
        return 0.5 * d_tau + k * sigma_n_max, d_tau
