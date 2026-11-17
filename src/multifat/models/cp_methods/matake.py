"""Matake (1977), Bull. JSME, doi:10.1299/jsme1958.20.257: critical-plane fatigue-limit criterion.

Eq. (1), p. 5, on the plane of maximum shear stress amplitude, with the coefficient of the
Appendix, p. 7 (line through (0, t_-1) and (f_-1/2, f_-1/2))::

    tau_a + mu sigma_n,max = t_-1,        mu = 2 t_-1 / f_-1 - 1

The plane damage parameter is DP = d_tau / 2 + mu sigma_max (MPa), equal to t_-1 at both
fully reversed limits; the selection metric is d_tau and MATAKE_ext takes the largest DP.

Repository choices:
  - sigma_max is clipped at zero (tension-only normal stress);
  - out-of-phase shear amplitude: half the longest chord of the shear path (``longest_chord``),
    i.e. the minimum circumscribed circle radius for one harmonic, as Papadopoulos et al. (1997)
    restate Matake; their Table 1 is reproduced only with it (maintainer decision, 17 Nov 2026;
    the rss rule was used before);
  - life from the shared master curve DP = A N^b.

User guide: docs/method-guides/MATAKE.md.
"""

from __future__ import annotations

from typing import Mapping, Tuple
import numpy as np
from multifat.models.utils.tensors import as33 as _as33
from multifat.models.utils.batch import check_planes, pypow2

from ..base import CPPlaneResult


def _matake_plane_value(
    *,
    S0r: np.ndarray,
    S1r: np.ndarray,
    f_m1: float,
    t_m1: float,
) -> Tuple[float, float]:
    """
    Matake criterion (stress form) on one candidate plane.

    Critical plane: maximum shear stress range (MSSR).
    Damage-like scalar (for correlation): tau_a + k * sigma_n,max
      where k = 2*Taum1/Sigm1 - 1  (common approximation using fatigue limits)

    Conventions:
      - candidate plane normal is local z-axis
      - shear stress components on plane: Sxz and Syz
      - normal stress: Szz (clipped at 0 when taking maximum)
      - two load steps: (0=max), (1=min)
    """
    S0r = _as33("S0r", S0r)
    S1r = _as33("S1r", S1r)

    f_m1 = float(f_m1)
    t_m1 = float(t_m1)

    if f_m1 <= 0.0:
        raise ValueError(f"Sigm1 (axial fatigue limit) must be > 0. Got {f_m1}")
    if t_m1 <= 0.0:
        raise ValueError(f"Taum1 (torsional fatigue limit) must be > 0. Got {t_m1}")

    # k from existing material properties (no new Excel fields)
    k = 2.0 * t_m1 / f_m1 - 1.0

    # Shear stress RANGE vector components on plane (xz, yz)
    d_tau_xz = float(S0r[0, 2] - S1r[0, 2])
    d_tau_yz = float(S0r[1, 2] - S1r[1, 2])
    d_tau = float(np.sqrt(d_tau_xz**2 + d_tau_yz**2))

    # shear stress amplitude
    tau_a = 0.5 * d_tau

    # max normal stress on plane (clip compression to 0, as in FS and FIN)
    sigma_z0 = max(float(S0r[2, 2]), 0.0)
    sigma_z1 = max(float(S1r[2, 2]), 0.0)
    sigma_n_max = max(sigma_z0, sigma_z1)

    damage = float(tau_a + k * sigma_n_max)

    # MSSR metric: maximize shear stress range (equivalently amplitude)
    metric = d_tau
    return damage, metric


class Matake:
    """
    Matake (MSSR-type) critical plane method.
    Uses fatigue limits to compute k internally:
      k = 2*Taum1/Sigm1 - 1
    """

    name = "MATAKE"
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_on_plane(
        self,
        *,
        S0r,
        E0r,  # unused (stress-based)
        S1r,
        E1r,  # unused (stress-based)
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        if params.get("Sigm1", None) is None:
            raise KeyError("MATAKE requires params['Sigm1'] (Tensile_fatigue_limit_MPa)")
        if params.get("Taum1", None) is None:
            raise KeyError("MATAKE requires params['Taum1'] (Torsional_fatigue_limit_MPa)")

        damage, metric = _matake_plane_value(
            S0r=S0r,
            S1r=S1r,
            f_m1=float(params["Sigm1"]),
            t_m1=float(params["Taum1"]),
        )
        return CPPlaneResult(damage=damage, metric=metric)

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        """All planes at once; bit-identical to ``evaluate_on_plane`` per plane."""
        check_planes(S0r, E0r, S1r, E1r)
        if params.get("Sigm1", None) is None:
            raise KeyError("MATAKE requires params['Sigm1'] (Tensile_fatigue_limit_MPa)")
        if params.get("Taum1", None) is None:
            raise KeyError("MATAKE requires params['Taum1'] (Torsional_fatigue_limit_MPa)")
        f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
        if f_m1 <= 0.0:
            raise ValueError(f"Sigm1 (axial fatigue limit) must be > 0. Got {f_m1}")
        if t_m1 <= 0.0:
            raise ValueError(f"Taum1 (torsional fatigue limit) must be > 0. Got {t_m1}")
        k = 2.0 * t_m1 / f_m1 - 1.0
        d_tau_xz = S0r[:, 0, 2] - S1r[:, 0, 2]
        d_tau_yz = S0r[:, 1, 2] - S1r[:, 1, 2]
        d_tau = np.sqrt(pypow2(d_tau_xz) + pypow2(d_tau_yz))
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], 0.0), np.maximum(S1r[:, 2, 2], 0.0))
        return 0.5 * d_tau + k * sigma_n_max, d_tau
