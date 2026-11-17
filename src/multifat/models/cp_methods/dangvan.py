"""Dang Van criterion in the formulation restated by Boehme, Papuga & Lange (2026), Eq. (3), p. 4.

    sigma_eq,DV = max over planes (a tau_sa) + b sigma_H,max,   a = kappa,  b = 3 - 3 kappa / 2,
    kappa = s_-1 / t_-1 (Eq. 4)

following the modified form of Papuga & Halama. This is a documented variant: Dang Van (1973)
defines the criterion via the mesoscopic shakedown stress and does not give closed-form expressions
for a and b. The calibration reproduces s_-1 in fully reversed tension and t_-1 in torsion, together
with the closed forms of Eq. (16) for equibiaxial and equitriaxial loading. tau_sa is the radius of
the smallest circle circumscribing the shear-stress path on the plane (exact for a single harmonic,
as in BP and CAIM).

This is a case-level model (it uses the whole single-harmonic history) and returns sigma_eq,DV / s_-1,
like the other criteria of Boehme et al. (2026).

Repository choices:
  - the amplitude-plus-maximum variant of Eq. (3), not the maximum over time of the mesoscopic
    sum; identical for proportional loading, conservative for single-harmonic out-of-phase
    loading when a = 3 (t_-1/s_-1 - 1/2) > 0;
  - the index is used as the damage parameter of the shared master curve DP = A N^b.

User guide: docs/method-guides/DANGVAN.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..integral_methods.bohme_papuga import (
    _plane_amplitudes,
    _sphere_grid,
    refine_plane_max,
)


def dang_van_equivalent(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    """Apply Eq. (3) to S(t) = m + p sin(wt) + q cos(wt)."""
    kappa = f_m1 / t_m1
    normals, _ = _sphere_grid()

    def tau_sa(normal):
        return _plane_amplitudes(S_mean, S_sin, S_cos, normal)[2]

    tau_sa_max = refine_plane_max(tau_sa, normals, tau_sa(normals))
    sigma_H_max = (np.trace(S_mean) + np.hypot(np.trace(S_sin), np.trace(S_cos))) / 3.0
    return float(kappa * tau_sa_max + (3.0 - 1.5 * kappa) * sigma_H_max)


class DangVan:
    name = "DANGVAN"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self,
        *,
        S0,
        E0,
        S1,
        E1,
        params: Mapping[str, float],
        R_list=None,
        harmonics=None,
    ) -> CaseResult:
        f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
        if not (np.isfinite(f_m1) and np.isfinite(t_m1) and f_m1 > 0 and t_m1 > 0):
            raise ValueError("DANGVAN needs positive, finite Sigm1 and Taum1")
        if harmonics is None:
            S0, S1 = np.asarray(S0, dtype=float), np.asarray(S1, dtype=float)
            S_mean, S_sin, S_cos = (
                0.5 * (S0 + S1),
                0.5 * (S0 - S1),
                np.zeros((3, 3)),
            )
        else:
            S_mean, S_sin, S_cos = (
                np.asarray(harmonics[key], dtype=float) for key in ("S_mean", "S_sin", "S_cos")
            )
        return CaseResult(values={"": dang_van_equivalent(S_mean, S_sin, S_cos, f_m1, t_m1) / f_m1})
