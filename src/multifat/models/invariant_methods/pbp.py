"""Projection-by-Projection (PbP) invariant fatigue criterion.

Published by Cristofori, Susmel & Tovo (2008), *International Journal of Fatigue* 30,
pp. 1646-1658, doi:10.1016/j.ijfatigue.2007.11.006 (roadmap N37), Eqs. (4)-(13)::

    Gamma = sqrt(sum_i sqrt(J_{2a,i})^2)
    rho_ref = sqrt(3) * sigma_{H,ref} / Gamma
    sqrt(J_{2A,rho}) = t_-1 + rho_ref * (f_-1 / sqrt(3) - t_-1)

The deviatoric stress 5-vector ``s(t)`` (``|s|^2 = J_2``, Eq. A.3) is projected along the
principal directions of its covariance matrix, cycle counted per axis, and combined via
the quadratic sum ``Gamma``.

Repository choices:
- Evaluation is limited to constant-amplitude single-harmonic histories, for which
  ``Gamma^2 = J_2(A) + J_2(B)`` holds identically in closed form independently of the
  covariance eigenbasis and ``sigma_{H,ref} = sigma_{H,a} + sigma_{H,m}``.
- ``rho_ref`` retains the ``sqrt(3)`` factor of Eq. (6) (omitted in Eq. B.5) so reversed
  axial loading yields ``rho_ref = 1``.
- The method returns the endurance-limit index ``Gamma / sqrt(J_{2A,rho})`` (``0.0`` when
  ``Gamma == 0``; ``inf`` when ``sqrt(J_{2A,rho}) <= 0``) for master-curve life prediction;
  Eqs. (9)-(11) life estimation with ``k_rho`` is not implemented.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from .stress_invariants import _deviatoric, harmonic_parts, hydrostatic_amplitude_mean


def deviatoric_vector(sigma: np.ndarray) -> np.ndarray:
    """Return the 5-vector whose squared norm equals J2 for stress tensors sigma (..., 3, 3), Eq. (A.3)."""
    sigma = np.asarray(sigma, dtype=float)
    hydro_part = np.trace(sigma, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    s = sigma - hydro_part
    return np.stack(
        [
            np.sqrt(3.0) / 2.0 * s[..., 0, 0],
            0.5 * (s[..., 1, 1] - s[..., 2, 2]),
            s[..., 0, 1],
            s[..., 0, 2],
            s[..., 1, 2],
        ],
        axis=-1,
    )


def projection_amplitudes(history: np.ndarray) -> np.ndarray:
    """Generic procedure applied to a sampled periodic history (n, 3, 3): build the covariance
    eigenbasis of the deviatoric path (Eqs. 4-5, with 1/T) and take the half range of each
    projection (a single cycle per projection for a harmonic signal)."""
    s = deviatoric_vector(history)
    cov = np.cov(s.T, bias=True)
    _, eigvecs = np.linalg.eigh(cov)
    projections = s @ eigvecs
    return 0.5 * (projections.max(axis=0) - projections.min(axis=0))


def pbp_amplitude(S_sin, S_cos) -> float:
    """Gamma for a single harmonic component: sqrt(J2(A) + J2(B))."""
    S_sin_dev, S_cos_dev = _deviatoric(S_sin), _deviatoric(S_cos)
    return float(np.sqrt(0.5 * (np.sum(S_sin_dev * S_sin_dev) + np.sum(S_cos_dev * S_cos_dev))))


def pbp_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    Gamma = pbp_amplitude(S_sin, S_cos)
    sigma_H_a, sigma_H_m = hydrostatic_amplitude_mean(S_mean, S_sin, S_cos)
    if Gamma == 0.0:
        return 0.0
    rho_ref = np.sqrt(3.0) * (sigma_H_a + sigma_H_m) / Gamma
    sqrt_J2A_rho = t_m1 + rho_ref * (f_m1 / np.sqrt(3.0) - t_m1)
    return float(Gamma / sqrt_J2A_rho) if sqrt_J2A_rho > 0.0 else float("inf")


class ProjectionByProjection:
    """Fatigue-limit index of the PbP criterion, Gamma / sqrt(J2A,rho_ref)."""

    name = "PBP"

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
            raise ValueError("PBP needs Sigm1 and Taum1 to be positive and finite")
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": pbp_index(S_mean, S_sin, S_cos, f_m1, t_m1)})
