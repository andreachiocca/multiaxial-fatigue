"""Papuga & Růžička (2008), Int J Fatigue 30:58-66, doi:10.1016/j.ijfatigue.2007.02.015:
the integral criterion PI (roadmap G02) and the critical-plane criterion PCr (G01).

PI, Eq. (16):   sqrt( 1/(4 pi) int int [a_I C_a^2 + b_I (N_a + (f_-1/t_-1) N_m)] sin psi dpsi dphi )
                <= f_-1,   a_I = 5 kappa^2 / 2,   b_I = f_-1 (3 - kappa^2)
PCr, Eq. (17):  max over planes sqrt( a_c C_a^2 + b_c (N_a + (t_-1/f_0) N_m) ) <= f_-1,
                kappa <  sqrt(4/3): a_c = (kappa^2 + sqrt(kappa^4 - kappa^2))/2,  b_c = f_-1
                kappa >= sqrt(4/3): a_c = (4 kappa^2/(4 + kappa^2))^2,
                                    b_c = 8 f_-1 kappa^2 (4 - kappa^2)/(4 + kappa^2)^2

Here kappa = f_-1/t_-1, and f_0 denotes the axial fatigue limit under repeated loading (R = 0),
taken as a maximum stress (in the paper it is obtained from Morrow's formula, Eq. 1). C_a is the
radius of the smallest circle enclosing the shear path on a given plane (p.4, following
Papadopoulos et al. 1997), while N_a and N_m are the amplitude and the mean of the normal stress.
Both criteria return LHS / f_-1 (equal to 1 on the fatigue limit; the fatigue index error of
Eq. (5) follows as 100 (index - 1)).

Repository choices: the PI sphere integral is evaluated with BP quadrature over a 48 x 192 grid
(weights summing to 1 = 1/(4 pi) int int sin psi); PCr is maximised on the default 24 x 48 grid of
``_sphere_grid``, refined by pattern search; a negative radicand (compressive mean normal stress) is clipped to 0; when Sig0 is
not provided, f_0 is estimated from Su through the Boehme et al. (2026) Eq. (1) relation, as done
for BP/CAIM.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..invariant_methods.stress_invariants import harmonic_parts
from .bohme_papuga import _plane_amplitudes, _sphere_grid, _strengths, refine_plane_max

# N_a = |n.A.n| is not smooth; 48 x 192 integrates it to roughly 4e-4 (24 x 48: 6e-3).
PI_GRID = (48, 192)


def pi_parameters(f_m1: float, t_m1: float) -> tuple[float, float]:
    kappa_sq = (f_m1 / t_m1) ** 2
    return 2.5 * kappa_sq, f_m1 * (3.0 - kappa_sq)


def pcr_parameters(f_m1: float, t_m1: float) -> tuple[float, float]:
    kappa_sq = (f_m1 / t_m1) ** 2
    if kappa_sq < 4.0 / 3.0:
        return 0.5 * (kappa_sq + np.sqrt(max(kappa_sq**2 - kappa_sq, 0.0))), f_m1
    return (4.0 * kappa_sq / (4.0 + kappa_sq)) ** 2, 8.0 * f_m1 * kappa_sq * (4.0 - kappa_sq) / (
        4.0 + kappa_sq
    ) ** 2


def pi_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    a_I, b_I = pi_parameters(f_m1, t_m1)
    normals, quad_weights = _sphere_grid(*PI_GRID)
    sigma_a, sigma_m, tau_a, _ = _plane_amplitudes(S_mean, S_sin, S_cos, normals)
    sphere_int = np.sum(quad_weights * (a_I * tau_a**2 + b_I * (sigma_a + f_m1 / t_m1 * sigma_m)))
    return float(np.sqrt(max(sphere_int, 0.0)) / f_m1)


def pcr_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float, f_0: float) -> float:
    a_c, b_c = pcr_parameters(f_m1, t_m1)

    def under_root(normals):
        sigma_a, sigma_m, tau_a, _ = _plane_amplitudes(S_mean, S_sin, S_cos, normals)
        return a_c * tau_a**2 + b_c * (sigma_a + t_m1 / f_0 * sigma_m)

    normals, _ = _sphere_grid()
    peak = refine_plane_max(under_root, normals, under_root(normals))
    return float(np.sqrt(max(peak, 0.0)) / f_m1)


class PapugaIntegral:
    """PI criterion, Eq. (16), given as LHS / f_-1."""

    name = "PI"

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
            raise ValueError("PI requires positive finite Sigm1 and Taum1")
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": pi_index(S_mean, S_sin, S_cos, f_m1, t_m1)})


class PapugaCriticalPlane:
    """PCr criterion, Eq. (17), given as LHS / f_-1."""

    name = "PCR"

    optional_params = {
        "Sig0": "Boehme et al. 2026 Eq. (1) from Su and Sigm1",
        "Su": "only to estimate Sig0",
    }

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
        f_m1, t_m1, f_0, _ = _strengths(params, name="PCR")
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": pcr_index(S_mean, S_sin, S_cos, f_m1, t_m1, f_0)})
