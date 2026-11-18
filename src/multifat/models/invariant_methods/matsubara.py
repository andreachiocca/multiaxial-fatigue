"""Matsubara, Hayashida & Kano (2018) unified equivalent shear stress amplitude (N09).

G. Matsubara, A. Hayashida, D. Kano, Int. J. Fatigue (2018), accepted manuscript,
doi:10.1016/j.ijfatigue.2017.12.001, Eqs. (7)-(9), (13)-(15) and Section 2.2.3.4.1::

    sqrt(J_2,amp) = sqrt(sum_i ((S_i,max - S_i,min)/2)^2),  i = 1..5                       (13)
    S_1 = sqrt(3)/2 s_xx, S_2 = (s_yy - s_zz)/2, S_3 = s_xy, S_4 = s_yz, S_5 = s_zx           (8)
    S_max = maximum over the cycle of the normal stress with the largest range              (2.2.3.2)
    sqrt(J_2,amp,eq) = sigma_T sqrt(J_2,amp) / (sigma_T - S_max + sqrt(3) sqrt(J_2,amp))   (15)
    sqrt(J_2,amp) + (sigma_w/sqrt3) / (sigma_T - sigma_w) S_max
        = (sigma_w/sqrt3) sigma_T / (sigma_T - sigma_w)                                     (9)

sqrt(J_2,amp,eq) maps a stress state onto the fully reversed uniaxial line
sqrt(J_2,amp)/S_max = 1/sqrt(3) at the same fatigue life (Fig. 4); there it equals sigma_a/sqrt(3).
The life follows from the unified equivalent S-N curve, the reversed axial S-N curve with the
ordinate sigma_a/sqrt(3) (Section 2.2.3.4.1). sqrt(J_2,amp,eq) = sigma_w/sqrt(3) is Eq. (9).

For one harmonic the half-ranges of S_i are the amplitudes of the components, so
J_2,amp = (|dev S_sin|^2 + |dev S_cos|^2)/2, the Mamiya and Li values coinciding (Section 2.2.3.1).
The normal stress on the plane of unit normal n is n.S(t).n, with amplitude
sqrt((n.S_sin.n)^2 + (n.S_cos.n)^2) and mean n.S_mean.n; S_max is mean plus amplitude on the plane
of the largest amplitude, found by a search over the unit sphere.

Repository choices:
  - sigma_T is the optional constant `sigma_T`; without it sigma_T = 1.84 Su, the estimate the paper
    states for 400 MPa <= Su <= 1000 MPa (its ref. [23]). Outside that range, without sigma_T, the
    material is out of scope (MaterialScopeError) (maintainer decision);
  - finite lives come from the shared master curve DP = A N^b with DP = sqrt(J_2,amp,eq), fitted on
    the R = -1 uniaxial tension failures: this is the paper's unified equivalent S-N curve from the
    reversed axial S-N curve (Section 2.2.3.4.1). The alternative estimate from Su and sigma_w
    (Eqs. 16-17) is not implemented;
  - S_max: when several planes share the largest normal stress range (within 1e-9), the one with
    the largest S_max is taken (the conservative choice; the paper does not discuss ties);
  - a state with S_max >= sigma_T + sqrt(3) sqrt(J_2,amp) leaves Eq. (15) undefined and raises
    LoadingScopeError.

User guide: docs/method-guides/MATSUBARA.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError, MaterialScopeError, MissingParamsError
from .stress_invariants import harmonic_parts

SQRT3 = np.sqrt(3.0)
SIGMA_T_FACTOR = 1.84  # sigma_T / Su, Section 2.2.3.3 (ref. [23]), for 400 <= Su <= 1000 MPa
SU_RANGE = (400.0, 1000.0)
TIE_RTOL = 1e-9


def true_fracture_strength(params: Mapping[str, float]) -> float:
    """sigma_T from the constant `sigma_T`, or 1.84 Su within 400-1000 MPa."""
    sigma_T = params.get("sigma_T")
    if sigma_T is not None:
        return float(sigma_T)
    su = params.get("Su")
    if su is None:
        raise MissingParamsError("MATSUBARA", ["sigma_T"], "or Su to estimate it (1.84 Su)")
    su = float(su)
    if not SU_RANGE[0] <= su <= SU_RANGE[1]:
        raise MaterialScopeError(
            f"MATSUBARA: sigma_T = 1.84 Su holds for 400 <= Su <= 1000 MPa (Su = {su:g} MPa); "
            "give sigma_T in the Summary"
        )
    return SIGMA_T_FACTOR * su


def equivalent_shear_amplitude(S_sin, S_cos) -> float:
    """sqrt(J_2,amp) of Eq. (13) for one harmonic."""
    total = 0.0
    for T in (np.asarray(S_sin, float), np.asarray(S_cos, float)):
        dev = T - np.trace(T) / 3.0 * np.eye(3)
        total += 0.5 * float(np.sum(dev**2))
    return float(np.sqrt(total))


def _sphere(n: int) -> np.ndarray:
    """n quasi-uniform unit normals on the upper hemisphere (Fibonacci lattice)."""
    k = np.arange(n) + 0.5
    z = k / n
    phi = np.pi * (1.0 + np.sqrt(5.0)) * k
    r = np.sqrt(1.0 - z * z)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def _normal_amplitude(normals, S_sin, S_cos) -> np.ndarray:
    a = np.einsum("ni,ij,nj->n", normals, S_sin, normals)
    b = np.einsum("ni,ij,nj->n", normals, S_cos, normals)
    return np.hypot(a, b)


def _tangents(n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    helper = np.where(np.abs(n[:, :1]) < 0.9, [[1.0, 0.0, 0.0]], [[0.0, 1.0, 0.0]])
    t1 = np.cross(n, helper)
    t1 /= np.linalg.norm(t1, axis=1)[:, None]
    return t1, np.cross(n, t1)


_GRID = _sphere(2000)
_OFFSETS = np.linspace(-1.0, 1.0, 7)


def max_range_normal_stress(S_mean, S_sin, S_cos) -> float:
    """S_max: maximum over the cycle of the normal stress whose range is the largest."""
    S_mean, S_sin, S_cos = (np.asarray(T, float) for T in (S_mean, S_sin, S_cos))
    amplitude = _normal_amplitude(_GRID, S_sin, S_cos)
    seeds = [_GRID[np.argsort(amplitude)[-12:]]]
    for T in (S_sin, S_cos, S_sin + S_cos, S_sin - S_cos):  # exact maxima of proportional paths
        seeds.append(np.linalg.eigh(T)[1].T)
    best = np.concatenate(seeds)
    best = best * np.where(best[:, 2:] < 0.0, -1.0, 1.0)  # n and -n are the same plane
    step = 0.05
    for _ in range(45):  # shrinking local searches around every seed at once
        t1, t2 = _tangents(best)
        g = step * _OFFSETS
        trial = (
            best[:, None, None, :]
            + g[None, :, None, None] * t1[:, None, None, :]
            + g[None, None, :, None] * t2[:, None, None, :]
        ).reshape(len(best), -1, 3)
        trial /= np.linalg.norm(trial, axis=2)[:, :, None]
        amp = _normal_amplitude(trial.reshape(-1, 3), S_sin, S_cos).reshape(len(best), -1)
        best = trial[np.arange(len(best)), np.argmax(amp, axis=1)]
        step *= 0.7
    amplitude = _normal_amplitude(best, S_sin, S_cos)
    mean = np.einsum("ni,ij,nj->n", best, S_mean, best)
    top = amplitude >= amplitude.max() * (1.0 - TIE_RTOL) - 1e-12
    return float(np.max(mean[top] + amplitude[top]))


def unified_equivalent_shear(S_mean, S_sin, S_cos, sigma_T: float) -> float:
    """sqrt(J_2,amp,eq) of Eq. (15)."""
    sqrt_j2 = equivalent_shear_amplitude(S_sin, S_cos)
    s_max = max_range_normal_stress(S_mean, S_sin, S_cos)
    denominator = float(sigma_T) - s_max + SQRT3 * sqrt_j2
    if not denominator > 0.0:
        raise LoadingScopeError(
            f"MATSUBARA: Eq. (15) undefined, S_max = {s_max:.1f} MPa reaches sigma_T + "
            f"sqrt(3) sqrt(J_2,amp) = {float(sigma_T) + SQRT3 * sqrt_j2:.1f} MPa"
        )
    return float(sigma_T) * sqrt_j2 / denominator


class Matsubara:
    name = "MATSUBARA"
    optional_params = {
        "sigma_T": "true fracture strength; default 1.84 Su for 400 <= Su <= 1000 MPa (ref. [23])",
        "Su": "used only to estimate sigma_T",
    }

    def required_params(self) -> set[str]:
        return set()

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
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        sigma_T = true_fracture_strength(params)
        return CaseResult(values={"": unified_equivalent_shear(S_mean, S_sin, S_cos, sigma_T)})
