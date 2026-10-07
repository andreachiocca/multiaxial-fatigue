"""Böhme–Papuga (BP) integral stress criterion.

Equation (6) and the BP coefficients in Böhme, Papuga & Lange (2026),
doi:10.1111/ffe.70244. The input is a single-harmonic stress history.
"""

from functools import lru_cache
from typing import Mapping

import numpy as np

from fatigue.models.base import CaseResult


@lru_cache(maxsize=8)
def _sphere_grid(n_polar: int = 24, n_azimuth: int = 48):
    """Normals and weights for the sphere average, including sin(phi) dphi."""
    z, wz = np.polynomial.legendre.leggauss(n_polar)
    theta = 2.0 * np.pi * np.arange(n_azimuth) / n_azimuth
    radius = np.sqrt(1.0 - z[:, None] ** 2)
    n = np.stack(np.broadcast_arrays(radius * np.cos(theta),
                                      radius * np.sin(theta), z[:, None]), axis=-1)
    normals = n.reshape(-1, 3)
    weights = np.repeat(wz / (2.0 * n_azimuth), n_azimuth)
    return normals, weights


def _strengths(params: Mapping[str, float]):
    s = float(params["Sigm1"])
    t = float(params["Taum1"])
    if not (np.isfinite(s) and np.isfinite(t) and s > 0 and t > 0):
        raise ValueError("BP requires positive finite Sigm1 and Taum1")

    s0 = params.get("Sig0")
    if s0 is None:
        su = params.get("Su")
        if su is None or not np.isfinite(su) or su <= 0:
            raise ValueError("BP requires Sig0, or Su to estimate it using Eq. (1) of Böhme et al. (2026)")
        s0 = 4.0 * s * float(su) / (s + 2.0 * float(su))
    s0 = float(s0)
    if not np.isfinite(s0) or s0 <= 0:
        raise ValueError("BP requires a positive finite Sig0")

    t0 = params.get("Tau0")
    if t0 is None:
        # Eq. (2): 4*t_-1/t_0 - 2*s_-1/s_0 = 1.
        t0 = 4.0 * t / (1.0 + 2.0 * s / s0)
    t0 = float(t0)
    if not np.isfinite(t0) or t0 <= 0:
        raise ValueError("BP requires a positive finite Tau0")
    return s, t, s0, t0


def _plane_components(stress: np.ndarray, normals: np.ndarray):
    traction = normals @ stress.T
    normal = np.sum(traction * normals, axis=1)
    shear = traction - normal[:, None] * normals
    return normal, shear


class BohmePapuga:
    """Orientation-averaged BP damage parameter, normalized by Sigm1."""

    name = "BP"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float],
        R_list=None, harmonics=None,
    ) -> CaseResult:
        s, t, s0, t0 = _strengths(params)
        k2 = (s / t) ** 2
        a = 4.5 * (k2 - 4.0 / 3.0)
        b = 3.0 * (3.0 - k2)
        c = 2.5 * k2 * ((2.0 * t / t0) ** 2 - 1.0)
        d = 5.0 * ((2.0 * s / s0) ** 2 - 1.0 - 2.0 * c / 15.0)

        if harmonics is None:
            mean = 0.5 * (np.asarray(S0, dtype=float) + np.asarray(S1, dtype=float))
            sin = 0.5 * (np.asarray(S0, dtype=float) - np.asarray(S1, dtype=float))
            cos = np.zeros((3, 3))
        else:
            mean, sin, cos = (np.asarray(harmonics[key], dtype=float)
                              for key in ("S_mean", "S_sin", "S_cos"))
        if any(x.shape != (3, 3) or not np.all(np.isfinite(x)) for x in (mean, sin, cos)):
            raise ValueError("BP requires finite 3x3 stress mean/sin/cos tensors")

        normals, weights = _sphere_grid()
        sigma_m, tau_m_vec = _plane_components(mean, normals)
        sigma_s, tau_s = _plane_components(sin, normals)
        sigma_c, tau_c = _plane_components(cos, normals)

        sigma_a = np.hypot(sigma_s, sigma_c)
        # Minimum circumscribed circle radius of the single-harmonic shear ellipse.
        ss = np.sum(tau_s * tau_s, axis=1)
        cc = np.sum(tau_c * tau_c, axis=1)
        sc = np.sum(tau_s * tau_c, axis=1)
        tau_a = np.sqrt(np.maximum(0.0, 0.5 * (ss + cc + np.hypot(ss - cc, 2.0 * sc))))
        tau_m = np.linalg.norm(tau_m_vec, axis=1)

        integral = np.dot(weights, a * tau_a**2 + b * sigma_a**2
                          + c * tau_a * tau_m + d * sigma_a * sigma_m)
        return CaseResult(values={"": float(np.sqrt(max(0.0, integral)) / s)})
