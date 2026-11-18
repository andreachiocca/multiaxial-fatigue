"""Boehme-Papuga (BP) integral stress criterion, Boehme & Papuga (2024), doi:10.1111/ffe.14281.

Equation (6) and the BP coefficients in Boehme, Papuga & Lange (2026), doi:10.1111/ffe.70244::

    sigma_eq^2 = <a tau_a^2 + b sigma_a^2 + c tau_a tau_m + d sigma_a sigma_m>   (sphere average)
    a = 9/2 (kappa^2 - 4/3),  b = 3 (3 - kappa^2),  kappa = s_-1 / t_-1
    c = 5/2 kappa^2 ((2 t_-1 / t_0)^2 - 1),  d = 5 ((2 s_-1 / s_0)^2 - 1 - 2 c / 15)

with s_0, t_0 the R = 0 strengths as maximum stresses, estimated by Eqs. (1)-(2) when missing.
The returned damage parameter is sigma_eq / s_-1. The input is a single-harmonic stress history.

Repository choices:
  - tau_a is the minimum-circumscribed-circle radius of the shear ellipse on each plane;
  - 24 x 48 Gauss-Legendre/azimuth quadrature instead of the MATLAB angular grid;
  - a negative integral is clipped to zero; the admissible kappa range is not enforced;
  - the index is used as the damage parameter of the shared master curve DP = A N^b.

User guide: docs/method-guides/BP.md.
"""

from functools import lru_cache
from typing import Mapping

import numpy as np

from multifat.models.errors import MissingParamsError

from multifat.models.base import CaseResult


@lru_cache(maxsize=8)
def _sphere_grid(n_polar: int = 24, n_azimuth: int = 48):
    """Normals and weights for the sphere average, including sin(phi) dphi."""
    z, wz = np.polynomial.legendre.leggauss(n_polar)
    theta = 2.0 * np.pi * np.arange(n_azimuth) / n_azimuth
    radius = np.sqrt(1.0 - z[:, None] ** 2)
    n = np.stack(
        np.broadcast_arrays(radius * np.cos(theta), radius * np.sin(theta), z[:, None]), axis=-1
    )
    normals = n.reshape(-1, 3)
    weights = np.repeat(wz / (2.0 * n_azimuth), n_azimuth)
    return normals, weights


def _strengths(params: Mapping[str, float], name="BP"):
    f_m1 = float(params["Sigm1"])
    t_m1 = float(params["Taum1"])
    if not (np.isfinite(f_m1) and np.isfinite(t_m1) and f_m1 > 0 and t_m1 > 0):
        raise ValueError(f"{name} requires positive finite Sigm1 and Taum1")

    f_0 = params.get("Sig0")
    if f_0 is None:
        su = params.get("Su")
        if su is None or not np.isfinite(su) or su <= 0:
            raise MissingParamsError(
                name, ["Sig0 or Su"], "Su estimates Sig0 by Eq. (1) of Boehme et al. (2026)"
            )
        f_0 = 4.0 * f_m1 * float(su) / (f_m1 + 2.0 * float(su))
    f_0 = float(f_0)
    if not np.isfinite(f_0) or f_0 <= 0:
        raise ValueError(f"{name} requires a positive finite Sig0")

    t_0 = params.get("Tau0")
    if t_0 is None:
        # Eq. (2): 4*t_-1/t_0 - 2*s_-1/s_0 = 1.
        t_0 = 4.0 * t_m1 / (1.0 + 2.0 * f_m1 / f_0)
    t_0 = float(t_0)
    if not np.isfinite(t_0) or t_0 <= 0:
        raise ValueError(f"{name} requires a positive finite Tau0")
    return f_m1, t_m1, f_0, t_0


def _plane_components(stress: np.ndarray, normals: np.ndarray):
    traction = normals @ stress.T
    normal = np.sum(traction * normals, axis=1)
    shear = traction - normal[:, None] * normals
    return normal, shear


def _plane_amplitudes(S_mean, S_sin, S_cos, normals):
    sigma_m, tau_m_vec = _plane_components(S_mean, normals)
    sigma_s, tau_s = _plane_components(S_sin, normals)
    sigma_c, tau_c = _plane_components(S_cos, normals)

    sigma_a = np.hypot(sigma_s, sigma_c)
    # Minimum circumscribed circle radius of the single-harmonic shear ellipse.
    ss = np.sum(tau_s * tau_s, axis=1)
    cc = np.sum(tau_c * tau_c, axis=1)
    sc = np.sum(tau_s * tau_c, axis=1)
    tau_a = np.sqrt(np.maximum(0.0, 0.5 * (ss + cc + np.hypot(ss - cc, 2.0 * sc))))
    tau_m = np.linalg.norm(tau_m_vec, axis=1)
    return sigma_a, sigma_m, tau_a, tau_m


def refine_plane_max(objective, normals, values, n_seeds: int = 12):
    """Largest objective over the grid, refined by pattern search from the n_seeds best normals
    and the three axes. objective maps (k, 3) unit normals to k values."""
    seeds = normals[np.argsort(values)[-n_seeds:]]
    seeds = np.concatenate((seeds, np.eye(3)))
    offsets = np.concatenate((np.eye(3), -np.eye(3), np.zeros((1, 3))))
    step = 0.15
    for _ in range(24):
        candidates = seeds[:, None, :] + step * offsets[None, :, :]
        candidates /= np.linalg.norm(candidates, axis=2, keepdims=True)
        refined = objective(candidates.reshape(-1, 3)).reshape(len(seeds), -1)
        seeds = candidates[np.arange(len(seeds)), np.argmax(refined, axis=1)]
        step *= 0.65
    return max(float(np.max(values)), float(np.max(refined)))


class BohmePapuga:
    """Orientation-averaged BP damage parameter, normalized by Sigm1."""

    name = "BP"
    optional_params = {
        "Sig0": "Boehme et al. 2026 Eq. (1) from Su and Sigm1",
        "Tau0": "Boehme et al. 2026 Eq. (2)",
        "Su": "only to estimate Sig0",
    }

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def effective_strengths(self, params: Mapping[str, float]) -> tuple[float, float, float, float]:
        """(Sigm1, Taum1, Sig0, Tau0), estimating Sig0 and Tau0 from Su if absent."""
        return _strengths(params, self.name)

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
        f_m1, t_m1, f_0, t_0 = _strengths(params, self.name)
        kappa_sq = (f_m1 / t_m1) ** 2
        a_BP, b_BP = self._amplitude_coefficients(kappa_sq)
        c_BP = 2.5 * kappa_sq * ((2.0 * t_m1 / t_0) ** 2 - 1.0)
        d_BP = 5.0 * ((2.0 * f_m1 / f_0) ** 2 - 1.0 - 2.0 * c_BP / 15.0)

        if harmonics is None:
            S_mean = 0.5 * (np.asarray(S0, dtype=float) + np.asarray(S1, dtype=float))
            S_sin = 0.5 * (np.asarray(S0, dtype=float) - np.asarray(S1, dtype=float))
            S_cos = np.zeros((3, 3))
        else:
            S_mean, S_sin, S_cos = (
                np.asarray(harmonics[key], dtype=float) for key in ("S_mean", "S_sin", "S_cos")
            )
        if any(x.shape != (3, 3) or not np.all(np.isfinite(x)) for x in (S_mean, S_sin, S_cos)):
            raise ValueError(f"{self.name} requires finite 3x3 stress mean/sin/cos tensors")

        normals, weights = _sphere_grid()
        sigma_a, sigma_m, tau_a, tau_m = _plane_amplitudes(S_mean, S_sin, S_cos, normals)
        integral = self._combine(
            a_BP,
            b_BP,
            c_BP,
            d_BP,
            sigma_a,
            sigma_m,
            tau_a,
            tau_m,
            weights,
            normals,
            S_mean,
            S_sin,
            S_cos,
        )
        return CaseResult(values={"": float(np.sqrt(max(0.0, integral)) / f_m1)})

    def _amplitude_coefficients(self, kappa_sq):
        return 4.5 * (kappa_sq - 4.0 / 3.0), 3.0 * (3.0 - kappa_sq)

    def _combine(
        self,
        a_BP,
        b_BP,
        c_BP,
        d_BP,
        sigma_a,
        sigma_m,
        tau_a,
        tau_m,
        weights,
        normals,
        S_mean,
        S_sin,
        S_cos,
    ):
        return np.dot(
            weights,
            a_BP * tau_a**2 + b_BP * sigma_a**2 + c_BP * tau_a * tau_m + d_BP * sigma_a * sigma_m,
        )
