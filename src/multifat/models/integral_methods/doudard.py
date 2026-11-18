"""Probabilistic two-scale fatigue criterion of Doudard, Hild & Calloch (2007).

Published in *Fatigue & Fracture of Engineering Materials & Structures* 30 (2007),
pp. 107-114, doi:10.1111/j.1460-2695.2006.01087.x (roadmap item N34), Eqs. (15)-(36)::

    P_F = 1 - exp(-(V * kappa / V0) * (Sigma_T / S0)^m)
    kappa(m) = int (2 * T_a / Sigma_T)^m * sin(zeta) dchi dzeta dxi   (measure 8 pi^2)
    N = A' / (<Sigma_T - Sigma_inf(P_F)> * Sigma_inf(P_F))

with Tresca amplitude ``Sigma_T`` (Eq. 26), effective volume ``V_eff = V * kappa``,
and statistical fatigue-limit moments (Eqs. 19-20).

Repository choices:
- Provided as uniform-stress point kernels rather than a registered model because the
  workbooks lack both the stressed volume ratio ``V / V0`` and Weibull scale ``S0``.
- Quadrature for ``kappa(m)`` uses ``KAPPA_GRID = (48, 64, 128)`` Gauss-Legendre nodes
  in ``cos(zeta)`` and uniform ``(xi, chi)`` with Eq. (36) ``sin(zeta)`` (fixing the
  ``cos(zeta)`` misprint in Eq. 14 and ``2/pi`` -> ``1/(2 pi)`` in Eqs. A4-A5). The
  printed ``kappa(16)`` values exceed direct quadrature by 5-10% (unresolved).
"""

from __future__ import annotations

from math import gamma

import numpy as np

KAPPA_GRID = (48, 64, 128)  # Gauss-Legendre nodes in cos(zeta), uniform xi, uniform chi


def tension_torsion_amplitudes(axial_stress: float, shear_stress: float, phase_shift: float):
    """Sine and cosine amplitude tensors of Eqs. (23)-(25), i.e. S(t) = A sin(wt) + B cos(wt)."""
    S_sin = np.zeros((3, 3))
    S_cos = np.zeros((3, 3))
    S_sin[0, 0] = axial_stress
    S_sin[0, 1] = S_sin[1, 0] = shear_stress * np.cos(phase_shift)
    S_cos[0, 1] = S_cos[1, 0] = shear_stress * np.sin(phase_shift)
    return S_sin, S_cos


def tresca_amplitude(axial_stress: float, shear_stress: float, phase_shift: float) -> float:
    """Sigma_T according to Eq. (26)."""
    invariant = axial_stress**2 + 4.0 * shear_stress**2
    disc = np.sqrt(
        max(
            invariant**2 - 16.0 * axial_stress**2 * shear_stress**2 * np.sin(phase_shift) ** 2,
            0.0,
        )
    )
    return float(np.sqrt((invariant + disc) / 2.0))


def kappa(
    m: float,
    S_sin,
    S_cos,
    sigma_T: float,
    grid: tuple[int, int, int] = KAPPA_GRID,
) -> float:
    """kappa(m) from Eq. (36) for the single-harmonic stress A sin + B cos."""
    n_cos_zeta, n_azimuth, n_in_plane = grid
    cos_zeta, gauss_weights = np.polynomial.legendre.leggauss(n_cos_zeta)
    azimuth = 2.0 * np.pi * np.arange(n_azimuth) / n_azimuth
    chi = 2.0 * np.pi * np.arange(n_in_plane) / n_in_plane
    cz_mesh, az_mesh = np.meshgrid(cos_zeta, azimuth, indexing="ij")
    sin_zeta = np.sqrt(1.0 - cz_mesh**2)
    normals = np.stack([sin_zeta * np.cos(az_mesh), sin_zeta * np.sin(az_mesh), cz_mesh], axis=-1)
    radial = np.stack([-cz_mesh * np.cos(az_mesh), -cz_mesh * np.sin(az_mesh), sin_zeta], axis=-1)
    tangent = np.stack([-np.sin(az_mesh), np.cos(az_mesh), np.zeros_like(az_mesh)], axis=-1)
    S_sin_proj = normals @ np.asarray(S_sin, float).T
    S_cos_proj = normals @ np.asarray(S_cos, float).T
    S_sin_tan, S_sin_rad = np.sum(tangent * S_sin_proj, -1), np.sum(radial * S_sin_proj, -1)
    S_cos_tan, S_cos_rad = np.sum(tangent * S_cos_proj, -1), np.sum(radial * S_cos_proj, -1)
    cos_chi, sin_chi = np.cos(chi), np.sin(chi)
    shear_sq = (S_sin_tan[..., None] * cos_chi + S_sin_rad[..., None] * sin_chi) ** 2 + (
        S_cos_tan[..., None] * cos_chi + S_cos_rad[..., None] * sin_chi
    ) ** 2
    integrand = (4.0 * shear_sq / sigma_T**2) ** (m / 2.0)
    quad_weights = (
        gauss_weights[:, None, None] * (2.0 * np.pi / n_azimuth) * (2.0 * np.pi / n_in_plane)
    )
    return float(np.sum(quad_weights * integrand))


def failure_probability(sigma_T: float, m: float, s0: float, volume_ratio: float) -> float:
    """P_F from Eq. (15), with V_eff / V0 = V H~_m / V0 (= V kappa / V0 under uniform stress)."""
    return float(1.0 - np.exp(-volume_ratio * (sigma_T / s0) ** m))


def fatigue_limit(probability: float, m: float, s0: float, volume_ratio: float) -> float:
    """Sigma_inf(P_F): Eq. (15) rearranged to give Sigma_T."""
    return float(s0 * (-np.log1p(-probability) / volume_ratio) ** (1.0 / m))


def mean_fatigue_limit(m: float, s0: float, volume_ratio: float) -> float:
    """First moment of Eq. (19)."""
    return float(s0 * volume_ratio ** (-1.0 / m) * gamma(1.0 + 1.0 / m))


def fatigue_limit_sd(m: float, s0: float, volume_ratio: float) -> float:
    """Second central moment of Eq. (20)."""
    variance_term = gamma(1.0 + 2.0 / m) - gamma(1.0 + 1.0 / m) ** 2
    return float(s0 * volume_ratio ** (-1.0 / m) * np.sqrt(variance_term))


def stromeyer_life(sigma_T: float, sigma_inf: float, A_prime: float) -> float:
    """Eq. (22): returns infinity when the amplitude does not exceed the fatigue limit."""
    margin = sigma_T - sigma_inf
    return float(A_prime / (margin * sigma_inf)) if margin > 0 else float("inf")
