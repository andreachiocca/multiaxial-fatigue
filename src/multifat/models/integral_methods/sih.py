"""Shear stress intensity hypothesis (SIH) of Zenner, Heidenreich & Richter, after Zenner et al. 2000 (N15).

H. Zenner, A. Simbuerger, J. Liu, Int. J. Fatigue 22 (2000) 137-145,
doi:10.1016/S0142-1123(99)00107-3, Eqs. (6)-(7) and Table 1::

    sigma_eq,a = { 15/(8 pi) int_0^pi int_0^2pi [ a tau_a^2 (1 + m tau_m^2)
                   + b sigma_a^2 (1 + n sigma_m) ] sin(gamma) dgamma dphi }^(1/2) = sigma_W      (6)
    a = (3 (sigma_W/tau_W)^2 - 4) / 5,     b = (6 - 2 (sigma_W/tau_W)^2) / 5
    a m = [sigma_W^2 - (sigma_W/tau_W)^2 (tau_Sch/2)^2] / [12/7 (tau_Sch/2)^4]
    b n = [sigma_W^2 - (sigma_Sch/2)^2 - 4/21 a m (sigma_Sch/2)^4] / [15/14 (sigma_Sch/2)^3]
    tau_Sch = 4 tau_W / (2 sigma_W / sigma_Sch + 1)                                            (7)

tau_a, tau_m, sigma_a, sigma_m are the shear and normal stress amplitudes and means on the plane
of polar angles (gamma, phi); sigma_W, tau_W are the reversed axial and torsion limits and
sigma_Sch, tau_Sch the pulsating (R = 0) limits as maximum stresses. The integral runs over the
whole sphere; 15/(8 pi) normalises it so that reversed tension gives sigma_a. a, b make the
reversed tension and torsion limits exact, m, n the pulsating ones (the factors 12/7, 4/21 and
15/14 are sphere integrals of the uniaxial and torsional plane stresses, checked in the tests).
Valid for 0.577 <= tau_W/sigma_W <= 0.866 (Section 4).

Repository choices:
  - shear amplitude on a plane: the source does not define it for a rotating shear vector (it
    cites Liu 1991 and Liu & Zenner 1993). Maintainer decision: half the longest chord of the
    shear path, i.e. the minimum circumscribed circle radius of the single-harmonic ellipse
    (``np_shear_rule = "longest_chord"``); ``"rss"`` (the quadratic amplitude) is available. With
    "rss", Eq. (6) without mean stresses is exactly independent of the phase between sigma_x and
    tau_xy, as Section 3 states for Eq. (5); with the longest chord it is not (-7 % at 90 deg for
    tau_a = sigma_a/2); the maintainer confirmed the longest chord with this evidence. tau_m is the magnitude of the mean shear vector;
  - mean stresses need sigma_Sch (workbook Sig0): without it a case with a mean stress raises
    LoadingScopeError. tau_Sch is the workbook Tau0, else Eq. (7);
  - tau_W/sigma_W outside [0.577, 0.866] raises MaterialScopeError;
  - a negative integrand sum (large compressive mean) is clipped to zero;
  - the sphere integral uses a GRID Gauss-Legendre x azimuth quadrature;
  - the model returns the index sigma_eq,a / sigma_W (1 on the fatigue limit); finite lives use the
    shared master curve on it (Section 6 notes that lives need damage sums per plane instead).

User guide: docs/method-guides/SIH.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError, MaterialScopeError
from ..invariant_methods.stress_invariants import harmonic_parts
from .bohme_papuga import _plane_components, _sphere_grid

GRID = (32, 64)
SPHERE = (
    15.0 / (8.0 * np.pi) * 4.0 * np.pi
)  # Eq. (6) factor times the sphere area (weights average)
MEAN_RTOL = 1e-12
RATIO_RANGE = (0.577, 0.866)  # tau_W / sigma_W, Section 4


def sih_constants(f_m1: float, t_m1: float, f_0=None, t_0=None) -> dict[str, float | None]:
    """a, b (Table 1) and the products a m, b n; None where the pulsating limits are missing."""
    kappa = f_m1 / t_m1
    if not RATIO_RANGE[0] <= 1.0 / kappa <= RATIO_RANGE[1]:
        raise MaterialScopeError(
            f"SIH covers 0.577 <= tau_W/sigma_W <= 0.866 (Zenner et al. 2000); got {1 / kappa:.3f}"
        )
    kappa_sq = kappa**2
    a, b = (3.0 * kappa_sq - 4.0) / 5.0, (6.0 - 2.0 * kappa_sq) / 5.0
    if t_0 is None and f_0 is not None:
        t_0 = 4.0 * t_m1 / (2.0 * f_m1 / float(f_0) + 1.0)  # Eq. (7)
    am = bn = None
    if t_0 is not None:
        th = float(t_0) / 2.0
        am = (f_m1**2 - kappa_sq * th**2) / (12.0 / 7.0 * th**4)
    if f_0 is not None and am is not None:
        sh = float(f_0) / 2.0
        bn = (f_m1**2 - sh**2 - 4.0 / 21.0 * am * sh**4) / (15.0 / 14.0 * sh**3)
    return {"a": a, "b": b, "am": am, "bn": bn}


def plane_amplitudes(S_mean, S_sin, S_cos, normals, shear_rule: str = "longest_chord"):
    """sigma_a, sigma_m, tau_a, tau_m on each normal for one harmonic."""
    sigma_m, tau_m_vec = _plane_components(np.asarray(S_mean, float), normals)
    sigma_s, tau_s = _plane_components(np.asarray(S_sin, float), normals)
    sigma_c, tau_c = _plane_components(np.asarray(S_cos, float), normals)
    ss, cc = np.sum(tau_s * tau_s, axis=1), np.sum(tau_c * tau_c, axis=1)
    if shear_rule == "longest_chord":
        sc = np.sum(tau_s * tau_c, axis=1)
        tau_a2 = np.maximum(0.0, 0.5 * (ss + cc + np.hypot(ss - cc, 2.0 * sc)))
    elif shear_rule == "rss":
        tau_a2 = ss + cc
    else:
        raise ValueError(f"unknown shear rule {shear_rule!r}")
    return np.hypot(sigma_s, sigma_c), sigma_m, np.sqrt(tau_a2), np.linalg.norm(tau_m_vec, axis=1)


def sih_index(
    S_mean, S_sin, S_cos, f_m1, t_m1, f_0=None, t_0=None, shear_rule: str = "longest_chord"
) -> float:
    """sigma_eq,a / sigma_W of Eq. (6) for one harmonic."""
    k = sih_constants(f_m1, t_m1, f_0, t_0)
    S_mean = np.asarray(S_mean, dtype=float)
    scale = max(np.max(np.abs(S_sin)), np.max(np.abs(S_cos)), np.max(np.abs(S_mean)), 1e-300)
    has_mean = np.max(np.abs(S_mean)) > MEAN_RTOL * scale
    if has_mean and k["bn"] is None:
        raise LoadingScopeError(
            "SIH mean-stress constants need the pulsating axial limit Sig0 (Table 1)"
        )
    normals, weights = _sphere_grid(*GRID)
    sigma_a, sigma_m, tau_a, tau_m = plane_amplitudes(S_mean, S_sin, S_cos, normals, shear_rule)
    integrand = k["a"] * tau_a**2 + k["b"] * sigma_a**2
    if has_mean:
        integrand = integrand + k["am"] * tau_a**2 * tau_m**2 + k["bn"] * sigma_a**2 * sigma_m
    return float(np.sqrt(max(SPHERE * np.sum(weights * integrand), 0.0)) / f_m1)


class SIH:
    """Shear stress intensity hypothesis, Eq. (6), returned as sigma_eq,a / sigma_W."""

    name = "SIH"
    np_shear_rule = "longest_chord"
    optional_params = {
        "Sig0": "required for mean stresses (Table 1)",
        "Tau0": "Eq. (7): 4 Taum1 / (2 Sigm1 / Sig0 + 1)",
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
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        value = sih_index(
            S_mean,
            S_sin,
            S_cos,
            float(params["Sigm1"]),
            float(params["Taum1"]),
            params.get("Sig0"),
            params.get("Tau0"),
            self.np_shear_rule,
        )
        return CaseResult(values={"": value})
