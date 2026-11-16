"""Wang, Shang, Sun & Chen (2018) mean-shear-sensitive weighted-plane criterion (N17).

X.-W. Wang, D.-G. Shang, Y.-J. Sun, H. Chen, Theor. Appl. Fract. Mech. (2018),
doi:10.1016/j.tafmec.2017.10.007, Eqs. (2)-(9), (11)-(20)::

    sigma_phi = sigma/2 + sigma/2 cos 2phi + tau sin 2phi,  tau_phi = sigma/2 sin 2phi - tau cos 2phi (2)
    tau_max(t) = sqrt((sigma/2)^2 + tau^2),  phi_45(t) = 1/2 arctan(-sigma / (2 tau))              (3)
    phi_hat = sum |phi_45(t_k)| W(t_k) / sum W(t_k),  sign of sum phi_45(t_k) W(t_k)              (4)-(6)
    W(t_k) = (tau_max(t_k) - min tau_max) / (max tau_max - min tau_max)                          (7)
    tau_eq = sqrt( [tau_a + k_tau |tau_m|]^2 + 1/3 [sigma_a + k_sigma sigma_m]^2 )                (18)
    k_tau = 2 t_-1 / t_0 - 1,  k_sigma = sqrt(48 t_-1^2 / f_0^2 - 12 t_-1^2 / t_0^2) - 1          (15)
    f_0 = sqrt(2) f_-1,  t_0 = 4 t_-1 / (2 f_-1 / f_0 + 1)                                        (19)

sigma is the axial and tau the shear stress of a thin-walled tube; phi is the angle between the
plane normal and the axis, on planes perpendicular to the surface (Case A cracks). The critical
plane is the weighted mean of the maximum-shear planes over one cycle (one block), phi_hat; of the
two candidates phi_hat and phi_hat + 90 deg, the one with the larger damage is critical (Section
3.1, step 3). The life follows from tau_eq on the fully reversed torsion S-N curve (Section 4.2).

Repository choices:
  - scope: tension-torsion, the source's loading. A case with a stress other than the axial normal
    stress and the axial-hoop shear (any mean, sine or cosine part) raises LoadingScopeError;
  - one cycle is sampled at N_T equally spaced instants (the paper uses the recorded time points);
    the amplitude and mean on a plane are the exact single-harmonic values;
  - the damage of the two candidate planes is compared through tau_eq (one S-N curve, so the larger
    tau_eq is the larger damage);
  - constant tau_max (W undefined, e.g. 90 deg out of phase with tau_a/sigma_a = 0.5, Eqs. 8-9):
    the critical plane is the plane of largest maximum normal stress max_t sigma_phi(t), searched
    on N_PHI planes (Section 3.2); a zero sum of signed angles (e.g. pure tension) gives phi_hat >= 0;
  - f_0, t_0 are the workbook Sig0, Tau0, else Eq. (19); t_0 by Eq. (19) uses f_0 as given or
    estimated. Constants violating Eq. (17) (k_tau <= 0 or k_sigma <= 0, or a negative root)
    raise MaterialScopeError;
  - f_-1, t_-1 are the workbook fatigue limits Sigm1, Taum1 (the paper takes them from the S-N
    curves at N_0 = 1e6, 2e6 or 1e7); the life law is the shared master curve tau_eq = A N^b fitted
    on the fully reversed pure torsion failures (``design_loading = "torsion"``), as in the paper;
  - the variable-amplitude procedure (Wang-Brown counting, Miner, Eqs. 21-24) is not implemented.

User guide: docs/method-guides/WANG2018.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError, MaterialScopeError
from ..invariant_methods.stress_invariants import harmonic_parts

N_T = 3600  # instants per cycle for Eqs. (3)-(7)
N_PHI = 3600  # planes over [0, pi) for the constant-tau_max fallback
SCOPE_RTOL = 1e-9
FLAT_RTOL = 1e-9  # tau_max variation below this fraction of its maximum counts as constant
_T = np.linspace(0.0, 2.0 * np.pi, N_T, endpoint=False)
_SIN, _COS = np.sin(_T), np.cos(_T)


def wang2018_constants(f_m1: float, t_m1: float, f_0=None, t_0=None) -> tuple[float, float]:
    """k_tau and k_sigma of Eq. (15), with f_0 and t_0 from Eq. (19) when not given."""
    f_0 = np.sqrt(2.0) * f_m1 if f_0 is None else float(f_0)
    t_0 = 4.0 * t_m1 / (2.0 * f_m1 / f_0 + 1.0) if t_0 is None else float(t_0)
    k_tau = 2.0 * t_m1 / t_0 - 1.0
    root = 48.0 * t_m1**2 / f_0**2 - 12.0 * t_m1**2 / t_0**2
    k_sigma = np.sqrt(root) - 1.0 if root >= 0.0 else np.nan
    if not (k_tau > 0.0 and k_sigma > 0.0):
        raise MaterialScopeError(
            f"WANG2018 requires k_tau > 0 and k_sigma > 0 (Eq. 17); got k_tau = {k_tau:.3f}, "
            f"k_sigma = {k_sigma:.3f} from f_-1 = {f_m1}, t_-1 = {t_m1}, f_0 = {f_0:.1f}, "
            f"t_0 = {t_0:.1f}"
        )
    return float(k_tau), float(k_sigma)


def _plane_values(phi, sigma, tau):
    """Normal and shear stress of Eq. (2) on planes phi (any shape) for (mean, sin, cos) parts."""
    phi = np.asarray(phi, dtype=float)
    c2, s2 = np.cos(2.0 * phi), np.sin(2.0 * phi)
    normal = [s / 2.0 + s / 2.0 * c2 + t * s2 for s, t in zip(sigma, tau)]
    shear = [s / 2.0 * s2 - t * c2 for s, t in zip(sigma, tau)]
    return normal, shear


def tau_eq_on_plane(phi, sigma, tau, k_tau: float, k_sigma: float):
    """tau_eq of Eq. (18) on plane phi; sigma, tau are (mean, sin, cos) of the applied stresses."""
    (sn_m, sn_s, sn_c), (tn_m, tn_s, tn_c) = _plane_values(phi, sigma, tau)
    tau_a, sigma_a = np.hypot(tn_s, tn_c), np.hypot(sn_s, sn_c)
    return np.sqrt((tau_a + k_tau * np.abs(tn_m)) ** 2 + (sigma_a + k_sigma * sn_m) ** 2 / 3.0)


def weighted_plane_angle(sigma, tau) -> float | None:
    """phi_hat of Eqs. (3)-(7) in radians, or None when tau_max is constant over the cycle."""
    s = sigma[0] + sigma[1] * _SIN + sigma[2] * _COS
    t = tau[0] + tau[1] * _SIN + tau[2] * _COS
    tau_max = np.hypot(s / 2.0, t)
    spread = tau_max.max() - tau_max.min()
    if spread <= FLAT_RTOL * max(tau_max.max(), np.finfo(float).tiny):
        return None
    weight = (tau_max - tau_max.min()) / spread
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = 0.5 * np.arctan(-s / (2.0 * t))  # in [-pi/4, pi/4]; t = 0 gives -+pi/4
    phi = np.where(np.isnan(phi), 0.0, phi)  # sigma = tau = 0 only where tau_max (weight) is 0
    magnitude = np.sum(np.abs(phi) * weight) / np.sum(weight)
    signed = np.sum(phi * weight) / np.sum(weight)
    return float(-magnitude if signed < 0.0 else magnitude)


def wang2018_tau_eq(sigma, tau, f_m1, t_m1, f_0=None, t_0=None) -> tuple[float, float]:
    """(tau_eq, critical-plane angle phi_c in rad) for one harmonic tension-torsion cycle."""
    k_tau, k_sigma = wang2018_constants(f_m1, t_m1, f_0, t_0)
    phi_hat = weighted_plane_angle(sigma, tau)
    if phi_hat is None:
        # Section 3.2: tau_max constant, the plane of largest normal stress is critical.
        phis = np.linspace(0.0, np.pi, N_PHI, endpoint=False)
        (sn_m, sn_s, sn_c), _ = _plane_values(phis, sigma, tau)
        phi_c = float(phis[np.argmax(sn_m + np.hypot(sn_s, sn_c))])
        return float(tau_eq_on_plane(phi_c, sigma, tau, k_tau, k_sigma)), phi_c
    candidates = np.array([phi_hat, phi_hat + np.pi / 2.0])
    values = tau_eq_on_plane(candidates, sigma, tau, k_tau, k_sigma)
    i = int(np.argmax(values))
    return float(values[i]), float(candidates[i])


def tension_torsion_parts(S_mean, S_sin, S_cos):
    """(sigma, tau) as (mean, sin, cos) triples; LoadingScopeError outside tension-torsion."""
    parts = [np.asarray(S, dtype=float) for S in (S_mean, S_sin, S_cos)]
    scale = max(max(np.max(np.abs(S)) for S in parts), np.finfo(float).tiny)
    for S in parts:
        other = S.copy()
        other[1, 1] = other[0, 1] = other[1, 0] = 0.0
        if np.max(np.abs(other)) > SCOPE_RTOL * scale:
            raise LoadingScopeError(
                "WANG2018 covers tension-torsion only (axial normal stress and axial-hoop shear); "
                "the case has another stress component"
            )
    # Repository frame: y is the axis, x the hoop direction; tau_xy is the shear of Eq. (1).
    return tuple(S[1, 1] for S in parts), tuple(S[0, 1] for S in parts)


class Wang2018:
    """tau_eq of Eq. (18) on the weighted critical plane; life on the torsion S-N curve."""

    name = "WANG2018"
    design_loading = "torsion"
    optional_params = {
        "Sig0": "Eq. (19): sqrt(2) Sigm1",
        "Tau0": "Eq. (19): 4 Taum1 / (2 Sigm1 / Sig0 + 1)",
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
        sigma, tau = tension_torsion_parts(S_mean, S_sin, S_cos)
        value, _ = wang2018_tau_eq(
            sigma,
            tau,
            float(params["Sigm1"]),
            float(params["Taum1"]),
            params.get("Sig0"),
            params.get("Tau0"),
        )
        return CaseResult(values={"": value})
