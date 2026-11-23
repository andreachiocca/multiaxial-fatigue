"""Bolchoun, Baumgartner & Kaufmann (2017): substitute amplitude and out-of-phase damage (N53).

A. Bolchoun, J. Baumgartner, H. Kaufmann, Int. J. Fatigue (2017), accepted manuscript,
doi:10.1016/j.ijfatigue.2017.01.012, Eqs. (1)-(15)::

    tau_n(t) = sigma(t) n - (n . sigma(t) n) n,   tau_n,0 = tau_n - mean(tau_n)        (4)-(6)
    F_OP = 1/2 int_t int_n |tau_n,0 x d(tau_n,0)/dt| dn dt                              (9)
    f_I = tau_max,a / sigma_max,a,   1/2 <= f_I <= 1                                     (11)
    X = 2 [X_ax (1 - f_I) + X_tors (f_I - 1/2)],   X = k, N_k or L_k                     (11)
    sigma_a,s = int_n sigma_n,a dn                                                       (12)
    D_OP(F_OP(sigma_a,s)) = sqrt(D_total^2(sigma_a,s) - D_IP^2(sigma_a,s))               (14)
    D_total,C = sqrt(D_IP,C^2 + (D_OP,C^(k_I,C / k_I,OP))^2)                             (15)

The integrals run over the unit hemisphere of plane normals (dn = cos psi dtheta dpsi).
sigma_a,s is the substitute stress amplitude in which every Woehler curve is expressed.
f_I places a cycle between pure tension (1/2) and pure torsion (1) and interpolates the slope k
and the knee point of the in-phase Woehler curve (Eq. 11). The damage of an in-phase cycle,
D_IP, is read on that curve. The extra damage of the out-of-phase behaviour, D_OP, is a function
of F_OP, the area swept by the rotating shear vectors over all planes. It is calibrated by
Eq. (14) on one experimental out-of-phase Woehler curve, by comparing its damage with the
in-phase damage at the same sigma_a,s.

Repository choices:
  - constant amplitude, single harmonic, one cycle per cycle (no counting, no Miner sum, no
    Haibach slope below the knee, D_theor = 1): N = 1 / D_total,C. For one harmonic
    sigma_n,a = sqrt(2 Var) is the normal stress amplitude, and
    F_OP = int_n pi |tau_s(n) x tau_c(n)| dn, the area of the shear ellipses (exact);
  - axial and torsion curves are the pipeline's Basquin fits on the R = -1 tension and torsion
    failures (sigma_f, b_fat, tau_f, b0), sigma_a = sigma_f (2N)^b, converted to sigma_a,s by the
    constant ratio of pure tension and pure torsion. There is no knee point: Eq. (11)
    interpolates k = -1/b and the amplitude at N_k = 1e6 (maintainer decision, 22 Nov 2026);
  - the out-of-phase curve is the Basquin fit on the most populated fully reversed out-of-phase
    tension-torsion path of the dataset (sigma_f_OP, b_OP, lambda_OP = tau_a/sigma_a,
    delta_OP; eval/material_fits.out_of_phase_fill), whose tests are design points. D_OP(F) is
    evaluated exactly on that curve at the axial amplitude giving the same F_OP (maintainer
    decision);
  - when the out-of-phase curve lies above the interpolated in-phase curve (a benefit, D_total <
    D_IP), the variant of the paper's Section 8, D_total^2 = D_IP^2 - D_OP^2, is used, decided
    at the amplitude of each cycle; a non-positive D_total^2 is an infinite life (maintainer
    decision);
  - mean stresses are out of scope, as in the paper (Section 8 leaves them open): a case with a
    mean stress raises LoadingScopeError (maintainer decision);
  - f_I below 1/2 uses the axial curve (paper); above 1, not expected by the paper, uses the
    torsion curve. The paper finds tau_max,a and sigma_max,a with the Maximum Variance Method;
    for one harmonic the largest resolved shear amplitude is the major semi-axis of the shear
    ellipse (its circumscribed-circle radius), so both are taken as the largest amplitudes over
    the GRID sphere with local refinement;
  - in the minus-sign variant, D_total,C^2 = D_IP,C^2 - D_OP,C^(2 k_I,C / k_I,OP) cancels when the
    out-of-phase curve lies far above the in-phase one: on the calibration path itself
    (k_I,C = k_I,OP) it returns the out-of-phase curve exactly, but a cycle whose f_I differs
    slightly can reach D^2 <= 0 (infinite life). This is a property of Eq. (15) with the
    minus sign, reported rather than regularised;
  - the life is native (``comparison_capabilities``); sigma_a,s is reported as the index.

User guide: docs/method-guides/BOLCHOUN.md.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..comparison import ComparisonCapabilities
from ..errors import LoadingScopeError, MissingParamsError
from ..invariant_methods.stress_invariants import harmonic_parts
from .bohme_papuga import _plane_amplitudes, _plane_components, _sphere_grid, refine_plane_max

GRID = (48, 96)
N_K = 1e6  # common knee-point life of Eq. (11) (no knee: reference life of the interpolation)
MEAN_RTOL = 1e-9
HEMISPHERE = 2.0 * np.pi  # area of the unit hemisphere: integral = 2 pi x sphere average
REQUIRED = ("sigma_f", "b_fat", "tau_f", "b0", "sigma_f_OP", "b_OP", "lambda_OP", "delta_OP")


def substitute_amplitude(S_sin, S_cos) -> float:
    """sigma_a,s of Eq. (12): hemisphere integral of the normal stress amplitude."""
    normals, weights = _sphere_grid(*GRID)
    sigma_a, _, _, _ = _plane_amplitudes(np.zeros((3, 3)), S_sin, S_cos, normals)
    return float(HEMISPHERE * np.sum(weights * sigma_a))


def out_of_phase_measure(S_sin, S_cos) -> float:
    """F_OP of Eq. (9) for one harmonic cycle: hemisphere integral of pi |tau_s x tau_c|."""
    normals, weights = _sphere_grid(*GRID)
    _, tau_s = _plane_components(np.asarray(S_sin, float), normals)
    _, tau_c = _plane_components(np.asarray(S_cos, float), normals)
    area = np.pi * np.linalg.norm(np.cross(tau_s, tau_c), axis=1)
    return float(HEMISPHERE * np.sum(weights * area))


def interpolation_factor(S_sin, S_cos) -> float:
    """f_I = tau_max,a / sigma_max,a of Eq. (11) (unclipped)."""
    zero = np.zeros((3, 3))
    normals, _ = _sphere_grid(*GRID)

    def sigma_on(n):
        return _plane_amplitudes(zero, S_sin, S_cos, n)[0]

    def tau_on(n):
        return _plane_amplitudes(zero, S_sin, S_cos, n)[2]

    sigma_max = refine_plane_max(sigma_on, normals, sigma_on(normals))
    tau_max = refine_plane_max(tau_on, normals, tau_on(normals))
    return float(tau_max / sigma_max) if sigma_max > 0 else 1.0


@lru_cache(maxsize=1)
def unit_factors() -> tuple[float, float]:
    """sigma_a,s of unit fully reversed tension and of unit fully reversed torsion."""
    tension, torsion = np.zeros((3, 3)), np.zeros((3, 3))
    tension[1, 1] = 1.0
    torsion[0, 1] = torsion[1, 0] = 1.0
    zero = np.zeros((3, 3))
    return substitute_amplitude(tension, zero), substitute_amplitude(torsion, zero)


def tension_torsion_path(sigma_a: float, lambda_OP: float, delta_OP: float):
    """(S_sin, S_cos) of sigma = sigma_a sin wt (axis y), tau = lambda sigma_a sin(wt + delta)."""
    d = np.radians(delta_OP)
    S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3))
    S_sin[1, 1] = sigma_a
    S_sin[0, 1] = S_sin[1, 0] = lambda_OP * sigma_a * np.cos(d)
    S_cos[0, 1] = S_cos[1, 0] = lambda_OP * sigma_a * np.sin(d)
    return S_sin, S_cos


class WoehlerSet:
    """Axial, torsion and out-of-phase curves in sigma_a,s and the Eq. (11) interpolation."""

    def __init__(self, params: Mapping[str, float]):
        missing = [k for k in REQUIRED if params.get(k) is None]
        if missing:
            raise MissingParamsError("BOLCHOUN", missing)
        c_ax, c_tor = unit_factors()
        # sigma_a,s = C (2N)^beta for each curve
        self.ax = (c_ax * float(params["sigma_f"]), float(params["b_fat"]))
        self.tor = (c_tor * float(params["tau_f"]), float(params["b0"]))
        self.op = (float(params["sigma_f_OP"]), float(params["b_OP"]))  # axial amplitude
        self.lambda_op, self.delta_op = float(params["lambda_OP"]), float(params["delta_OP"])
        S_sin, S_cos = tension_torsion_path(1.0, self.lambda_op, self.delta_op)
        self.c_op = substitute_amplitude(S_sin, S_cos)  # sigma_a,s per unit axial amplitude
        self.c_f = out_of_phase_measure(S_sin, S_cos)  # F_OP per unit axial amplitude squared
        self.f_op = interpolation_factor(S_sin, S_cos)
        self.k_op = self.interpolated(self.f_op)[0]

    def interpolated(self, f_I: float) -> tuple[float, float]:
        """(k, L_k) of Eq. (11) at f_I, clipped to the axial and torsion curves."""
        f = min(max(f_I, 0.5), 1.0)
        k_ax, k_tor = -1.0 / self.ax[1], -1.0 / self.tor[1]
        l_ax = self.ax[0] * (2.0 * N_K) ** self.ax[1]
        l_tor = self.tor[0] * (2.0 * N_K) ** self.tor[1]
        k = 2.0 * (k_ax * (1.0 - f) + k_tor * (f - 0.5))
        l_k = 2.0 * (l_ax * (1.0 - f) + l_tor * (f - 0.5))
        return k, l_k

    def in_phase_damage(self, amplitude: float, f_I: float) -> float:
        """D_IP = 1/N on the interpolated in-phase curve N = N_k (L/L_k)^(-k)."""
        k, l_k = self.interpolated(f_I)
        return 1.0 / (N_K * (amplitude / l_k) ** (-k))

    def out_of_phase_damage(self, F_OP: float) -> tuple[float, float]:
        """(D_OP, sign) of Eq. (14) at F_OP, on the calibration path with that F_OP."""
        sigma_a = np.sqrt(F_OP / self.c_f)
        coef, exp = self.op
        d_total = 1.0 / (0.5 * (sigma_a / coef) ** (1.0 / exp))
        d_ip = self.in_phase_damage(self.c_op * sigma_a, self.f_op)
        diff = d_total**2 - d_ip**2
        return float(np.sqrt(abs(diff))), (1.0 if diff >= 0 else -1.0)

    def life(self, amplitude: float, F_OP: float, f_I: float) -> float:
        """N = 1 / D_total,C of Eq. (15) for one constant-amplitude cycle."""
        if amplitude <= 0:
            return float("inf")
        d_ip = self.in_phase_damage(amplitude, f_I)
        if F_OP <= 0:
            return 1.0 / d_ip
        d_op, sign = self.out_of_phase_damage(F_OP)
        k_c = self.interpolated(f_I)[0]
        d2 = d_ip**2 + sign * d_op ** (2.0 * k_c / self.k_op)
        return float("inf") if d2 <= 0 else float(1.0 / np.sqrt(d2))


@lru_cache(maxsize=32)
def woehler_set(values: tuple[float, ...]) -> WoehlerSet:
    """Cached WoehlerSet for the REQUIRED constants in that order."""
    return WoehlerSet(dict(zip(REQUIRED, values)))


def bolchoun_quantities(S_mean, S_sin, S_cos) -> dict:
    """sigma_a,s, F_OP and f_I of one zero-mean harmonic cycle."""
    S_mean = np.asarray(S_mean, dtype=float)
    scale = max(np.max(np.abs(S_sin)), np.max(np.abs(S_cos)), 1e-300)
    if np.max(np.abs(S_mean)) > MEAN_RTOL * scale:
        raise LoadingScopeError(
            "BOLCHOUN has no mean-stress treatment (Section 8 of the paper): fully reversed only"
        )
    return {
        "sigma_as": substitute_amplitude(S_sin, S_cos),
        "F_OP": out_of_phase_measure(S_sin, S_cos),
        "f_I": interpolation_factor(S_sin, S_cos),
    }


class Bolchoun:
    """Substitute amplitude (Eq. 12) with out-of-phase damage (Eqs. 14-15); native life."""

    name = "BOLCHOUN"
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_comparison=False, dp_reference="none"
    )
    fitted_constants = REQUIRED

    def required_params(self) -> set[str]:
        return set(REQUIRED)

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        q = bolchoun_quantities(S_mean, S_sin, S_cos)
        return CaseResult(values={"": q["sigma_as"]}, aux={"F_OP": q["F_OP"], "f_I": q["f_I"]})

    def native_life(self, aux, params: Mapping[str, float], fit=None):
        """N from sigma_a,s (aux["dp"]), F_OP and f_I through Eqs. (11), (14), (15).

        The master-curve design rows carry no aux: they are fully reversed uniaxial tension
        tests, for which F_OP = 0 and f_I = 1/2 exactly.
        """
        missing = [k for k in REQUIRED if params.get(k) is None]
        if missing:
            raise MissingParamsError("BOLCHOUN", missing)
        curves = woehler_set(tuple(float(params[k]) for k in REQUIRED))
        return curves.life(
            float(aux["dp"]), float(aux.get("F_OP", 0.0)), float(aux.get("f_I", 0.5))
        )
