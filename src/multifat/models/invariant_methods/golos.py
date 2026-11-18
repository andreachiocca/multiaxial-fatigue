"""Golos, Debski & Debski (2014) modified Gough limit curve for in-phase bending-torsion (N13).

K. M. Golos, D. K. Debski, M. A. Debski, Theor. Appl. Fract. Mech. 73 (2014) 3-8,
doi:10.1016/j.tafmec.2014.07.005, Eqs. (29)-(32), for fully reversed in-phase bending (sigma) and
torsion (tau) amplitudes::

    A (sigma/f_-1)^2 + B (sigma/f_-1) + C (tau/t_-1)^2 = 1                               (29)
    A = 2/(3 sqrt3 - 4) kappa + (sqrt3 - 4)/(3 sqrt3 - 4)
    B = -2/(3 sqrt3 - 4) kappa + 2 sqrt3/(3 sqrt3 - 4),   C = 1,   kappa = f_-1/t_-1     (30)
    A = 1.67 kappa - 1.90,   B = -1.67 kappa + 2.90,   C = 1                             (31)-(32)

A + B = 1 in both forms, so the bending limit lies on the curve for every kappa; kappa = sqrt3 gives
the Gough ellipse (A = 1, B = 0, Eq. 26) and kappa = 4/3 the curve of Debski's equivalent stress
(A = 1/3, B = 2/3, Eqs. 27-28). The returned index is the left-hand side of Eq. (32) (1 on the
limit curve). The paper also states A = 1, B = 0 for kappa = 2, which Eq. (30) does not give
(A = 1.45 there): Eq. (30) is the straight line through the kappa = 4/3 and kappa = sqrt3 points.

Scope: the criterion is defined for in-phase, fully reversed bending-torsion only. A fully
reversed proportional plane stress whose principal amplitudes satisfy sigma_1 >= 0 >= sigma_2 is
exactly a bending-torsion state in a rotated frame, with sigma = sigma_1 + sigma_2 and
tau^2 = -sigma_1 sigma_2; every other case raises LoadingScopeError and the model is skipped.

Repository choices:
  - Eq. (32) with the approximated coefficients (31), as the paper's final criterion;
  - the frame-invariant reduction above (no extension beyond the source's loading scope);
  - the index is used as the damage parameter of the shared master curve DP = A N^b for
    finite-life comparisons (not part of the criterion; the index is not homogeneous in the load
    because of the linear term).

User guide: docs/method-guides/GOLOS.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError
from .stress_invariants import harmonic_parts

SQRT3 = np.sqrt(3.0)
SCOPE_TOL = 1e-9  # relative tolerance of the scope checks (mean, phase, plane stress)


def golos_coefficients(f_m1: float, t_m1: float, *, exact: bool = False) -> tuple[float, float]:
    """(A, B) of Eq. (31), or of Eq. (30) with ``exact=True``; C = 1."""
    kappa = float(f_m1) / float(t_m1)
    if exact:
        d = 3.0 * SQRT3 - 4.0
        return 2.0 / d * kappa + (SQRT3 - 4.0) / d, -2.0 / d * kappa + 2.0 * SQRT3 / d
    return 1.67 * kappa - 1.90, -1.67 * kappa + 2.90


def golos_index(sigma_a: float, tau_a: float, f_m1: float, t_m1: float) -> float:
    """Left-hand side of Eq. (32) for bending amplitude sigma_a and torsion amplitude tau_a."""
    A, B = golos_coefficients(f_m1, t_m1)
    s, t = abs(float(sigma_a)) / float(f_m1), abs(float(tau_a)) / float(t_m1)
    return A * s**2 + B * s + t**2


def bending_torsion_amplitudes(S_mean, S_sin, S_cos) -> tuple[float, float]:
    """(sigma_a, tau_a) of a fully reversed in-phase plane stress, or LoadingScopeError."""
    S_mean, S_sin, S_cos = (np.asarray(x, dtype=float) for x in (S_mean, S_sin, S_cos))
    scale = max(np.abs(S_sin).max(), np.abs(S_cos).max(), np.abs(S_mean).max(), 1e-300)
    if np.abs(S_mean).max() > SCOPE_TOL * scale:
        raise LoadingScopeError("GOLOS: Eq. (32) is defined for fully reversed loading only")
    # in phase: the two harmonics are parallel tensors (one of them may be zero)
    a, b = S_sin.ravel(), S_cos.ravel()
    aa, bb = float(np.dot(a, a)), float(np.dot(b, b))
    if aa * bb > 0.0 and (aa * bb - float(np.dot(a, b)) ** 2) > SCOPE_TOL**2 * aa * bb:
        raise LoadingScopeError("GOLOS: Eq. (32) is defined for in-phase loading only")
    principal = np.linalg.eigvalsh(S_sin if aa >= bb else S_cos)
    # plane stress: one principal amplitude is zero; the other two must have opposite signs
    zero = int(np.argmin(np.abs(principal)))
    if abs(principal[zero]) > SCOPE_TOL * scale:
        raise LoadingScopeError("GOLOS: Eq. (32) is defined for plane stress only")
    p1, p2 = np.delete(principal, zero)
    if p1 * p2 > (SCOPE_TOL * scale) ** 2:
        raise LoadingScopeError(
            "GOLOS: principal amplitudes of equal sign are not a bending-torsion state"
        )
    return abs(p1 + p2), float(np.sqrt(max(-p1 * p2, 0.0)))


class Golos:
    name = "GOLOS"

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
        sigma_a, tau_a = bending_torsion_amplitudes(*harmonic_parts(S0, S1, harmonics))
        index = golos_index(sigma_a, tau_a, float(params["Sigm1"]), float(params["Taum1"]))
        return CaseResult(values={"": index})
