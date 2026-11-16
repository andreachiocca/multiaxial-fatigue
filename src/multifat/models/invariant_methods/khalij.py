"""Khalij, Pagnacco & Troian (2017) improved Gough / Nishihara-Kawamoto criterion, OP setting (N19).

L. Khalij, E. Pagnacco, R. Troian, Int. J. Fatigue 99 (2017) 137-150,
doi:10.1016/j.ijfatigue.2017.02.016, Eqs. (12), (22), (23) and (26), zero-mean plane stress::

    D = alpha(r) [ (s_xx,a/s_-1)^2 + (s_yy,a/s_-1)^2 - (s_xx,a/s_-1)(s_yy,a/s_-1) cos phi_yy ]
        + (s_xy,a/t_-1)^2 + (1 - alpha(r)) I_1,a / s_-1 <= 1                               (23)
    I_1,a = sqrt(s_xx,a^2 + s_yy,a^2 + 2 s_xx,a s_yy,a cos phi_yy)                       (22)
    alpha(r) = 1 / (3 r^2),   r = t_-1 / s_-1                                              (26)

s_-1 and t_-1 are the fully reversed tension (or bending) and torsion limits, phi_yy the phase
shift between the two normal stresses; the phase of the shear stress does not enter (Eq. 12, the
prismatic-hull amplitude in the principal coordinate system). With alpha = 1/(3 r^2) the first two
groups are J_2,a / t_-1^2, J_2,a = (s_xx,a^2 + s_yy,a^2 - s_xx,a s_yy,a cos phi_yy)/3 + s_xy,a^2
(Eq. 12), so the criterion is frame-invariant::

    D = J_2,a / t_-1^2 + (1 - 1/(3 r^2)) I_1,a / s_-1

which is the form evaluated here. For one harmonic, J_2,a = (|dev S_sin|^2 + |dev S_cos|^2)/2 and
I_1,a = sqrt(tr(S_sin)^2 + tr(S_cos)^2); in plane stress both equal the paper's expressions. The
returned damage is D (1 on the limit), the quantity tabulated as "Damage" in Table 6.

Repository choices:
  - the OP setting alpha = 1/(3 r^2), the authors' proposal; the Gough and Nishihara-Kawamoto
    settings of Section 4.2 are not registered;
  - zero-mean loading only: the mean-stress terms of Eqs. (29)-(31) are not implemented, because the
    printed term gamma_2 s_ms / s_-1^2 is not dimensionless (s_ms is linear in the mean stresses),
    Eq. (30) for gamma_2 follows from a different base criterion, and no reading of them reproduces
    Table 7 (maintainer decision). A case with any nonzero mean stress raises LoadingScopeError;
  - plane stress as in the source: a case with an out-of-plane (z) stress component raises
    LoadingScopeError;
  - r is t_-1 / s_-1 of the material (Table 6 rows are reproduced with r from the limits of
    Table 1, not with the two-digit r printed beside them);
  - D is used as the damage parameter of the shared master curve DP = A N^b for finite-life
    comparisons (not part of the criterion; D is not homogeneous in the load because of the
    linear I_1 term).

User guide: docs/method-guides/KHALIJ.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError
from .stress_invariants import harmonic_parts

SCOPE_TOL = 1e-9  # relative tolerance of the scope checks (mean stress, plane stress)


def _deviatoric(T: np.ndarray) -> np.ndarray:
    return T - np.trace(T) / 3.0 * np.eye(3)


def khalij_alpha(f_m1: float, t_m1: float) -> float:
    """alpha(r) = 1/(3 r^2) of Eq. (26), r = t_-1/f_-1."""
    f_m1, t_m1 = float(f_m1), float(t_m1)
    if not (f_m1 > 0.0 and t_m1 > 0.0):
        raise ValueError(f"KHALIJ needs f_-1 > 0 and t_-1 > 0, got f_-1 = {f_m1}, t_m1 = {t_m1}")
    return f_m1**2 / (3.0 * t_m1**2)


def khalij_amplitudes(S_sin, S_cos) -> tuple[float, float]:
    """(J_2,a, I_1,a) of Eqs. (12) and (22) for one harmonic."""
    S_sin, S_cos = np.asarray(S_sin, float), np.asarray(S_cos, float)
    j2a = 0.5 * (np.sum(_deviatoric(S_sin) ** 2) + np.sum(_deviatoric(S_cos) ** 2))
    i1a = float(np.hypot(np.trace(S_sin), np.trace(S_cos)))
    return float(j2a), i1a


def khalij_damage(S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    """Damage D of Eq. (23) with alpha = 1/(3 r^2), zero mean stress."""
    j2a, i1a = khalij_amplitudes(S_sin, S_cos)
    return j2a / float(t_m1) ** 2 + (1.0 - khalij_alpha(f_m1, t_m1)) * i1a / float(f_m1)


def _check_scope(S_mean, S_sin, S_cos) -> None:
    scale = max(float(np.max(np.abs(S_sin))), float(np.max(np.abs(S_cos))), 1.0)
    if np.max(np.abs(S_mean)) > SCOPE_TOL * scale:
        raise LoadingScopeError(
            "KHALIJ: the source's mean-stress terms (Eqs. 29-31) are not implemented; "
            "only zero-mean loading is in scope"
        )
    out_of_plane = [T[2, :] for T in (S_sin, S_cos)] + [T[:, 2] for T in (S_sin, S_cos)]
    if max(float(np.max(np.abs(v))) for v in out_of_plane) > SCOPE_TOL * scale:
        raise LoadingScopeError("KHALIJ: plane stress (no z component) only, as in the source")


class Khalij:
    name = "KHALIJ"

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
        _check_scope(S_mean, S_sin, S_cos)
        damage = khalij_damage(S_sin, S_cos, float(params["Sigm1"]), float(params["Taum1"]))
        return CaseResult(values={"": damage})
