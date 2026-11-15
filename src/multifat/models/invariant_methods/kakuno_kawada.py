"""Kakuno & Kawada (1979) octahedral criterion for combined static and repeated stresses (G04).

H. Kakuno, Y. Kawada, Fatigue Engng Mater. Struct. 2 (1979) 229-236,
doi:10.1111/j.1460-2695.1979.tb01358.x, Eqs. (4) and (6)-(10)::

    1/3 sqrt((s_a1 - s_a2)^2 + (s_a2 - s_a3)^2 + (s_a3 - s_a1)^2)
        = A - B (s_a1 + s_a2 + s_a3) - C (s_m1 + s_m2 + s_m3)                             (4)
    A = sqrt(2/3) t_w                                                                      (7)
    B = sqrt(2)/3 (sqrt(3) t_w / s_wb - 1)                                                 (8)
    C = sqrt(2/3) (1 - k_1)/k_1 (t_w / s_wb),   k_1 = s_up / (2 s_wb)                      (9)

s_ai are the principal stress amplitudes, s_mi the principal mean stresses, s_wb and t_w the
fully reversed bending and torsion limits and s_up the pulsating (R = 0) bending limit, given as
the maximum stress of its cycle. The left-hand side of Eq. (4) is the octahedral shear stress
amplitude sqrt(2/3) sqrt(J_2,a); the sum of the principal amplitudes is the amplitude of the first
invariant. For bending sigma and torsion tau, Eq. (4) becomes Eq. (10), and Eq. (11) without mean
stresses. A torsional mean stress does not enter. The returned fatigue index is the left-hand side
plus the two normal-stress terms, divided by A (1 on the limit).

Repository choices:
  - proportional loading only, as the source's principal stress amplitudes require: a case whose
    harmonic parts S_sin and S_cos are not parallel raises LoadingScopeError. The mean stress may
    have any direction (static stresses superimposed on a repeated proportional stress);
  - the amplitude of the first invariant is |tr S_a|, positive like every amplitude, so the sign
    convention of the amplitude tensor does not matter;
  - s_wb, t_w and s_up are the workbook limits Sigm1, Taum1 and Sig0 (tension or bending as the
    workbook records them; the paper uses plane bending);
  - without Sig0, C is undefined: a case with a nonzero mean normal stress (tr S_m != 0) raises
    LoadingScopeError, while cases without it are evaluated (maintainer decision). The paper's
    average k_1 = 0.83 for structural carbon and alloy steels is not used as a default;
  - the index is used as the damage parameter of the shared master curve DP = A N^b for
    finite-life comparisons (not part of the criterion).

User guide: docs/method-guides/KAKUNO.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError
from .stress_invariants import harmonic_parts

SQRT2, SQRT3 = np.sqrt(2.0), np.sqrt(3.0)
SCOPE_TOL = 1e-9  # relative tolerance of the proportionality and mean-stress checks


def kakuno_constants(f_m1: float, t_m1: float, f_0: float | None = None):
    """(A, B, C) of Eqs. (7)-(9); C is None without the pulsating limit f_0 (maximum stress)."""
    f_m1, t_m1 = float(f_m1), float(t_m1)
    if not (f_m1 > 0.0 and t_m1 > 0.0):
        raise ValueError(f"KAKUNO needs f_-1 > 0 and t_-1 > 0, got f_-1 = {f_m1}, t_m1 = {t_m1}")
    A = np.sqrt(2.0 / 3.0) * t_m1
    B = SQRT2 / 3.0 * (SQRT3 * t_m1 / f_m1 - 1.0)
    if f_0 is None:
        return float(A), float(B), None
    k_1 = float(f_0) / (2.0 * f_m1)
    if not k_1 > 0.0:
        raise ValueError(f"KAKUNO needs f_0 > 0, got f_0 = {f_0}")
    C = np.sqrt(2.0 / 3.0) * (1.0 - k_1) / k_1 * t_m1 / f_m1
    return float(A), float(B), float(C)


def amplitude_tensor(S_sin, S_cos) -> np.ndarray:
    """Amplitude tensor S_a of a proportional harmonic; LoadingScopeError otherwise."""
    S_sin, S_cos = np.asarray(S_sin, float), np.asarray(S_cos, float)
    norm_sin, norm_cos = np.linalg.norm(S_sin), np.linalg.norm(S_cos)
    if norm_sin == 0.0 and norm_cos == 0.0:
        return np.zeros((3, 3))
    base = S_sin if norm_sin >= norm_cos else S_cos
    unit = base / np.linalg.norm(base)
    a, b = float(np.sum(S_sin * unit)), float(np.sum(S_cos * unit))
    residual = np.linalg.norm(S_sin - a * unit) + np.linalg.norm(S_cos - b * unit)
    if residual > SCOPE_TOL * max(norm_sin, norm_cos):
        raise LoadingScopeError(
            "KAKUNO: Eq. (4) uses principal stress amplitudes, defined for proportional loading "
            "only (S_sin and S_cos not parallel)"
        )
    return float(np.hypot(a, b)) * unit


def kakuno_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float, f_0: float | None = None) -> float:
    """Fatigue index of Eq. (4): (tau_oct,a + B |tr S_a| + C tr S_m) / A."""
    A, B, C = kakuno_constants(f_m1, t_m1, f_0)
    S_a = amplitude_tensor(S_sin, S_cos)
    dev = S_a - np.trace(S_a) / 3.0 * np.eye(3)
    tau_oct = np.sqrt(2.0 / 3.0) * np.sqrt(0.5 * np.sum(dev**2))
    S_mean = np.asarray(S_mean, float)
    mean_term = 0.0
    trace_mean = float(np.trace(S_mean))
    scale = max(float(np.max(np.abs(S_a))), float(np.max(np.abs(S_mean))), 1.0)
    if abs(trace_mean) > SCOPE_TOL * scale:
        if C is None:
            raise LoadingScopeError(
                "KAKUNO: a mean normal stress needs C (Eq. 9), i.e. the pulsating limit Sig0, "
                "absent from the workbook"
            )
        mean_term = C * trace_mean
    return float((tau_oct + B * abs(np.trace(S_a)) + mean_term) / A)


class KakunoKawada:
    name = "KAKUNO"
    optional_params = {
        "Sig0": "pulsating (R = 0) limit, maximum stress; needed for mean normal stresses"
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
        f_0 = params.get("Sig0")
        index = kakuno_index(
            S_mean,
            S_sin,
            S_cos,
            float(params["Sigm1"]),
            float(params["Taum1"]),
            None if f_0 is None else float(f_0),
        )
        return CaseResult(values={"": index})
