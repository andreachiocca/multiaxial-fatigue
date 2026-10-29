"""Morel (1998), Fatigue Fract Engng Mater Struct 21:241-256, doi:10.1046/j.1460-2695.1998.00452.x
(roadmap item G23): mesoscopic constant-amplitude life model, provided as kernel functions.

Endurance criterion (Eqs. 9-11): T_S + a P_max <= b, where T_S is the maximum over all planes of
T_s = sqrt(int_0^2pi tau_a(psi)^2 dpsi) = sqrt(pi) T_a as defined by Papadopoulos (2001), with
a = sqrt(pi)(t - s/2)/(s/3) and b = sqrt(pi) t (Eq. 10). For a single harmonic, T_s^2 = pi (|t_s|^2 + |t_c|^2).
For the similar-loading group, Eq. (33) gives the limit loading T_Slim = f (b - a P_m)/(a + f), f = T_S/P_a.
The phase-difference coefficient is H = T_S / C_A (Eq. 35), where C_A is half the longest chord of the
shear path on the critical plane; tau_lim = T_Slim / H (Eq. 37). The life follows Eq. (38):

    N_i = p ln(C_A / (C_A - tau_lim)) + q tau_lim / (C_A - tau_lim) - r / C_A

The coefficients p, q, r are fitted to a single S-N curve in the paper (no numerical values are printed).
Because T_S / T_Slim coincides with the index of Eq. (41) in Papadopoulos (2001), the model is identical
to PAPADOPOULOS when used as a damage parameter. The MOREL case model is registered for its native life:
p, q, r are obtained on the design points (uniaxial tension, R = -1) via linear least squares on the
relative life error, since Eq. (38) is linear in p, q, r (repository choice: the paper fits them to one
S-N curve without stating an objective). Erratum: the worked example reports b = 775 MPa, whereas
Eq. (10) yields sqrt(pi) 426 = 755 MPa; only 755 reproduces the printed T_Slim = 633 MPa.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .bohme_papuga import _plane_amplitudes, _sphere_grid
from ..base import CaseResult
from ..comparison import ComparisonCapabilities
from ..errors import NativeLifeError
from ..invariant_methods.stress_invariants import harmonic_parts
from .papadopoulos import generalised_shear_amplitude

SQRT_PI = np.sqrt(np.pi)


@dataclass(frozen=True)
class MorelLimit:
    T_S: float
    C_A: float
    H: float
    P_a: float
    P_m: float
    T_Slim: float
    tau_lim: float


def morel_parameters(f_m1: float, t_m1: float) -> tuple[float, float]:
    """Return a and b of Eq. (10) from the reversed axial and torsional fatigue limits."""
    return (
        SQRT_PI * (t_m1 - f_m1 / 2.0) / (f_m1 / 3.0),
        SQRT_PI * t_m1,
    )


def _max_normal(S_sin, S_cos) -> np.ndarray:
    """Normal direction that maximises T_s: sphere grid, then pattern search from the 12 best normals."""

    def shear_amp(normal):
        return generalised_shear_amplitude(S_sin, S_cos, normal)

    normals, _ = _sphere_grid()
    current = normals[np.argsort(shear_amp(normals))[-12:]]
    moves = np.concatenate((np.eye(3), -np.eye(3), np.zeros((1, 3))))
    delta = 0.15
    for _ in range(30):
        trial = current[:, None, :] + delta * moves[None]
        trial /= np.linalg.norm(trial, axis=2, keepdims=True)
        score = shear_amp(trial.reshape(-1, 3)).reshape(len(current), -1)
        current = trial[np.arange(len(current)), np.argmax(score, axis=1)]
        delta *= 0.65
    return current[np.argmax(shear_amp(current))]


def morel_limit(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> MorelLimit:
    a, b = morel_parameters(f_m1, t_m1)
    normal = _max_normal(S_sin, S_cos)[None]
    T_S = float(SQRT_PI * generalised_shear_amplitude(S_sin, S_cos, normal)[0])
    C_A = float(_plane_amplitudes(S_mean, S_sin, S_cos, normal)[2][0])
    P_a = float(np.hypot(np.trace(S_sin), np.trace(S_cos)) / 3.0)
    P_m = float(np.trace(S_mean) / 3.0)
    f = T_S / P_a if P_a > 0.0 else np.inf
    T_Slim = (b - a * P_m) if not np.isfinite(f) else f * (b - a * P_m) / (a + f)
    H = T_S / C_A
    return MorelLimit(
        T_S=T_S,
        C_A=C_A,
        H=H,
        P_a=P_a,
        P_m=P_m,
        T_Slim=float(T_Slim),
        tau_lim=float(T_Slim / H),
    )


def morel_life(C_A, T_Slim, p: float, q: float, r: float):
    """N_i of Eq. (38); infinite when C_A <= tau_lim."""
    C_A = np.asarray(C_A, dtype=float)
    T_Slim = np.asarray(T_Slim, dtype=float)
    above = C_A > T_Slim
    gap = np.where(above, C_A - T_Slim, 1.0)
    N_i = p * np.log(np.where(above, C_A, 2.0) / gap) + q * T_Slim / gap - r / C_A
    return np.where(above, N_i, np.inf)


def fit_morel_life(C_A, T_Slim, N) -> tuple[float, float, float]:
    """p, q, r of Eq. (38) minimising the relative life error over the tests above the limit."""
    C_A, T_Slim, N = (np.asarray(v, dtype=float) for v in (C_A, T_Slim, N))
    valid = (C_A > T_Slim) & np.isfinite(N) & (N > 0)
    if valid.sum() < 3:
        raise NativeLifeError(
            f"Morel life fit needs >= 3 design tests above the limit (got {valid.sum()})"
        )
    C_A, T_Slim, N_f = C_A[valid], T_Slim[valid], N[valid]
    design = (
        np.column_stack((np.log(C_A / (C_A - T_Slim)), T_Slim / (C_A - T_Slim), -1.0 / C_A))
        / N_f[:, None]
    )
    (p, q, r), *_ = np.linalg.lstsq(design, np.ones_like(N_f), rcond=None)
    return float(p), float(q), float(r)


class Morel:
    """Morel (1998) case model: damage T_S / T_Slim (Eqs. 9-11 and 33), native life from Eq. (38)."""

    name = "MOREL"
    comparison_capabilities = ComparisonCapabilities(life_prediction="native")

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        limit = morel_limit(S_mean, S_sin, S_cos, float(params["Sigm1"]), float(params["Taum1"]))
        return CaseResult(
            values={"": limit.T_S / limit.T_Slim},
            aux={"C_A": limit.C_A, "tau_lim": limit.tau_lim},
        )

    def life_fit(self, aux, N, params: Mapping[str, float]) -> dict:
        p, q, r = fit_morel_life(
            [entry["C_A"] for entry in aux], [entry["tau_lim"] for entry in aux], N
        )
        return {"p": p, "q": q, "r": r}

    def native_life(self, aux, params: Mapping[str, float], fit):
        return morel_life(aux["C_A"], aux["tau_lim"], fit["p"], fit["q"], fit["r"])
