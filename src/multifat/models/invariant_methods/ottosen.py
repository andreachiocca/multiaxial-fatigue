"""Continuum criterion of Ottosen et al. (2008), cast in the cycle-based reduced form of Tveit et al. (2024).

Source: Tveit, Reyes & Erduran, Results in Engineering 22 (2024) 102171, Eqs. (11)-(17), p. 6::

    sigma_e,a + A * I1,m <= S0                                   (11)
    I1,m = (min_P I1 + max_P I1) / 2                             (12)
    sigma_e,a = (sum |s_i,a - s_j,a|^m / 2)^(1/m)  (Hershey-Hosford, proportional cycles)  (13)
    S0 = b_-1,  A = 2 b_-1/b_0 - 1,  kappa = t_-1/b_-1 = (1 + 2^(m-1))^(-1/m)   (14)-(17)

Repository choices: The whole left-hand side of Eq. (11) is the damage parameter
correlated with life via ``"OTT"`` in the registry. Here ``b_0`` is the R = 0 normal
fatigue strength, obtained from ``A_OTT`` when supplied, otherwise from ``Sig0``, or
from the Boehme et al. (2026) Eq. (1) estimate built from ``Su``. When ``m_OTT`` is
omitted, it is inverted from ``kappa = Taum1 / Sigm1`` and clamped to ``[1.0, 2.7670]``
(``M_CAP``). For non-proportional cycles the paper calls for the smallest hypersphere
enclosing the deviatoric path (p. 6) without giving a construction; because the
Hershey-Hosford norm is centrally symmetric, the minimal ball for a single-harmonic
history is centred at ``S_mean`` and its radius is maximized over 720 phase angles,
reproducing in-phase Table A.1 rows exactly and incremental out-of-phase rows within
about 0.8 percentage points.
"""

from __future__ import annotations

import math
from typing import Mapping

import numpy as np
from ..base import CaseResult

M_CAP = 2.7670  # Tveit p.6: the value of m for kappa >= 0.5852 (Eq. 17 peaks at kappa = 0.58524)
_PHASE_ANGLES = np.linspace(0.0, 2.0 * math.pi, 720, endpoint=False)


def kappa_of_m(m: float) -> float:
    """Tveit Eq. (17)."""
    return (1.0 + 2.0 ** (m - 1.0)) ** (-1.0 / m)


def m_of_kappa(kappa: float) -> float:
    """Inverse of Eq. (17) restricted to 1 <= m <= 2.767, the branch used in Tveit Table A.1.

    For kappa <= 0.5 the result is m = 1 (Tresca, the weak end of the model); for kappa >= kappa(2.767)
    the result is the paper's cap m = 2.7670, which lies outside the calibrated range
    0.5 <= kappa <= 0.5852.
    """
    kappa = float(kappa)
    if not math.isfinite(kappa):
        raise ValueError(f"OTT needs a finite kappa = Taum1/Sigm1, but received {kappa}")
    if kappa <= 0.5:
        return 1.0
    if kappa >= kappa_of_m(M_CAP):
        return M_CAP
    lower, upper = 1.0, M_CAP  # kappa_of_m rises monotonically along this branch
    for _ in range(60):
        middle = 0.5 * (lower + upper)
        lower, upper = (middle, upper) if kappa_of_m(middle) < kappa else (lower, middle)
    return 0.5 * (lower + upper)


def hershey_hosford(S: np.ndarray, m: float) -> float:
    """Hershey-Hosford effective stress evaluated on the deviatoric part of S (Tveit Eqs. 7, 13)."""
    stress = np.asarray(S, dtype=float)
    principal = np.linalg.eigvalsh(stress - np.trace(stress) / 3.0 * np.eye(3))
    gaps = np.abs(
        [
            principal[0] - principal[1],
            principal[1] - principal[2],
            principal[2] - principal[0],
        ]
    )
    return float((np.sum(gaps**m) / 2.0) ** (1.0 / m))


def ott_parameter(S_mean, S_sin, S_cos, m: float, A: float) -> float:
    """Value of sigma_e,a + A I1,m for the path S(t) = S_mean + S_sin sin t + S_cos cos t."""
    S_sin = np.asarray(S_sin, dtype=float)
    S_cos = np.asarray(S_cos, dtype=float)
    if np.any(S_cos != 0.0):
        sigma_e_a = max(
            hershey_hosford(S_sin * math.sin(phase) + S_cos * math.cos(phase), m)
            for phase in _PHASE_ANGLES
        )
    else:
        sigma_e_a = hershey_hosford(S_sin, m)
    # I1 varies affinely in sin t and cos t, hence its mid-range across one period is trace(S_mean).
    return float(sigma_e_a + A * np.trace(np.asarray(S_mean, dtype=float)))


class OttosenReduced:
    """Fatigue parameter of the reduced Ottosen/Tveit formulation ("OTT")."""

    name = "OTT"
    optional_params = {
        "m_OTT": "m(kappa = Taum1/Sigm1), not allowed to exceed 2.767",
        "A_OTT": "2 Sigm1/Sig0 - 1",
        "Sig0": "Boehme et al. 2026 Eq. (1), obtained from Su and Sigm1",
        "Su": "used only for estimating Sig0",
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
        if harmonics is None:
            first = np.asarray(S0, dtype=float)
            second = np.asarray(S1, dtype=float)
            S_mean, S_sin, S_cos = (
                0.5 * (first + second),
                0.5 * (first - second),
                np.zeros((3, 3)),
            )
        else:
            S_mean, S_sin, S_cos = (
                np.asarray(harmonics[key], dtype=float) for key in ("S_mean", "S_sin", "S_cos")
            )

        m = params.get("m_OTT")
        if m is None:
            m = m_of_kappa(float(params["Taum1"]) / float(params["Sigm1"]))
        A = params.get("A_OTT")
        if A is None:
            from ..integral_methods.bohme_papuga import _strengths

            f_m1_ref, _, f_0_ref, _ = _strengths(params, self.name)
            A = 2.0 * f_m1_ref / f_0_ref - 1.0
        return CaseResult(
            values={
                "": ott_parameter(
                    S_mean,
                    S_sin,
                    S_cos,
                    float(m),
                    float(A),
                )
            }
        )
