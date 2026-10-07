# fatigue/models/invariant_methods/ottosen.py
from __future__ import annotations

"""Ottosen (2008) continuum approach – *reduced* cycle-based form.

This implementation intentionally follows the **cycle-based reduction** of the
Ottosen continuum model to an invariant criterion for proportional loading.

For proportional cyclic loading the model can be written as (see Eq. (11)–(13)
in Tveit et al., *Results in Engineering* 22 (2024) 102171):

    sigma_e,a + A * I1,m <= S0

Here we use the left-hand side as a scalar **fatigue parameter** that can be
correlated to fatigue life (Basquin fit) consistently with the other methods in
this repository.

Notes / limitations (by design for v5.4):
  - No incremental backstress/damage evolution is implemented.
  - Non-proportional loading is handled by sampling the harmonic stress history
    and extracting principal deviatoric stress amplitude estimates.
  - Mean-stress sensitivity A is optional; if not provided it defaults to 0.

The method name exposed through the registry is: "OTT".
"""

from typing import Mapping
import math
import numpy as np

from ..base import CaseResult


def _solve_hershey_hosford_m_from_kappa(
    kappa: float,
    *,
    m_lo: float = 1.0,
    m_hi: float = 5.0,
    it: int = 80,
) -> float:
    """Solve kappa = (1 + 2^(m-1))^(-1/m) for m by bisection.

    The paper reports the calibration interval 0.5 <= kappa <= 0.5852.
    """
    kappa = float(kappa)
    if not (math.isfinite(kappa) and kappa > 0.0):
        return 2.0

    # Clamp to the interval where the equation is meaningful for the intended model.
    if kappa <= 0.5:
        return 1.1218  # lower bound mentioned in the paper
    if kappa >= 0.5852:
        return 2.7670  # conservative upper-bound approximation (paper)

    def f(m: float) -> float:
        return (1.0 + 2.0 ** (m - 1.0)) ** (-1.0 / m) - kappa

    a, b = float(m_lo), float(m_hi)
    fa, fb = f(a), f(b)
    # Ensure bracketing; if not, fall back to m=2.
    if fa * fb > 0:
        return 2.0

    for _ in range(int(it)):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < 1e-12:
            return float(c)
        if fa * fc <= 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return float(0.5 * (a + b))


def _principal_deviatoric_eigs(S: np.ndarray) -> np.ndarray:
    """Return sorted eigenvalues (descending) of deviatoric stress tensor."""
    S = np.asarray(S, dtype=float)
    I1 = float(np.trace(S))
    dev = S - (I1 / 3.0) * np.eye(3)
    w = np.linalg.eigvalsh(dev)
    w = np.sort(w)[::-1]
    return w


def _hh_effective_from_deviatoric_eigs(w: np.ndarray, m: float) -> float:
    """Hershey–Hosford effective stress for a deviatoric tensor (Eq. (7)-type form)."""
    w1, w2, w3 = map(float, np.asarray(w, dtype=float).ravel()[:3])
    m = float(m)
    if not (math.isfinite(m) and m > 0.0):
        m = 2.0
    t1 = abs(w1 - w2) ** m
    t2 = abs(w2 - w3) ** m
    t3 = abs(w3 - w1) ** m
    return float(((t1 + t2 + t3) / 2.0) ** (1.0 / m))


def _sigma_e_amplitude_from_principal_deviatoric_amplitudes(
    s_a: np.ndarray,
    m: float,
) -> float:
    """Hershey–Hosford effective stress amplitude (Eq. (13))."""
    s1a, s2a, s3a = map(float, np.asarray(s_a, dtype=float).ravel()[:3])
    m = float(m)
    if not (math.isfinite(m) and m > 0.0):
        m = 2.0
    t1 = abs(s1a - s2a) ** m
    t2 = abs(s2a - s3a) ** m
    t3 = abs(s3a - s1a) ** m
    return float(((t1 + t2 + t3) / 2.0) ** (1.0 / m))


class OttosenReduced:
    """Ottosen reduced cycle-based fatigue parameter ("OTT")."""

    name = "OTT"

    def required_params(self) -> set[str]:
        # No strictly required parameters for the reduced form.
        # Optional:
        #   - A_OTT : mean-stress coefficient A
        #   - m_OTT : Hershey–Hosford exponent m
        #   - Sigm1, Taum1 : if m_OTT not provided, m can be derived from kappa = Taum1/Sigm1
        return set()

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
        from fatigue.models.references import ensure_enabled
        ensure_enabled("OTT")
        # ---------
        # Harmonics
        # ---------
        if harmonics is None:
            S0 = np.asarray(S0, dtype=float)
            S1 = np.asarray(S1, dtype=float)
            S_mean = 0.5 * (S0 + S1)
            S_sin = 0.5 * (S0 - S1)
            S_cos = np.zeros((3, 3), dtype=float)
        else:
            S_mean = np.asarray(harmonics.get("S_mean"), dtype=float)
            S_sin = np.asarray(harmonics.get("S_sin"), dtype=float)
            S_cos = np.asarray(harmonics.get("S_cos"), dtype=float)

        # ------------------------
        # Exponent m (optional)
        # ------------------------
        m = params.get("m_OTT", None)
        if m is None:
            Sigm1 = params.get("Sigm1", None)
            Taum1 = params.get("Taum1", None)
            if Sigm1 is not None and Taum1 is not None:
                try:
                    kappa = float(Taum1) / float(Sigm1)
                    m = _solve_hershey_hosford_m_from_kappa(kappa)
                except Exception:
                    m = 2.0
            else:
                m = 2.0
        m = float(m)

        # ------------------------
        # Mean stress coefficient A
        # ------------------------
        A = params.get("A_OTT", None)

        # Optional: if repeated bending fatigue limit b0 is available, compute A as
        # A = 2*b_-1/b0 - 1 (paper Eq. around (15)), using Sigm1 as b_-1.
        if A is None:
            b0 = params.get("b0", None)
            Sigm1 = params.get("Sigm1", None)
            if b0 is not None and Sigm1 is not None:
                try:
                    A = 2.0 * float(Sigm1) / float(b0) - 1.0
                except Exception:
                    A = 0.0
            else:
                A = 0.0
        A = float(A)

        # ------------------------
        # Sample one cycle
        # ------------------------
        # IMPORTANT:
        #   - I1,m is computed on the *full* stress history.
        #   - sigma_e,a must represent the *alternating* deviatoric part (exclude mean).
        tt = np.linspace(0.0, 2.0 * math.pi, 721)

        I1_vals = np.empty(tt.shape[0], dtype=float)
        sigma_e_alt = np.empty(tt.shape[0], dtype=float)

        for i, t in enumerate(tt):
            s = math.sin(t)
            c = math.cos(t)

            # Full stress for I1,m
            S_full = S_mean + S_sin * s + S_cos * c
            I1_vals[i] = float(np.trace(S_full))

            # Alternating part only (exclude mean)
            S_alt = S_sin * s + S_cos * c
            w_alt = _principal_deviatoric_eigs(S_alt)
            sigma_e_alt[i] = _hh_effective_from_deviatoric_eigs(w_alt, m)

        # I1 mean over the cycle (Eq. (12))
        I1_m = 0.5 * (float(np.nanmin(I1_vals)) + float(np.nanmax(I1_vals)))

        # Effective deviatoric amplitude (use maximum over the alternating cycle)
        sigma_e_a = float(np.nanmax(sigma_e_alt))

        # Fatigue parameter (left-hand side of Eq. (11))
        F = float(sigma_e_a + A * I1_m)
        if not math.isfinite(F):
            F = float("nan")
        F = max(1e-12, F)
        return CaseResult(values={"": F})
