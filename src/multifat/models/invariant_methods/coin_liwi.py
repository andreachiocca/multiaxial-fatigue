"""Equivalent stress amplitude of the COIN LiWI criterion of Schmidt & Kraft (2023).

Published in *International Journal of Fatigue* 168 (2023), 107474,
doi:10.1016/j.ijfatigue.2022.107474 (roadmap G08), single-frequency Eqs. (21)-(32)::

    sigma_eq,ka^2 = I_1^2 - k_a^2 * I_2,   k_a = s_-1 / t_-1
    H = |tr(Z)|^2 - k_a^2 * (|tr(Z)|^2 - tr(Z conj(Z))) / 2
    C = tr(Z)^2 - k_a^2 * (tr(Z)^2 - tr(Z Z)) / 2
    sigma_COIN = sqrt(k_phi * H + (1 - k_phi) * |C|)

where ``Z = S_cos - 1j * S_sin`` is the complex stress amplitude tensor of
``sigma(t) = Re(Z e^{j w t})`` (Eq. 23), ``H = 2 mean(Q)``, and ``|C| = 2 (max Q - mean Q)``
for the quadratic invariant form ``Q(t) = I_1^2 - k_a^2 I_2``.

Repository choices:
- Plane-stress Eqs. (21)-(32) extend directly to 3-D stress tensors through ``tr(Z)`` and
  ``tr(Z Z)``.
- Because no universal formula is given for ``k_phi`` (Table 3 fits 0.85-1.10 per material;
  ``1`` gives invariant RMS, ``0.5`` static maximum, ``0`` cancels 90-deg shifts),
  optional ``k_phi`` is used when provided and defaults to ``1.0`` otherwise.
- No mean-stress term is included in ``sigma_COIN`` (mean stress is handled via PRAM/PRAJ
  in the source); the normalized fatigue index ``sigma_COIN / Sigm1`` is returned.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from .stress_invariants import harmonic_parts


def coin_amplitude(S_sin, S_cos, k_a: float, k_phi: float = 1.0) -> float:
    Z = np.asarray(S_cos, dtype=float) - 1j * np.asarray(S_sin, dtype=float)
    tr_Z = np.trace(Z)
    H = abs(tr_Z) ** 2 - 0.5 * k_a**2 * (abs(tr_Z) ** 2 - np.real(np.trace(Z @ Z.conj().T)))
    C = tr_Z**2 - 0.5 * k_a**2 * (tr_Z**2 - np.trace(Z @ Z))
    return float(np.sqrt(max(k_phi * H + (1.0 - k_phi) * abs(C), 0.0)))


class CoinLiwi:
    """COIN LiWI amplitude sigma_COIN / s_-1, with k_a = s_-1/t_-1 (Eq. 22)."""

    name = "COIN"
    optional_params = {"k_phi": "1.0 (RMS formulation)"}

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
        f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
        if not (np.isfinite(f_m1) and np.isfinite(t_m1) and f_m1 > 0 and t_m1 > 0):
            raise ValueError("COIN needs Sigm1 and Taum1 to be positive and finite")
        k_phi = params.get("k_phi")
        k_phi = 1.0 if k_phi is None else float(k_phi)
        _, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": coin_amplitude(S_sin, S_cos, f_m1 / t_m1, k_phi) / f_m1})
