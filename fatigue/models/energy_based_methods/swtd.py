# fatigue/models/energy_based_methods/swtd.py
# SWTD model: Deviatoric version of SWT (SWTd), Kujawski (2014) http://dx.doi.org/10.1016/j.ijfatigue.2013.12.002

from __future__ import annotations
from typing import Mapping
import numpy as np
from ..base import CaseResult


def _sym(A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    return 0.5 * (A + A.T)


def _dev(A: np.ndarray) -> np.ndarray:
    A = _sym(A)
    return A - np.trace(A) / 3.0 * np.eye(3)


def _eig_sym(A: np.ndarray) -> np.ndarray:
    # symmetric => real eigenvalues; sorted ascending
    return np.linalg.eigvalsh(_sym(A))


class SWTD:
    name = "SWTD"

    def required_params(self) -> set[str]:
        return set()

    def evaluate_case(
        self,
        *,
        S0: np.ndarray,
        E0: np.ndarray,
        S1: np.ndarray,
        E1: np.ndarray,
        params: Mapping[str, float],
        R_list=None,
        harmonics=None,
    ) -> CaseResult:
        from fatigue.models.references import ensure_enabled
        ensure_enabled("SWTD")
        # IMPORTANT:
        # Here we assume: loadstep 0 = MAX, loadstep 1 = MIN (your stated convention).
        Smax_dev = _dev(S0)
        Smin_dev = _dev(S1)

        Emax_dev = _dev(E0)
        Emin_dev = _dev(E1)

        # deviatoric strain RANGE tensor
        dE_dev = Emax_dev - Emin_dev

        # eigenvalues
        Smax = _eig_sym(Smax_dev)      # [Smax1, Smax2, Smax3]
        Smin = _eig_sym(Smin_dev)      # [Smin1, Smin2, Smin3]
        dE = _eig_sym(dE_dev)          # [ΔE1, ΔE2, ΔE3]

        # strain amplitude principal values are |ΔEi|/2
        ea = 0.5 * np.abs(dE)

        # Eq.(25): MAX of 6 candidates with same index i
        candidates = np.concatenate([Smax * ea, Smin * ea])
        swtd = float(np.max(candidates))

        return CaseResult(values={"": swtd})