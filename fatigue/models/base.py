# fatigue/models/base.py
from __future__ import annotations

"""Shared interfaces for fatigue models.

This project distinguishes two levels:

1) **Case-level models** ("models")
   - Consume tensors for two load steps (S0/E0 and S1/E1)
   - Optionally consume a list of candidate plane rotations (R_list)
   - Return one or more scalar variants through :class:`CaseResult`

2) **Plane-level critical-plane methods** ("CP methods")
   - Consume already-rotated tensors for ONE candidate plane
   - Return a damage value and a plane-selection metric through :class:`CPPlaneResult`

CP methods are wrapped into case-level models via :class:`fatigue.models.cp_adapter.CPModelAdapter`.
"""

from dataclasses import dataclass
from typing import Mapping, Protocol

import numpy as np


# -------------------------
# Case-level model interface
# -------------------------


@dataclass(frozen=True)
class CaseResult:
    """Returned by a case-level model.

    The `values` dictionary stores one or more variants for the same model.

    Convention used across the codebase:
      - ""     : standard value
      - "_ext" : extreme/max value (optional)
    """

    values: dict[str, float]


class FatigueModel(Protocol):
    """Case-level model protocol."""

    name: str

    def required_params(self) -> set[str]:
        """Material parameters required by the model."""

        ...

    def evaluate_case(
        self,
        *,
        S0: np.ndarray,
        E0: np.ndarray,
        S1: np.ndarray,
        E1: np.ndarray,
        params: Mapping[str, float],
        R_list: np.ndarray | None = None,
        harmonics: Mapping[str, np.ndarray] | None = None,
    ) -> CaseResult:
        """Evaluate the model on a single load case."""

        ...


# ---------------------------------
# Plane-level critical-plane methods
# ---------------------------------


@dataclass(frozen=True)
class CPPlaneResult:
    """Returned by a critical-plane method evaluated on a single rotated plane."""

    damage: float  # value to correlate with fatigue life
    metric: float  # value used to select the plane (max metric plane)


class CriticalPlaneMethod(Protocol):
    """Plane-level critical plane method protocol."""

    name: str

    def required_params(self) -> set[str]:
        """Material parameters required by this method (e.g. Sy, k_FS, k_FI)."""

        ...

    def evaluate_on_plane(
        self,
        *,
        S0r: np.ndarray,
        E0r: np.ndarray,
        S1r: np.ndarray,
        E1r: np.ndarray,
        params: Mapping[str, float],
    ) -> CPPlaneResult:
        """Evaluate the method on ONE rotated candidate plane."""

        ...
