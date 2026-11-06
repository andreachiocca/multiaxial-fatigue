"""Range strategies: how the stress and strain ranges are taken from the load cases of a location.

A range strategy maps the ``N`` load cases ``(sigma_i, epsilon_i)`` of one location to the ranges
``X[a] - X[b]`` of a selected pair of load cases (Grossi, Chiocca and Sgamma, exact-cp, Sec. 2.2).
Peak quantities such as the maximum normal stress are maxima over all load cases and do not come
from here. The longest-chord distances are provisional definitions: their scientific review, and
the remaining rules (rainflow-based ones), belong to the variable-amplitude work (issue #20).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np


@dataclass(frozen=True)
class RangeResult:
    """Ranges ``X[a] - X[b]`` of the selected pair ``(a, b)`` of load cases."""

    S_range: np.ndarray
    E_range: np.ndarray | None
    pair: tuple[int, int]


class RangeStrategy(Protocol):
    """Callable mapping ``(N, 3, 3)`` stress (and optional strain) histories to a range."""

    def __call__(self, S_hist: np.ndarray, E_hist: np.ndarray | None) -> RangeResult: ...


def _result(S_hist, E_hist, a: int, b: int) -> RangeResult:
    E_range = None if E_hist is None else E_hist[a] - E_hist[b]
    return RangeResult(S_hist[a] - S_hist[b], E_range, (a, b))


@dataclass(frozen=True)
class ExplicitPair:
    """The user, or a cycle counter, supplies the pair of load cases."""

    a: int
    b: int

    def __call__(self, S_hist: np.ndarray, E_hist: np.ndarray | None) -> RangeResult:
        n = len(S_hist)
        for index in (self.a, self.b):
            if not 0 <= index < n:
                raise IndexError(f"load case {index} outside 0..{n - 1}")
        return _result(S_hist, E_hist, self.a, self.b)


def _distance(D: np.ndarray, metric: str) -> np.ndarray:
    """Metric of each difference tensor ``D`` of shape ``(K, 3, 3)``."""
    if metric == "deviatoric_norm":
        dev = D - np.trace(D, axis1=-2, axis2=-1)[:, None, None] / 3.0 * np.eye(3)
        return np.sqrt(np.sum(dev * dev, axis=(-2, -1)))
    w = np.linalg.eigvalsh(D)
    return w[:, -1] - w[:, 0]


@dataclass(frozen=True)
class LongestChord:
    """The pair of load cases farthest apart in the six-dimensional space of symmetric tensors.

    ``metric`` is the distance: the Frobenius norm of the deviatoric part of the difference
    (``"deviatoric_norm"``) or the spread of its principal values (``"principal_spread"``). The
    distance is taken on stress or strain (``on``), and both ranges are returned for the selected
    pair. Both distances ignore hydrostatic change, so a pure pressure cycle ties on every pair.
    Exact ties go to the lexicographically smallest pair ``(a, b)`` with ``a < b``. A single load
    case gives pair ``(0, 0)`` and a zero range. Memory is ``O(N)``, time ``O(N^2)``.
    """

    metric: Literal["deviatoric_norm", "principal_spread"] = "deviatoric_norm"
    on: Literal["stress", "strain"] = "stress"

    def __post_init__(self):
        if self.metric not in ("deviatoric_norm", "principal_spread"):
            raise ValueError(f"unknown metric {self.metric!r}")
        if self.on not in ("stress", "strain"):
            raise ValueError(f"on must be 'stress' or 'strain', got {self.on!r}")

    def __call__(self, S_hist: np.ndarray, E_hist: np.ndarray | None) -> RangeResult:
        if self.on == "strain" and E_hist is None:
            raise ValueError("a longest chord on strain needs a strain history")
        X = S_hist if self.on == "stress" else E_hist
        n = len(X)
        if n == 1:
            return _result(S_hist, E_hist, 0, 0)
        best, pair = -1.0, (0, 1)
        for a in range(n - 1):
            values = _distance(X[a] - X[a + 1 :], self.metric)
            k = int(np.argmax(values))
            if values[k] > best:
                best, pair = float(values[k]), (a, a + 1 + k)
        return _result(S_hist, E_hist, *pair)
