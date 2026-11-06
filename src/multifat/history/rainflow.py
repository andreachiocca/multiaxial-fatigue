"""Scalar rainflow cycle counting (ASTM E1049 three-point form) with half-cycle residue.

``turning_points`` extracts the reversals of a sampled signal; ``rainflow`` counts them with the
three-point algorithm of Downing and Socie (1982), the form of ASTM E1049: a range ``Y`` that contains
the starting point is a half cycle, any other closed range a full cycle, and the residue left at
the end is counted as half cycles. For a history that repeats back to back, ``residue="rearrange"``
rotates the reversal sequence to start and end at its highest peak, so that no residue remains and
every cycle is a full one (exact for strictly periodic repetition). Cycles keep the sample indices of
their two reversals so that the full tensors at those instants can be recovered.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class Cycle:
    """A counted cycle: ``span`` (range), ``mean``, ``count`` (1.0 or 0.5), and the sample indices
    ``i_start < i_end`` of its two reversals (sample order; with ``residue="rearrange"`` a cycle can
    pair a reversal near the end of the block with one near its start)."""

    span: float
    mean: float
    count: float
    i_start: int
    i_end: int


def _zigzag(x: np.ndarray, gate: float) -> list[int]:
    """Positions of the reversals of ``x``; reversals smaller than ``gate`` are ignored.

    A plateau keeps its first sample. The first and last samples are always reversals (the last one
    only if it differs from the previous reversal).
    """
    out = [0]
    cand, direction = 0, 0
    for i in range(1, len(x)):
        if direction == 0:
            move = x[i] - x[out[0]]
            if move != 0.0 and abs(move) >= gate:
                direction, cand = (1 if move > 0.0 else -1), i
        elif direction > 0:
            if x[i] > x[cand]:
                cand = i
            elif x[cand] - x[i] > 0.0 and x[cand] - x[i] >= gate:
                out.append(cand)
                direction, cand = -1, i
        else:
            if x[i] < x[cand]:
                cand = i
            elif x[i] - x[cand] > 0.0 and x[i] - x[cand] >= gate:
                out.append(cand)
                direction, cand = 1, i
    if direction != 0 and cand != out[-1]:
        out.append(cand)
    return out


def turning_points(x: np.ndarray, gate: float = 0.0) -> np.ndarray:
    """Sample indices of the reversals of the signal ``x`` (endpoints included).

    Parameters
    ----------
    x : ndarray
        Sampled scalar signal.
    gate : float
        Reversals smaller than this are ignored (hysteresis filter); zero keeps every reversal.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError("x must be a non-empty 1-D signal")
    if not np.isfinite(gate) or gate < 0.0:
        raise ValueError(f"gate must be finite and not negative, got {gate}")
    return np.array(_zigzag(x, gate), dtype=int)


def _count(values: np.ndarray, indices: np.ndarray, *, closed: bool) -> list[Cycle]:
    stack: list[int] = []
    cycles: list[Cycle] = []

    def make(a: int, b: int, count: float) -> Cycle:
        lo, hi = sorted((int(indices[a]), int(indices[b])))
        return Cycle(abs(values[b] - values[a]), 0.5 * (values[a] + values[b]), count, lo, hi)

    for k in range(len(values)):
        stack.append(k)
        while len(stack) >= 3:
            a, b, c = stack[-3], stack[-2], stack[-1]
            if abs(values[c] - values[b]) < abs(values[b] - values[a]):
                break
            if len(stack) == 3 and not closed:
                cycles.append(make(a, b, 0.5))
                del stack[0]
            else:
                cycles.append(make(a, b, 1.0))
                del stack[-3:-1]
    if not closed:
        cycles.extend(make(a, b, 0.5) for a, b in zip(stack, stack[1:]))
    return cycles


def rainflow(
    x: np.ndarray, *, gate: float = 0.0, residue: Literal["half", "rearrange"] = "half"
) -> list[Cycle]:
    """Count the cycles of a sampled scalar signal.

    Parameters
    ----------
    x : ndarray
        Sampled scalar signal (a channel of a load history).
    gate : float
        Hysteresis gate of ``turning_points``.
    residue : {"half", "rearrange"}
        ``"half"``: ASTM E1049, the residue and the ranges that contain the start are half cycles.
        ``"rearrange"``: the history is taken as repeating back to back (its end connects to its
        start); the reversal sequence is rotated to start at its highest peak, so every cycle closes
        and counts one.

    Returns
    -------
    list of Cycle
        In counting order. Fewer than two reversals give an empty list.
    """
    if residue not in ("half", "rearrange"):
        raise ValueError(f"residue must be 'half' or 'rearrange', got {residue!r}")
    x = np.asarray(x, dtype=float)
    idx = turning_points(x, gate)
    if len(idx) < 2:
        return []
    values = x[idx]
    if residue == "half":
        return _count(values, idx, closed=False)
    ring = list(range(len(idx) - 1 if values[0] == values[-1] else len(idx)))
    if len(ring) < 2:
        return []
    start = int(np.argmax(values[ring]))
    order = ring[start:] + ring[:start] + [ring[start]]
    sequence = values[order]
    keep = _zigzag(sequence, 0.0)  # the wrap-around can make a reversal a point on a monotone run
    order = [order[k] for k in keep]
    if len(order) < 2:
        return []
    return _count(values[order], idx[order], closed=True)
