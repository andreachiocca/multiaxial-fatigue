"""Events: the pairs of load cases a reduction of a history hands to the damage criteria.

An event is a pair ``(a, b)`` of sample indices with a count, and the stress (and strain) ranges
between the two samples. The full tensors at both samples stay available through the history, which
is what the criteria need (the ranges alone do not carry the mean, so two cycles of equal span at
different means are different events). Events come from counted rainflow cycles or from the longest
chord of the whole history.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from multifat.fe.ranges import LongestChord

from .rainflow import Cycle
from .signal import TimeHistory


@dataclass(frozen=True, eq=False)
class Event:
    """A pair of samples with a count; ``span`` is the scalar-channel range when it was counted."""

    pair: tuple[int, int]
    count: float
    S_range: np.ndarray
    E_range: np.ndarray | None
    span: float | None


def _event(history: TimeHistory, a: int, b: int, count: float, span: float | None) -> Event:
    E = history.E_hist
    return Event(
        (int(a), int(b)),
        float(count),
        history.S_hist[a] - history.S_hist[b],
        None if E is None else E[a] - E[b],
        None if span is None else float(span),
    )


def events_from_cycles(history: TimeHistory, cycles: list[Cycle]) -> list[Event]:
    """One event per counted cycle, with the tensors at its two reversals."""
    return [_event(history, c.i_start, c.i_end, c.count, c.span) for c in cycles]


def events_from_longest_chord(
    history: TimeHistory, *, metric: str = "deviatoric_norm", on: str = "stress"
) -> list[Event]:
    """The single equivalent cycle given by the longest chord of the whole history (count 1)."""
    result = LongestChord(metric=metric, on=on)(history.S_hist, history.E_hist)
    a, b = result.pair
    return [_event(history, a, b, 1.0, None)]


def _key(history: TimeHistory, event: Event) -> tuple:
    ends = []
    for i in event.pair:
        strain = b"" if history.E_hist is None else history.E_hist[i].tobytes()
        ends.append((history.S_hist[i].tobytes(), strain))
    return tuple(sorted(ends))  # orientation (a, b) or (b, a) is the same event


def merge_identical(history: TimeHistory, events: list[Event]) -> list[Event]:
    """Merge events whose two end tensors (stress and strain) are equal, summing their counts.

    The first pair of each group is kept; the order of first appearance is preserved.
    """
    groups: dict[tuple, list[Event]] = {}
    for event in events:
        groups.setdefault(_key(history, event), []).append(event)
    merged = []
    for members in groups.values():
        first = members[0]
        merged.append(
            Event(
                first.pair, sum(e.count for e in members), first.S_range, first.E_range, first.span
            )
        )
    return merged


def check_event_count(events: list[Event], max_events: int = 2000) -> None:
    """Raise ``ValueError`` naming the limit when there are more than ``max_events`` events."""
    if len(events) > max_events:
        raise ValueError(
            f"{len(events)} distinct events exceed the limit of {max_events}; raise max_events or "
            "reduce the history (gate, longest-chord counting)"
        )
