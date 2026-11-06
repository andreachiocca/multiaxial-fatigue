"""Events: pairs of load cases with a count, from counted cycles or the longest chord."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.fe.ranges import LongestChord
from multifat.history.channels import dominant
from multifat.history.events import (
    Event,
    check_event_count,
    events_from_cycles,
    events_from_longest_chord,
    merge_identical,
)
from multifat.history.rainflow import Cycle, rainflow
from multifat.history.signal import TimeHistory


def _tension_torsion(levels, *, strain=False, shear=0.4, mean=30.0):
    """Proportional history that visits the given load-factor levels in order (piecewise linear)."""
    lam = np.concatenate(
        [np.linspace(a, b, 5, endpoint=False) for a, b in zip(levels[:-1], levels[1:])]
        + [[levels[-1]]]
    )
    S = np.zeros((len(lam), 3, 3))
    S[:, 2, 2] = mean + 100.0 * lam
    S[:, 0, 2] = S[:, 2, 0] = 100.0 * shear * lam
    E = S / 2.1e5 if strain else None
    return TimeHistory.from_arrays(np.arange(len(lam), dtype=float), S, E, stress_unit="MPa")


def test_events_carry_full_tensors_at_the_reversals():
    h = _tension_torsion([0.0, 1.0, -0.5, 0.8, -1.0, 0.2, 0.0])
    cycles = rainflow(dominant(h.S_hist))
    events = events_from_cycles(h, cycles)
    assert len(events) == len(cycles) and all(isinstance(e, Event) for e in events)
    for e, c in zip(events, cycles):
        assert e.pair == (c.i_start, c.i_end) and e.count == c.count and e.span == c.span
        assert np.array_equal(e.S_range, h.S_hist[e.pair[0]] - h.S_hist[e.pair[1]])
        assert e.E_range is None


def test_ranges_match_difference_of_end_tensors_with_strain():
    h = _tension_torsion([0.0, 1.0, -1.0, 0.5, 0.0], strain=True)
    for e in events_from_cycles(h, rainflow(dominant(h.S_hist))):
        a, b = e.pair
        assert np.array_equal(e.E_range, h.E_hist[a] - h.E_hist[b])


def test_longest_chord_event_matches_fe_longest_chord():
    h = _tension_torsion([0.0, 1.0, -1.0, 0.5, 0.0], strain=True)
    (event,) = events_from_longest_chord(h, metric="principal_spread", on="strain")
    ref = LongestChord(metric="principal_spread", on="strain")(h.S_hist, h.E_hist)
    assert event.pair == ref.pair and event.count == 1.0 and event.span is None
    assert np.array_equal(event.S_range, ref.S_range) and np.array_equal(event.E_range, ref.E_range)


def test_merge_identical_sums_counts_and_keeps_distinct_means():
    # two cycles of equal span at different means, each repeated
    h = _tension_torsion([0.0, 1.0] * 3 + [0.0] + [2.0, 3.0] * 2 + [2.0])
    events = events_from_cycles(h, rainflow(dominant(h.S_hist), residue="rearrange"))
    merged = merge_identical(h, events)
    assert len(merged) < len(events)
    assert sum(e.count for e in merged) == pytest.approx(sum(e.count for e in events))
    means = sorted(
        {round(float(np.trace(h.S_hist[e.pair[0]] + h.S_hist[e.pair[1]]) / 6.0), 6) for e in merged}
    )
    assert len(means) >= 2


def test_merge_treats_swapped_orientation_as_the_same_event():
    h = _tension_torsion([0.0, 1.0, 0.0, 1.0, 0.0])
    cycles = [Cycle(1.0, 0.5, 0.5, 0, 5), Cycle(1.0, 0.5, 0.5, 5, 10)]  # 0 -> 1 and 1 -> 0
    events = events_from_cycles(h, cycles)
    merged = merge_identical(h, events)
    assert len(merged) == 1 and merged[0].count == 1.0 and merged[0].pair == (0, 5)


def test_merge_is_idempotent_and_conserves_total_count():
    h = _tension_torsion([0.0, 1.0, -1.0, 1.0, -1.0, 0.5, -0.5, 0.5, -0.5, 0.0])
    events = events_from_cycles(h, rainflow(dominant(h.S_hist)))
    once = merge_identical(h, events)
    twice = merge_identical(h, once)
    assert [(e.pair, e.count) for e in once] == [(e.pair, e.count) for e in twice]
    assert sum(e.count for e in once) == pytest.approx(sum(e.count for e in events))


def test_event_count_guard_names_the_limit():
    h = _tension_torsion([0.0, 1.0, -1.0, 0.5, 0.0])
    events = events_from_cycles(h, rainflow(dominant(h.S_hist)))
    check_event_count(events, max_events=len(events))
    with pytest.raises(ValueError, match=str(len(events) - 1)):
        check_event_count(events, max_events=len(events) - 1)
