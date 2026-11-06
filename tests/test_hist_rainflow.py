"""Rainflow counting against independent oracles: the standard's example, a four-point algorithm,
Rychlik's definition of full cycles, conservation of variation and invariances."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pytest

from multifat.history.rainflow import Cycle, rainflow, turning_points

ASTM_E1049_EXAMPLE = [-2, 1, -3, 5, -1, 3, -4, 4, -2]
# The worked example of ASTM E1049: the series is the one in the MathWorks rainflow documentation
# (which cites E1049) and the counts are those reported in github.com/iamlikeme/rainflow issue 4.
# Secondary sources; the standard's own text was not available. Lee and Tjhung, "Rainflow Cycle
# Counting Techniques" (Zotero A8XFUU86) should be checked when its PDF is fetched.
ASTM_E1049_COUNTS = {3.0: 0.5, 4.0: 1.5, 6.0: 0.5, 8.0: 1.0, 9.0: 0.5}


def _multiset(cycles):
    return sorted((round(c.span, 9), round(c.mean, 9), c.count) for c in cycles)


def _by_span(cycles):
    out = defaultdict(float)
    for c in cycles:
        out[round(c.span, 9)] += c.count
    return dict(out)


def _random_reversals(rng, n):
    """Alternating reversal values: random magnitudes, alternating direction."""
    steps = rng.uniform(0.1, 10.0, size=n - 1) * np.where(np.arange(n - 1) % 2 == 0, 1.0, -1.0)
    steps *= rng.choice([-1.0, 1.0])
    return np.concatenate([[0.0], np.cumsum(steps)])


def _four_point(values):
    """Four-point rainflow (Amzallag et al. 1994) with the residue counted as half cycles."""
    stack, cycles = [], []
    for value in values:
        stack.append(value)
        while len(stack) >= 4:
            p1, p2, p3, p4 = stack[-4:]
            inner = abs(p3 - p2)
            if inner <= abs(p2 - p1) and inner <= abs(p4 - p3):
                cycles.append(Cycle(inner, (p2 + p3) / 2, 1.0, 0, 0))
                del stack[-3:-1]
            else:
                break
    for a, b in zip(stack, stack[1:]):
        cycles.append(Cycle(abs(b - a), (a + b) / 2, 0.5, 0, 0))
    return cycles


def _rychlik_periodic(ring):
    """Rainflow cycles of a periodic history by Rychlik's definition: for every local maximum, the
    pair (higher of the lowest values between it and the next higher point on each side, the maximum);
    the highest maximum, which has no higher point, pairs with the global minimum."""
    ring = np.asarray(ring, dtype=float)
    n = len(ring)
    v = np.tile(ring, 3)
    out = []
    for i in range(n, 2 * n):
        if not (v[i] > v[i - 1] and v[i] > v[i + 1]):
            continue
        right = next((j for j in range(i + 1, len(v)) if v[j] > v[i]), None)
        left = next((j for j in range(i - 1, -1, -1) if v[j] > v[i]), None)
        if right is None or left is None:
            low = ring.min()
        else:
            low = max(v[i : right + 1].min(), v[left : i + 1].min())
        out.append((round(v[i] - low, 9), round((v[i] + low) / 2, 9), 1.0))
    return sorted(out)


def test_astm_e1049_example_sequence():
    cycles = rainflow(np.array(ASTM_E1049_EXAMPLE, dtype=float))
    assert _by_span(cycles) == ASTM_E1049_COUNTS


def test_three_point_equals_independent_four_point_on_random_sequences():
    rng = np.random.default_rng(11)
    for _ in range(200):
        values = _random_reversals(rng, int(rng.integers(3, 41)))
        assert _multiset(rainflow(values)) == _multiset(_four_point(values))


def test_conservation_of_variation():
    """Twice the sum of span * count equals the total variation of the reversal sequence."""
    rng = np.random.default_rng(13)
    for _ in range(100):
        values = _random_reversals(rng, int(rng.integers(2, 40)))
        total = float(np.sum(np.abs(np.diff(values))))
        assert 2.0 * sum(c.span * c.count for c in rainflow(values)) == pytest.approx(
            total, rel=1e-12
        )


def _closed_history(rng, n=12):
    values = _random_reversals(rng, n)
    values[-1] = values[0]
    if (values[1] - values[0]) * (values[-1] - values[-2]) > 0.0:  # keep alternation at the closure
        values = values[:-1]
        values = np.concatenate([values, [values[0]]])
    return values


def test_rearranged_closed_history_has_only_full_cycles_and_closes_the_variation():
    rng = np.random.default_rng(14)
    for _ in range(60):
        x = _closed_history(rng, int(rng.integers(5, 25)))
        cycles = rainflow(x, residue="rearrange")
        assert cycles and all(c.count == 1.0 for c in cycles)
        total = float(np.sum(np.abs(np.diff(x))))
        assert 2.0 * sum(c.span for c in cycles) == pytest.approx(total, rel=1e-12)


def test_rearranged_cycles_equal_rychlik_definition_on_periodic_histories():
    rng = np.random.default_rng(12)
    for _ in range(200):
        x = _closed_history(rng, int(rng.integers(5, 31)))
        assert _multiset(rainflow(x, residue="rearrange")) == _rychlik_periodic(x[:-1])


def test_rearranged_cycle_multiset_is_invariant_to_the_starting_reversal():
    rng = np.random.default_rng(15)
    for _ in range(40):
        x = _closed_history(rng, int(rng.integers(6, 25)))
        base = _multiset(rainflow(x, residue="rearrange"))
        ring = x[:-1]
        for k in (1, 2, 3):
            rotated = np.concatenate([ring[k:], ring[:k], [ring[k]]])
            assert _multiset(rainflow(rotated, residue="rearrange")) == base


def test_indices_refer_to_original_samples():
    reversals = np.array(ASTM_E1049_EXAMPLE, dtype=float)
    x = np.concatenate(
        [np.linspace(a, b, 4, endpoint=False) for a, b in zip(reversals[:-1], reversals[1:])]
        + [[reversals[-1]]]
    )
    assert list(turning_points(x)) == [4 * k for k in range(len(reversals))]
    for c in rainflow(x):
        ends = sorted((x[c.i_start], x[c.i_end]))
        assert ends == pytest.approx([c.mean - c.span / 2, c.mean + c.span / 2])
        assert c.i_start < c.i_end


def test_constant_series_gives_no_cycles():
    assert (
        rainflow(np.full(10, 3.0)) == [] and rainflow(np.full(10, 3.0), residue="rearrange") == []
    )
    assert list(turning_points(np.full(10, 3.0))) == [0]


def test_two_samples_give_one_half_cycle():
    (c,) = rainflow(np.array([1.0, 4.0]))
    assert (c.span, c.mean, c.count, c.i_start, c.i_end) == (3.0, 2.5, 0.5, 0, 1)


def test_plateaus_and_equal_neighbours_are_collapsed():
    x = np.array([0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 1.0, 3.0, 3.0, 0.0])
    assert list(turning_points(x)) == [0, 2, 5, 7, 9]
    assert _by_span(rainflow(x)) == _by_span(rainflow(np.array([0.0, 2.0, 1.0, 3.0, 0.0])))


def test_gate_removes_small_reversals():
    x = np.array([0.0, 10.0, 9.5, 10.2, 0.0, 6.0, 5.9, 6.1, 0.0])
    kept = turning_points(x, gate=1.0)
    assert x[kept].tolist() == [0.0, 10.2, 0.0, 6.1, 0.0]
    assert turning_points(x, gate=0.0).size == 9


def test_ca_block_counts_m_cycles_in_total_and_m_full_cycles_when_rearranged():
    m = 4
    x = np.array([0.0, 1.0] * m + [0.0])
    half = rainflow(x)
    assert sum(c.count for c in half) == m and {c.span for c in half} == {1.0}
    rearranged = rainflow(x, residue="rearrange")
    assert len(rearranged) == m and all(c.count == 1.0 and c.span == 1.0 for c in rearranged)


def test_nonfinite_gate_is_rejected():
    with pytest.raises(ValueError, match="gate"):
        rainflow(np.array([0.0, 1.0, 0.0, 2.0]), gate=float("nan"))
    with pytest.raises(ValueError, match="gate"):
        turning_points(np.array([0.0, 1.0, 0.0]), gate=float("inf"))


def test_rearranged_open_history_equals_the_extra_cycles_of_one_more_tiled_block():
    """The block does not close (its end differs from its start): tile it and count with the half rule."""
    rng = np.random.default_rng(21)
    checked = 0
    for _ in range(60):
        block = _random_reversals(rng, int(rng.integers(4, 14)))

        def tally(cycles):
            out = defaultdict(float)
            for c in cycles:
                out[(round(c.span, 9), round(c.mean, 9))] += c.count
            return out

        a, b = tally(rainflow(np.tile(block, 8))), tally(rainflow(np.tile(block, 9)))
        step = {k: b.get(k, 0.0) - a.get(k, 0.0) for k in set(a) | set(b)}
        step = {k: v for k, v in step.items() if abs(v) > 1e-12}
        got = {k: v for k, v in tally(rainflow(block, residue="rearrange")).items() if v}
        assert got == pytest.approx(step)
        checked += 1
    assert checked == 60
