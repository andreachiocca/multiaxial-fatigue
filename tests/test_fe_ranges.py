"""Range strategies: explicit pair and longest chord over the load cases of one location."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from multifat.eval.evaluator import _apply_harmonic_overrides
from multifat.fe.ranges import ExplicitPair, LongestChord, RangeResult


def _sym(rng, shape):
    a = rng.normal(size=shape + (3, 3)) * 100.0
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def _metric(D, metric):
    if metric == "deviatoric_norm":
        dev = D - np.trace(D) / 3.0 * np.eye(3)
        return float(np.sqrt(np.sum(dev * dev)))
    w = np.linalg.eigvalsh(D)
    return float(w[-1] - w[0])


def test_explicit_pair_range_and_bounds():
    rng = np.random.default_rng(0)
    S, E = _sym(rng, (6,)), _sym(rng, (6,)) * 1e-5
    got = ExplicitPair(1, 4)(S, E)
    assert isinstance(got, RangeResult) and got.pair == (1, 4)
    assert np.array_equal(got.S_range, S[1] - S[4]) and np.array_equal(got.E_range, E[1] - E[4])
    assert ExplicitPair(2, 3)(S, None).E_range is None
    with pytest.raises(IndexError, match="load case"):
        ExplicitPair(0, 6)(S, None)


def test_longest_chord_two_states_returns_that_pair():
    rng = np.random.default_rng(1)
    S = _sym(rng, (2,))
    got = LongestChord()(S, None)
    assert got.pair == (0, 1) and np.array_equal(got.S_range, S[0] - S[1])


@pytest.mark.parametrize("metric", ["deviatoric_norm", "principal_spread"])
@pytest.mark.parametrize("on", ["stress", "strain"])
def test_longest_chord_matches_bruteforce(metric, on):
    rng = np.random.default_rng(2)
    S, E = _sym(rng, (12,)), _sym(rng, (12,)) * 1e-5
    X = S if on == "stress" else E
    best = max(
        itertools.combinations(range(12), 2), key=lambda p: _metric(X[p[0]] - X[p[1]], metric)
    )
    got = LongestChord(metric=metric, on=on)(S, E)
    assert got.pair == best
    assert np.array_equal(got.S_range, S[best[0]] - S[best[1]])
    assert np.array_equal(got.E_range, E[best[0]] - E[best[1]])


def test_single_load_case_zero_range():
    S = _sym(np.random.default_rng(3), (1,))
    got = LongestChord()(S, None)
    assert got.pair == (0, 0) and not got.S_range.any()


def test_constant_history_zero_range_no_nan():
    S = np.tile(_sym(np.random.default_rng(4), ())[None], (5, 1, 1))
    for metric in ("deviatoric_norm", "principal_spread"):
        got = LongestChord(metric=metric)(S, None)
        assert got.pair == (0, 1) and not got.S_range.any() and np.all(np.isfinite(got.S_range))


def test_tie_break_is_smallest_pair():
    A = _sym(np.random.default_rng(5), ())
    S = np.stack([0 * A, A, 0 * A, A])  # pairs (0,1), (0,3), (1,2), (2,3) tie exactly
    assert LongestChord()(S, None).pair == (0, 1)


def test_longest_chord_matches_d069_closed_form():
    rng = np.random.default_rng(6)
    s_vec, c_vec = rng.normal(size=2) * 80.0, rng.normal(size=2) * 50.0  # (xz, yz) coefficients
    t = np.arange(720) * (2.0 * np.pi / 720)
    S = np.zeros((720, 3, 3))
    for k, comp in enumerate(((0, 2), (1, 2))):
        i, j = comp
        S[:, i, j] = S[:, j, i] = s_vec[k] * np.sin(t) + c_vec[k] * np.cos(t)
    got = LongestChord()(S, None)
    span_here = np.linalg.norm(got.S_range) / np.sqrt(
        2.0
    )  # pure shear: deviatoric norm = sqrt(2) |d_tau|
    s_sin, s_cos = np.zeros((1, 3, 3)), np.zeros((1, 3, 3))
    for k, (i, j) in enumerate(((0, 2), (1, 2))):
        s_sin[0, i, j] = s_sin[0, j, i] = s_vec[k]
        s_cos[0, i, j] = s_cos[0, j, i] = c_vec[k]
    zero = np.zeros((1, 3, 3))
    T0, T1 = _apply_harmonic_overrides(zero, zero, s_sin, s_cos, shear_rule="longest_chord")
    span_d069 = np.linalg.norm((T0 - T1)[0, :2, 2])
    assert span_here == pytest.approx(span_d069, rel=1e-4)


def test_strain_missing_with_on_strain_raises():
    S = _sym(np.random.default_rng(7), (4,))
    with pytest.raises(ValueError, match="strain"):
        LongestChord(on="strain")(S, None)


def test_pure_hydrostatic_cycle_is_a_tie_on_both_metrics():
    """Both provisional distances ignore hydrostatic change: the documented blind spot."""
    S = np.stack([np.eye(3) * p for p in (0.0, 50.0, 120.0)])
    for metric in ("deviatoric_norm", "principal_spread"):
        got = LongestChord(metric=metric)(S, None)
        assert got.pair == (0, 1) and not np.any(got.S_range - got.S_range[0, 0] * np.eye(3))
        assert got.S_range[0, 0] == -50.0
