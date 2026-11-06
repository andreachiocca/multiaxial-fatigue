"""Basquin curve of a damage parameter and linear Miner accumulation."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from multifat.history.accumulate import Endurance, MinerResult, SNCurve, miner

CURVE = SNCurve(A=1000.0, b=-0.1, method="FS")


def test_basquin_roundtrip_dp_and_cycles():
    for n in (1e3, 1e5, 1e7):
        assert CURVE.cycles_to_failure(CURVE.dp_at(n)) == pytest.approx(n, rel=1e-9)
    assert CURVE.cycles_to_failure(np.array([500.0]))[0] == pytest.approx(0.5 ** (-10.0))


@pytest.mark.parametrize(
    "a,b",
    [(1000.0, 0.05), (1000.0, 0.0), (-1.0, -0.1), (np.nan, -0.1), (1000.0, np.nan), (np.inf, -0.1)],
)
def test_curve_rejects_flat_positive_or_nonfinite_constants(a, b):
    with pytest.raises(ValueError, match="no usable life|finite|positive"):
        SNCurve(A=a, b=b, method="FS")


def test_from_fit_uses_a_surv_and_records_provenance():
    fit = {"A": 2000.0, "A_surv": 1500.0, "b": -0.12}
    curve = SNCurve.from_fit(
        fit, method="FIN", variant="_ext", design_loading="torsion", source="test"
    )
    assert (curve.A, curve.b, curve.method, curve.variant) == (1500.0, -0.12, "FIN", "_ext")
    assert curve.design_loading == "torsion" and curve.source == "test"
    with pytest.raises(KeyError):
        SNCurve.from_fit({"A": 1.0, "b": -0.1}, method="FS")


def test_miner_equals_hand_sum_on_analytic_curve():
    dp = np.array([500.0, 250.0, 125.0])
    counts = np.array([0.5, 1.0, 2.0])
    result = miner(dp, counts, CURVE)
    assert isinstance(result, MinerResult)
    n_f = np.array([0.5**-10, 0.25**-10, 0.125**-10])
    assert np.allclose(result.n_f, n_f, rtol=1e-12)
    assert np.allclose(result.damage, counts / n_f, rtol=1e-12)
    assert result.total == pytest.approx(float(np.sum(counts / n_f)), rel=1e-12)


def test_zero_and_negative_dp_give_zero_damage_without_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = miner(np.array([0.0, -5.0, 400.0]), np.array([1.0, 1.0, 1.0]), CURVE)
    assert result.damage[0] == 0.0 and result.damage[1] == 0.0 and result.damage[2] > 0.0
    assert np.isinf(result.n_f[0]) and np.isinf(result.n_f[1])


def test_ignore_below_differs_from_extrapolate_exactly_by_the_events_below_the_limit():
    dp = np.array([800.0, 500.0, 100.0, 40.0])
    counts = np.array([1.0, 3.0, 50.0, 400.0])
    full = miner(dp, counts, CURVE, Endurance())
    cut = miner(dp, counts, CURVE, Endurance.ignore_below(120.0))
    assert np.array_equal(cut.damage[:2], full.damage[:2])
    assert np.all(cut.damage[2:] == 0.0)
    assert full.total - cut.total == pytest.approx(float(np.sum(full.damage[2:])), rel=1e-12)
    assert np.array_equal(cut.n_f, full.n_f)  # lives are still reported


def test_from_knee_limit_matches_curve_and_event_at_the_limit_counts():
    endurance = Endurance.from_knee(CURVE, 1e6)
    assert endurance.mode == "ignore_below" and endurance.dp_limit == CURVE.dp_at(1e6)
    at_limit = miner(np.array([endurance.dp_limit]), np.array([1.0]), CURVE, endurance)
    assert at_limit.damage[0] == pytest.approx(1e-6, rel=1e-9)  # dp == limit is not below it


def test_endurance_validation():
    with pytest.raises(ValueError, match="dp_limit"):
        Endurance(mode="ignore_below")
    with pytest.raises(ValueError, match="mode"):
        Endurance(mode="nonsense")
    with pytest.raises(ValueError, match="dp_limit"):
        Endurance.ignore_below(-1.0)


def test_total_damage_zero_is_finite_and_reported_as_zero():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = miner(np.array([0.0, 0.0]), np.array([1.0, 0.5]), CURVE)
    assert result.total == 0.0 and np.isfinite(result.total)


def test_shape_and_count_validation():
    with pytest.raises(ValueError, match="shape"):
        miner(np.array([1.0, 2.0]), np.array([1.0]), CURVE)
    with pytest.raises(ValueError, match="count"):
        miner(np.array([1.0]), np.array([-1.0]), CURVE)


def test_tiny_damage_parameter_gives_zero_damage_without_overflow_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = miner(np.array([1e-300, 400.0]), np.array([1.0, 1.0]), CURVE)
    assert np.isinf(result.n_f[0]) and result.damage[0] == 0.0 and result.damage[1] > 0.0
