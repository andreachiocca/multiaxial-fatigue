"""Cross-batch statistics (examples/batch_comparison.py, issue #71).

The statistics are checked on synthetic errors with known structure; the last test fails when a
method or a workbook of a declared family changed since the committed results were generated.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from examples import batch_comparison as bc  # noqa: E402
from examples import batch_reference as br  # noqa: E402

BATCHES = ["B1", "B2", "B3"]
METHODS = ["M1", "M2", "M3", "M4"]


def _assessed(bias, n=12, noise=0.02, seed=0):
    """Assessed errors with a given method x batch bias and specimen noise shared by methods."""
    rng = np.random.default_rng(seed)
    rows = []
    for j, b in enumerate(BATCHES):
        shared = rng.normal(0.0, noise, n)
        shared -= shared.mean()  # the cell means are then exactly the given biases
        for k in range(n):
            for i, m in enumerate(METHODS):
                rows.append(dict(Method=m, Batch=b, Path="T", Observation_id=f"{b}-{k}",
                                 Error_load=bias[i, j] + shared[k]))  # fmt: skip
    return pd.DataFrame(rows)


def test_two_way_recovers_additive_offsets_with_intervals():
    alpha = np.array([0.1, -0.2, 0.05, 0.3])
    alpha -= alpha.mean()  # effects of an additive model are centred
    beta = np.array([0.12, 0.04, -0.16])
    summary, offsets, residuals = bc.two_way(
        _assessed(0.07 + alpha[:, None] + beta[None, :]), BATCHES
    )
    full = summary[summary["batches"] == "B1+B2+B3"].iloc[0]
    assert full["mu"] == pytest.approx(0.07)
    assert full["method_sd"] == pytest.approx(np.std(alpha, ddof=1))
    assert full["interaction_rms"] == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(residuals["residual"], 0.0)
    rows = offsets[offsets["batches"] == "B1+B2+B3"].set_index("Batch")
    for b, value in zip(BATCHES, beta):
        assert rows.loc[b, "offset"] == pytest.approx(value)
        assert rows.loc[b, "offset_lo"] < value < rows.loc[b, "offset_hi"]
        assert rows.loc[b, "offset_hi"] - rows.loc[b, "offset_lo"] < 0.06  # 12 tests, SD 0.02
    # every subset of >= 2 batches is decomposed as well
    assert set(summary["batches"]) == {"B1+B2+B3", "B1+B2", "B1+B3", "B2+B3"}


def test_two_way_interaction_is_not_sampling_noise():
    bias = np.zeros((4, 3))
    bias[0, 0] = 0.3  # one method reacts to one batch
    summary, _, _ = bc.two_way(_assessed(bias), BATCHES)
    full = summary[summary["batches"] == "B1+B2+B3"].iloc[0]
    assert full["interaction_rms"] > 0.05
    # noise shared by the methods cancels in the interaction
    assert full["noise_rms"] < 0.01
    assert full["interaction_excess"] == pytest.approx(full["interaction_rms"], rel=0.05)


def test_two_way_drops_methods_with_a_missing_test():
    data = _assessed(np.zeros((4, 3)))
    data.loc[(data["Method"] == "M4") & (data["Observation_id"] == "B2-3"), "Error_load"] = np.nan
    summary, _, _ = bc.two_way(data, BATCHES)
    assert summary[summary["batches"] == "B1+B2+B3"].iloc[0]["n_methods"] == 3
    assert summary[summary["batches"] == "B1+B3"].iloc[0]["n_methods"] == 4


def test_offsets_are_predicted_by_the_relative_strengths():
    offsets = pd.DataFrame(dict(Path="T", batches="B1+B2", Batch=["B1", "B2"], offset=[0.1, -0.1]))
    relative = pd.DataFrame(
        dict(Design=["B1", "B2", "B1"], Target=["B1", "B2", "B2"], Path="T",
             q=[-0.2, -0.5, 9.0])
    )  # fmt: skip
    out = bc.predicted_offsets(offsets, relative).set_index("Batch")
    assert out.loc["B1", "predicted"] == pytest.approx(0.15)
    assert out.loc["B2", "predicted"] == pytest.approx(-0.15)


def test_transfer_shift_of_a_degree_one_method_is_method_free():
    """bias(A -> B) - bias(A) = q(B | A) - q(A | A), whatever the method constant c_m."""
    relative = pd.DataFrame(
        dict(Design="A", Target=["A", "B", "A", "B"], Path=["T", "T", "A", "A"],
             q=[-0.20, -0.05, 0.0, 0.30])
    )  # fmt: skip
    own, moved = [], []
    for k, c in enumerate([0.0, 0.3, -0.1]):  # ln c_m of three methods
        own.append(dict(Method=f"M{k}", Batch="A", Path="T", bias=c - 0.20))
        moved.append(dict(Method=f"M{k}", Design="A", Target="B", Path="T", n=9, bias=c - 0.05))
        moved.append(dict(Method=f"M{k}", Design="A", Target="B", Path="A", n=9, bias=0.30))
    table = bc.transfer_prediction(pd.DataFrame(own), pd.DataFrame(moved), relative)
    table = table.set_index("Path")
    assert table.loc["T", "shift"] == pytest.approx(0.15)
    assert table.loc["T", "predicted"] == pytest.approx(0.15)
    assert table.loc["A", "shift"] == pytest.approx(0.30)  # design path: own bias is zero
    assert table.loc["A", "predicted"] == pytest.approx(0.30)


def test_censored_mean_recovers_the_mean_of_a_truncated_sample():
    rng = np.random.default_rng(5)
    x = rng.normal(0.2, 0.1, 4000)
    stop = 0.25  # values above it are only known to exceed it (run-outs)
    seen, bounds = x[x <= stop], np.full(int((x > stop).sum()), stop)
    assert seen.mean() < 0.17  # the failures alone are biased low
    assert bc.censored_mean(seen, bounds) == pytest.approx(0.2, abs=0.01)
    assert bc.censored_mean(seen, np.array([])) == pytest.approx(seen.mean())


def test_bootstrap_interval_is_reproducible_and_brackets_the_mean():
    x = np.random.default_rng(1).normal(0.1, 0.2, 30)
    a = bc.bootstrap_mean_ci(x, np.random.default_rng(bc.SEED))
    b = bc.bootstrap_mean_ci(x, np.random.default_rng(bc.SEED))
    assert a == b and a[0] < x.mean() < a[1]
    assert np.isnan(bc.bootstrap_mean_ci(x[:2], np.random.default_rng(0))).all()


def test_designer_factor_is_one_for_an_unbiased_method_without_scatter():
    own = pd.DataFrame(dict(Method="M", Batch="B", Flag=1, Observation_id=list("abcdef"),
                            Error_load=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))  # fmt: skip
    moved = own.assign(Design="A", Error_load=[-0.1] * 6)
    (row,) = bc.designer_summary(own, moved).to_dict("records")
    assert row["gamma_own"] == pytest.approx(1.0)
    assert row["gamma_transfer"] == pytest.approx(np.exp(0.1))  # non-conservative by 0.1


def test_basquin_strength_and_fabest_benchmark():
    N = np.array([1e4, 1e5, 1e6])
    S, b = br.basquin_strength(500.0 * (N / 1e5) ** -0.08, N)
    assert S == pytest.approx(500.0) and b == pytest.approx(-0.08)
    table, summary = br.fabest_benchmark()
    axial = summary[(summary["Load"] == "Axial") & (summary["group"] == "all unnotched series")]
    assert int(axial["series"].iloc[0]) >= 20 and int(axial["labs"].iloc[0]) >= 4
    # the PBS series of the workbook are among them, with the same strength as from the workbook
    series = br.series_strengths(["42CrMo4_FABEST"])
    for curve in ("FAB060", "FAB070"):
        a = table.set_index("Curve").loc[curve, "S_ref"]
        w = series[(series["Path"] == "A") & (series["Series"] == curve)]["S_ref"].iloc[0]
        assert a == pytest.approx(w, rel=1e-9)


def test_transfer_penalty_pairs_the_same_tests():
    own = pd.DataFrame(
        dict(Method="M", Batch="B2", Flag=[0, 1, 1, 1, 1], Observation_id=list("abcde"),
             Error_load=[0.0, 0.1, -0.1, 0.1, -0.1])
    )  # fmt: skip
    moved = pd.DataFrame(
        dict(Method="M", Batch="B2", Design="B1", Flag=1, Observation_id=list("abcdef"),
             Error_load=[0.5, 0.3, 0.3, 0.3, 0.3, 9.0])
    )  # fmt: skip
    (row,) = bc.transfer_penalty(own, moved).to_dict("records")
    # 'a' is a design test of the own run and 'f' is not in it: both are left out
    assert row["n"] == 4
    assert row["rms_own"] == pytest.approx(0.1)
    assert row["rms_transfer"] == pytest.approx(0.3)
    assert row["penalty"] == pytest.approx(0.2)


def test_committed_batch_comparison_is_current():
    """A new method, a changed workbook or a new revision needs the results regenerated."""
    problems = bc.check()
    assert not problems, (
        "; ".join(problems)
        + " -> uv run python examples/batch_comparison.py all <tmp_dir> --jobs 4"
    )


def test_variance_components_subtract_the_series_scatter():
    """s_batch^2 = s_datasets^2 - s_series^2, floored at zero (d105)."""
    series = pd.DataFrame(
        dict(Batch=["X", "Y", "Z"], Path="A", Series="all", S_ref=np.exp([0.0, 0.3, 0.6]),
             S_ref_over_Su=np.exp([0.0, 0.01, 0.02]))
    )  # fmt: skip
    yardstick = pd.DataFrame(
        dict(Load="Axial", group="all unnotched series", series=[20], sd_ln=[0.1])
    )
    out = br.variance_components(series, yardstick).set_index("scale")
    assert out.loc["raw", "datasets_sd"] == pytest.approx(0.3)
    assert out.loc["raw", "batch_sd"] == pytest.approx(np.sqrt(0.09 - 0.01))
    assert out.loc["over_Su", "batch_sd"] == 0.0  # spread below the series scatter
    assert out.loc["raw", "datasets_sd_hi"] > 3 * out.loc["raw", "datasets_sd"]  # 2 dof


@pytest.mark.parametrize(
    "p, df, expected",
    [  # reference quantiles (standard tables)
        (0.025, 1, 0.000982069),
        (0.025, 2, 0.0506356),
        (0.975, 2, 7.37776),
        (0.5, 10, 9.34182),
        (0.95, 30, 43.7730),
    ],
)
def test_chi2_quantile_without_scipy(p, df, expected):
    assert br.chi2_ppf(p, df) == pytest.approx(expected, rel=1e-5)
