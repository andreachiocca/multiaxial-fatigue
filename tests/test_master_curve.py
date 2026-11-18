"""Master curve DP = A N^b: unbiased scatter and identifiability (issue #74)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multifat.cli import pipeline
from multifat.cli.run import check_fit
from multifat.eval.fitting import MIN_DECADES, MIN_POINTS, fit_power_law_with_survival_std

DATA = Path(__file__).resolve().parents[1] / "Experimental_Data"


def test_scatter_has_n_minus_2_degrees_of_freedom():
    """Mean of STD_ln_DP^2 over many fits equals the true variance (two fitted parameters)."""
    rng = np.random.default_rng(74)
    sigma, n = 0.1, 5
    N = np.logspace(4, 6, n)
    variances = [
        fit_power_law_with_survival_std(N, 1e3 * N**-0.1 * np.exp(rng.normal(0, sigma, n)))[
            "STD_ln_DP"
        ]
        ** 2
        for _ in range(4000)
    ]
    assert np.mean(variances) / sigma**2 == pytest.approx(1.0, abs=0.05)


def test_two_points_have_no_scatter_and_no_identifiable_curve():
    fit = fit_power_law_with_survival_std([1e5, 1e6], [300.0, 250.0])
    assert np.isnan(fit["STD_ln_DP"]) and np.isnan(fit["b_se"])
    assert fit["A_surv"] == pytest.approx(fit["A"])  # median curve, no shift needed
    assert not fit["identifiable"] and "2 < 5 design tests" in fit["reason"]


def test_identifiable_curve():
    N = np.logspace(4, 6, 8)
    fit = fit_power_law_with_survival_std(N, 1e3 * N**-0.1 * (1 + 0.01 * np.sin(np.arange(8))))
    assert fit["identifiable"] and fit["reason"] == ""
    assert fit["n"] >= MIN_POINTS and fit["decades"] >= MIN_DECADES
    assert fit["b"] == pytest.approx(-0.1, abs=0.01) and fit["b_se"] < 0.01


@pytest.mark.parametrize(
    "N, DP, expected",
    [
        (np.logspace(4.5, 4.9, 6), 1e3 * np.logspace(4.5, 4.9, 6) ** -0.1, "decades of life"),
        (np.logspace(4, 6, 6), np.array([1.0, 1.1, 0.9, 1.05, 0.95, 1.0]), "not distinct from 0"),
        (np.logspace(4, 6, 6), 1e3 * np.logspace(4, 6, 6) ** 0.1, "not distinct from 0"),
    ],
)
def test_curves_that_are_not_identifiable(N, DP, expected):
    fit = fit_power_law_with_survival_std(N, DP)
    assert not fit["identifiable"] and expected in fit["reason"]
    with pytest.raises(ValueError, match="not identifiable"):
        check_fit("X", "", {"fit": {**fit, "A_surv": 1.0, "b": -0.1}, "fit_ext": None})


def test_pipeline_skips_and_reports_a_curve_that_is_not_identifiable(tmp_path):
    """Ti6Al4V_Meneghetti: 6 axial design tests over 0.39 decades of life."""
    outcomes = pipeline.run_workbook(
        DATA / "1_Ti6Al4V_Meneghetti.xlsx", models=["CROSSLAND"], results_dir=tmp_path
    )
    assert [o.status for o in outcomes] == ["skipped"]
    assert "master curve not identifiable" in outcomes[0].reason
    table = pd.read_csv(tmp_path / "Ti6Al4V_Meneghetti_master_curves.csv")
    row = table.set_index("model").loc["CROSSLAND"]
    assert row["status"] == "skipped" and not row["identifiable"]
    assert row["n"] == 6 and row["decades"] == pytest.approx(0.389, abs=1e-3)
    assert not (tmp_path / "CROSSLAND_Ti6Al4V_Meneghetti.csv").exists()
