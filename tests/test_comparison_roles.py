"""Damage and life error summaries keep design and assessment rows distinct."""

import numpy as np
import pandas as pd
import pytest

from multifat.cli import compare
from multifat.eval.metrics import log10_dp_ratio
from multifat.models.references import IMPLEMENTATION_REVISION


def _rows():
    data = [
        (0, 100.0, 10.0, 0.2),
        (0, 100.0, 100.0, 0.4),
        (1, 10.0, 100.0, -0.2),
        (1, 100.0, 0.0, -0.4),
        (2, 100.0, 1.0, 99.0),
    ]
    return pd.DataFrame(
        [
            {
                "Material": "Demo",
                "Method": "FIN",
                "Method_family": "cp_methods",
                "Metric_mode": "Nf",
                "Flag": flag,
                "Nf_exp": 1000.0,
                "Nf_expected": 900.0,
                "CP_value": observed,
                "DP_fit": expected,
                "Error_ln": life_error,
                "Implementation_revision": IMPLEMENTATION_REVISION,
            }
            for flag, expected, observed, life_error in data
        ]
    )


def test_dp_error_sign_base_and_invalid_values():
    values = log10_dp_ratio(
        np.array([100.0, 10.0, 0.0, np.inf, 10.0]),
        np.array([10.0, 100.0, 1.0, 1.0, -1.0]),
    )
    assert values[:2] == pytest.approx([1.0, -1.0])
    assert np.isnan(values[2:]).all()
    assert np.log(10.0 / 100.0) == pytest.approx(-np.log(10.0) * values[0])


def test_role_statistics_exclude_runouts_and_mask_undefined_dp():
    rows = _rows()
    unknown = rows.iloc[[0]].copy()
    unknown["Flag"] = np.nan
    stats = compare.compute_stats_table_by_role(pd.concat([rows, unknown], ignore_index=True))
    lookup = stats.set_index(["Observation_role", "Metric"])
    assert lookup.loc[("design", "Error_log10_dp"), "n"] == 2
    assert lookup.loc[("design", "Error_log10_dp"), "Error_mean"] == pytest.approx(0.5)
    assert lookup.loc[("assessment", "Error_log10_dp"), "n"] == 1
    assert lookup.loc[("assessment", "Error_log10_dp"), "Error_mean"] == pytest.approx(-1.0)
    assert lookup.loc[("assessment", "Error_ln"), "Error_mean"] == pytest.approx(-0.3)
    assert lookup.loc[("unclassified", "Error_log10_dp"), "n"] == 1
    assert stats.n.sum() == 9  # The run-out's two finite-looking errors are excluded.


def test_compare_writes_role_summary_without_changing_aggregate(tmp_path, monkeypatch):
    monkeypatch.setattr(compare, "collect_all_results", lambda path: _rows())
    monkeypatch.setattr(compare, "write_faber_style_file", lambda *args: None)
    monkeypatch.setattr(compare, "write_stats_table_with_styles", lambda **kwargs: None)
    for name in (
        "plot_boxplot_errors",
        "plot_normal_curves",
        "plot_expected_vs_experimental",
        "plot_dp_point_vs_fit",
    ):
        monkeypatch.setattr(compare, name, lambda **kwargs: None)
    out = tmp_path / "Comparison"
    compare.main(["--results_dir", str(tmp_path), "--out_dir", str(out), "--plot-kind", "damage"])
    aggregate = pd.read_csv(out / "summary_error_stats.csv")
    roles = pd.read_csv(out / "summary_error_stats_by_role.csv")
    assert set(aggregate.Metric) == {"Error_log10_dp"}
    assert aggregate.n.iloc[0] == 3
    assert set(roles.Observation_role) == {"design", "assessment"}
    assert roles.set_index("Observation_role").loc["design", "n"] == 2
    assert roles.set_index("Observation_role").loc["assessment", "n"] == 1
