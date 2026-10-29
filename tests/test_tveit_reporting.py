"""Comparison capabilities select valid life and DP outputs per method."""

import numpy as np
import pandas as pd
import pytest

from multifat.cli import compare
from multifat.cli.run import Evaluation, write_variant
from multifat.models.comparison import (
    ComparisonCapabilities,
    capabilities_for_name,
)
from multifat.models.references import IMPLEMENTATION_REVISION


class _NativeExample:
    name = "TVEIT2025"
    comparison_capabilities = ComparisonCapabilities(life_prediction="native", dp_comparison=False)

    def life_fit(self, aux, N, params):
        return {}

    def native_life(self, aux, params, fit):
        return aux["Nf"]


def _rows():
    return pd.DataFrame(
        [
            {
                "Material": "Demo",
                "Method": "TVEIT2025",
                "Method_family": "invariant_methods",
                "Metric_mode": "Nf_direct",
                "Flag": 1,
                "Nf_exp": 1000.0,
                "Nf_expected": np.inf,
                "Nf_native": np.inf,
                "Error_ln": np.nan,
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 1.2,
                "DP_fit": np.nan,
            },
            {
                "Material": "Demo",
                "Method": "TVEIT2025",
                "Method_family": "invariant_methods",
                "Metric_mode": "Nf_direct",
                "Flag": 2,
                "Nf_exp": 1000000.0,
                "Nf_expected": np.inf,
                "Nf_native": np.inf,
                "Error_ln": np.nan,
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 0.8,
                "DP_fit": np.nan,
            },
            {
                "Material": "Demo",
                "Method": "FIN",
                "Method_family": "cp_methods",
                "Metric_mode": "Nf",
                "Flag": 1,
                "Nf_exp": 1000.0,
                "Nf_expected": 900.0,
                "Nf_native": np.nan,
                "Error_ln": np.log(1000.0 / 900.0),
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 100.0,
                "DP_fit": 110.0,
            },
        ]
    )


def test_declared_method_capabilities():
    tveit = capabilities_for_name("TVEIT2025")
    assert tveit.life_prediction == "native"
    assert not tveit.dp_comparison
    assert capabilities_for_name("FIN").life_prediction == "master_curve"
    assert capabilities_for_name("FIN_ext").dp_comparison
    with pytest.raises(ValueError, match="life prediction mechanism"):
        ComparisonCapabilities(life_prediction="unsupported")


def test_old_tveit_auxiliary_index_not_accepted_for_life_or_dp():
    legacy = _rows().iloc[[0]].copy()
    legacy["Metric_mode"] = "Nf"
    support = compare.annotate_comparison_capabilities(legacy)
    assert not support.iloc[0]["Life_comparison_supported"]
    assert not support.iloc[0]["DP_comparison_supported"]


def test_native_life_export_does_not_require_a_master_curve(tmp_path):
    ev = Evaluation(
        material="Demo",
        base_name="TVEIT2025",
        params={},
        R_list=[],
        cal={
            "series": {
                "N": np.array([600.0]),
                "CP": np.array([1.1]),
                "CP_ext": np.array([]),
            },
            "fit": None,
            "fit_ext": None,
            "design_case_ids": [3],
        },
        model=_NativeExample(),
    )
    ev.design_aux = [{"Nf": 640.0}]
    ev.values = [{"": 1.4}, {"": 1.2}]
    ev.aux = [{"Nf": np.inf}, {"Nf": 1200.0}]
    ev.Nexp = [1000.0, 1000.0]
    ev.loading_type = ["Torsion", "Tensile"]
    ev.runout["values"] = [{"": 0.8}]
    ev.runout["aux"] = [{"Nf": np.inf}]
    ev.runout["Nexp"] = [2000000.0]
    ev.runout["lt"] = ["Tensile"]
    output = pd.read_csv(write_variant(ev, "", tmp_path))
    assert set(output.Metric_mode) == {"Nf_direct"}
    assert output.DP_fit.isna().all()
    assert np.isinf(output.loc[output.Flag == 1, "Nf_expected"]).sum() == 1
    assert np.isinf(output.loc[output.Flag == 2, "Nf_native"]).all()
    assert output.loc[output.Flag == 1, "Error_log10_dp"].isna().all()


@pytest.mark.parametrize(
    ("kind", "life_methods", "damage_methods"),
    [
        ("life", {"TVEIT2025", "FIN"}, set()),
        ("damage", set(), {"FIN"}),
        ("both", {"TVEIT2025", "FIN"}, {"FIN"}),
    ],
)
def test_plot_types_respect_capabilities(tmp_path, monkeypatch, kind, life_methods, damage_methods):
    seen = {"life": [], "damage": []}
    monkeypatch.setattr(compare, "collect_all_results", lambda path: _rows())
    monkeypatch.setattr(compare, "write_faber_style_file", lambda *args: None)
    monkeypatch.setattr(compare, "write_stats_table_with_styles", lambda **kwargs: None)

    def plot_life(**kwargs):
        seen["life"].append(set(kwargs["df_mat"].Method))

    def plot_damage(**kwargs):
        seen["damage"].append(set(kwargs["df_mat"].Method))

    monkeypatch.setattr(compare, "plot_expected_vs_experimental", plot_life)
    monkeypatch.setattr(compare, "plot_dp_point_vs_fit", plot_damage)
    monkeypatch.setattr(compare, "plot_boxplot_errors", lambda **kwargs: None)
    monkeypatch.setattr(compare, "plot_normal_curves", lambda **kwargs: None)
    compare.main(
        [
            "--results_dir",
            str(tmp_path),
            "--out_dir",
            str(tmp_path / "Comparison"),
            "--plot-kind",
            kind,
        ]
    )
    assert (set.union(*seen["life"]) if seen["life"] else set()) == life_methods
    assert (set.union(*seen["damage"]) if seen["damage"] else set()) == damage_methods
    counts = pd.read_csv(tmp_path / "Comparison" / "native_life_prediction_summary.csv")
    assert counts.loc[0, "Failed_specimens_predicted_infinite_life"] == 1
    assert counts.loc[0, "Runouts_predicted_infinite_life"] == 1
    stats = pd.read_csv(tmp_path / "Comparison" / "summary_error_stats.csv")
    assert set(stats.Metric) == (
        {"Error_ln"}
        if kind == "life"
        else {"Error_log10_dp"} if kind == "damage" else {"Error_ln", "Error_log10_dp"}
    )
    assert not ((stats.Method == "TVEIT2025") & (stats.Metric == "Error_log10_dp")).any()


def test_damage_common_subset_rejects_native_only_method(tmp_path, monkeypatch):
    monkeypatch.setattr(compare, "collect_all_results", lambda path: _rows())
    with pytest.raises(ValueError, match="valid damage comparison"):
        compare.main(
            [
                "--results_dir",
                str(tmp_path),
                "--out_dir",
                str(tmp_path / "Comparison"),
                "--comparison-metric",
                "damage",
                "--common-methods",
                "FIN",
                "TVEIT2025",
            ]
        )


def test_invalid_native_method_has_no_comparison_in_legacy_results():
    invalid = _rows().iloc[[0]].copy()
    invalid["Metric_mode"] = "Nf"
    selected = compare.annotate_comparison_capabilities(invalid)
    assert selected["Life_comparison_supported"].eq(False).all()
    assert selected["DP_comparison_supported"].eq(False).all()
