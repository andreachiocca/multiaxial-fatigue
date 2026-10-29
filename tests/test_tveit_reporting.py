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
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_comparison=False, dp_reference="none"
    )

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


@pytest.mark.parametrize("with_design", [False, True])
def test_native_life_export_does_not_require_a_master_curve(tmp_path, with_design):
    ev = Evaluation(
        material="Demo",
        base_name="TVEIT2025",
        params={},
        R_list=[],
        cal={
            "series": {
                "N": np.array([600.0]) if with_design else np.array([]),
                "CP": np.array([1.1]) if with_design else np.array([]),
                "CP_ext": np.array([]),
            },
            "fit": None,
            "fit_ext": None,
            "design_case_ids": [3] if with_design else [],
        },
        model=_NativeExample(),
    )
    ev.design_aux = [{"Nf": 640.0}] if with_design else []
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


def test_native_life_summary_separates_design_and_assessment_rows(tmp_path):
    rows = _rows()
    design = rows.iloc[[0]].copy()
    design["Flag"] = 0
    annotated = compare.annotate_comparison_capabilities(pd.concat([design, rows]))
    compare.write_native_life_summary(annotated, tmp_path)
    summary = pd.read_csv(tmp_path / "native_life_prediction_summary.csv")
    assert summary.loc[0, "Design_observations"] == 1
    assert summary.loc[0, "Design_points_predicted_infinite_life"] == 1
    assert summary.loc[0, "Failure_observations"] == 1
    assert summary.loc[0, "Failed_specimens_predicted_infinite_life"] == 1
    assert summary.loc[0, "Runout_observations"] == 1
    assert summary.loc[0, "Runouts_predicted_infinite_life"] == 1


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


def test_all_registered_methods_have_explicitly_verified_comparison_semantics():
    from multifat.models.registry import (
        _available_cp_methods,
        _available_direct_models,
        get_models,
    )
    from multifat.models.comparison import capabilities_for_model

    names = sorted({**_available_cp_methods(), **_available_direct_models()})
    assert len(names) == 35
    native_methods = {"TVEIT2025", "MOREL", "MDC", "MMCCM"}
    for name in names:
        model = get_models([name])[0]
        caps = capabilities_for_model(model)
        assert model.comparison_capabilities == caps
        assert caps.life_prediction == ("native" if name in native_methods else "master_curve")
        assert bool(getattr(model, "native_life", None)) == (name in native_methods)
        assert caps.dp_comparison == (name != "TVEIT2025")
        if caps.dp_reference == "model_law":
            assert callable(getattr(model, "expected_dp_at_life", None))
    assert capabilities_for_name("MOREL").dp_reference == "master_curve"
    assert capabilities_for_name("MDC").dp_reference == "model_law"
    assert capabilities_for_name("MMCCM").dp_reference == "model_law"
    with pytest.raises(ValueError, match="DP comparison"):
        ComparisonCapabilities(dp_comparison=False)
    with pytest.raises(ValueError, match="DP reference"):
        ComparisonCapabilities(dp_reference="incorrect")


@pytest.mark.parametrize(
    ("method", "params", "damage", "reference"),
    [
        (
            "MDC",
            {"sigma_f": 1200.0, "b_fat": -0.1, "tau_f": 600.0, "b0": -0.09},
            300.0,
            600.0 * (2.0 * 1000.0) ** -0.09,
        ),
        ("MMCCM", {}, 1.0 / 1200.0, 1.0 / 1000.0),
    ],
)
def test_direct_models_preserve_theoretical_dp_reference(
    method, params, damage, reference, tmp_path
):
    from multifat.models.registry import get_models

    ev = Evaluation(
        material="Demo",
        base_name=method,
        params=params,
        R_list=[],
        cal={
            "series": {
                "N": np.array([]),
                "CP": np.array([]),
                "CP_ext": np.array([]),
            },
            "fit": None,
            "fit_ext": None,
            "design_case_ids": [],
        },
        model=get_models([method])[0],
    )
    ev.values = [{"": damage}]
    ev.aux = [{}]
    ev.Nexp = [1000.0]
    ev.loading_type = ["Torsion"]
    output = pd.read_csv(write_variant(ev, "", tmp_path))
    assert output.Metric_mode.iloc[0] == "Nf_direct"
    assert output.DP_fit.iloc[0] == pytest.approx(reference)
    assert np.isfinite(output.Error_log10_dp.iloc[0])
    assert output.Nf_expected.iloc[0] == pytest.approx(output.Nf_native.iloc[0])
    annotation = compare.annotate_comparison_capabilities(output)
    assert bool(annotation.DP_comparison_supported.iloc[0])
    assert bool(annotation.Life_comparison_supported.iloc[0])


def test_morel_supports_native_life_and_independently_fitted_dp_reference(tmp_path):
    from multifat.models.registry import get_models

    model = get_models(["MOREL"])[0]
    ev = Evaluation(
        material="Demo",
        base_name="MOREL",
        params={"Sigm1": 300.0, "Taum1": 170.0},
        R_list=[],
        cal={
            "series": {
                "N": np.array([1000.0, 2000.0, 3000.0]),
                "CP": np.array([1.6, 1.4, 1.2]),
                "CP_ext": np.array([]),
            },
            "fit": {"A_surv": 3.0, "b": -0.1},
            "fit_ext": None,
            "design_case_ids": [1, 2, 3],
        },
        model=model,
    )
    ev.design_aux = [
        {"C_A": 50.0, "tau_lim": 20.0},
        {"C_A": 45.0, "tau_lim": 20.0},
        {"C_A": 40.0, "tau_lim": 20.0},
    ]
    ev.values = [{"": 1.5}]
    ev.aux = [{"C_A": 48.0, "tau_lim": 20.0}]
    ev.Nexp = [2500.0]
    ev.loading_type = ["Multiaxial proportional"]
    output = pd.read_csv(write_variant(ev, "", tmp_path))
    actual = output.loc[output.Flag == 1].iloc[0]
    assert actual.Metric_mode == "Nf_direct"
    assert actual.DP_fit == pytest.approx(3.0 * 2500.0**-0.1)
    assert np.isfinite(actual.Nf_native)
    assert actual.Nf_expected == pytest.approx(actual.Nf_native)
    support = compare.annotate_comparison_capabilities(output)
    assert support.DP_comparison_supported.all()
    assert support.Life_comparison_supported.all()
