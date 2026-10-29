"""Regression coverage for native-life methods pending a dedicated comparison route."""

import numpy as np
import pandas as pd
import pytest

from multifat.cli import compare
from multifat.models.references import IMPLEMENTATION_REVISION


def _rows():
    return pd.DataFrame(
        [
            {
                "Material": "Demo",
                "Method": "TVEIT2025",
                "Flag": 1,
                "Nf_exp": 1000.0,
                "Nf_expected": 900.0,
                "Nf_native": np.inf,
                "Error_ln": 0.1,
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 1.2,
            },
            {
                "Material": "Demo",
                "Method": "TVEIT2025",
                "Flag": 2,
                "Nf_exp": 1000000.0,
                "Nf_expected": 90000.0,
                "Nf_native": np.inf,
                "Error_ln": 0.1,
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 0.8,
            },
            {
                "Material": "Demo",
                "Method": "FIN",
                "Flag": 1,
                "Nf_exp": 1000.0,
                "Nf_expected": 900.0,
                "Nf_native": np.nan,
                "Error_ln": 0.1,
                "Implementation_revision": IMPLEMENTATION_REVISION,
                "CP_value": 100.0,
            },
        ]
    )


def test_native_life_omitted_from_generic_ranking_and_infinite_life_counted(tmp_path):
    result = compare.exclude_native_life_pending(_rows(), tmp_path)
    assert result["Method"].tolist() == ["FIN"]
    audit = pd.read_csv(tmp_path / "native_life_comparison_pending.csv")
    assert audit.loc[0, "Rows_excluded"] == 2
    assert audit.loc[0, "Failed_specimens_predicted_infinite_life"] == 1
    assert audit.loc[0, "Runouts_predicted_infinite_life"] == 1
    assert "Nf_native" in audit.loc[0, "Reason"]


def test_no_pending_native_models_leaves_frame_unchanged(tmp_path):
    frame = _rows().query("Method == 'FIN'")
    assert compare.exclude_native_life_pending(frame, tmp_path).equals(frame)
    assert not (tmp_path / "native_life_comparison_pending.csv").exists()


def test_native_only_comparison_exits_cleanly_with_audit(tmp_path):
    results = tmp_path / "Results"
    results.mkdir()
    _rows().query("Method == 'TVEIT2025'").to_csv(
        results / "Demo_TVEIT2025.csv", index=False
    )
    compare.main(["--results_dir", str(results), "--out_dir", str(tmp_path / "out")])
    audit = pd.read_csv(tmp_path / "out" / "native_life_comparison_pending.csv")
    assert audit["Rows_excluded"].sum() == 2
    assert not (tmp_path / "out" / "summary_error_stats.csv").exists()


def test_native_method_cannot_join_common_subset(tmp_path):
    results = tmp_path / "Results"
    results.mkdir()
    _rows().to_csv(results / "Demo.csv", index=False)
    with pytest.raises(ValueError, match="native-life comparison"):
        compare.main(
            [
                "--results_dir",
                str(results),
                "--out_dir",
                str(tmp_path / "out"),
                "--common-methods",
                "TVEIT2025",
                "FIN",
            ]
        )
