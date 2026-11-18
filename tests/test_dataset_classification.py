"""Workflow classification of imported datasets (issue #18): limit-only, life, mixed, ambiguous."""

from __future__ import annotations

import json

import pytest

from multifat.imp.classification import classify_dataset
from multifat.imp.experimental_data import load_material_workbook
from multifat.imp.schema import WorkbookError
from workbook_factory import SUMMARY_HEADERS, SUMMARY_VALUES, basquin_rows, write_workbook


def _classify(path):
    ((material, block),) = load_material_workbook(path).items()
    return classify_dataset(block, material)


def _limit_rows(cycles=1e6):
    rows = [list(r) for r in basquin_rows()]
    for r in rows:
        r[0], r[1] = cycles, "Yes"
    return rows


def _without_limits():
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if "fatigue limit" not in h]
    return dict(
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )


def test_life_dataset_with_runouts_keeps_them_censored(tmp_path):
    result = _classify(write_workbook(tmp_path / "Life.xlsx"))
    assert result["workflow"] == "fatigue_life"
    roles = result["roles"]
    assert len(roles["design"]) == 6 and len(roles["censored"]) == 1
    assert not roles["limit_state"]
    assert len(roles["design"]) + len(roles["assessment"]) == result["counts"]["failures"]
    assert any("censored life observations" in d for d in result["diagnostics"])
    ids = {r["Observation_id"] for group in roles.values() for r in group}
    assert len(ids) == sum(len(group) for group in roles.values())  # one role per observation


def test_limit_only_dataset(tmp_path):
    result = _classify(write_workbook(tmp_path / "Lim.xlsx", rows=_limit_rows()))
    assert result["workflow"] == "fatigue_limit" and result["fatigue_limits_usable"]
    assert len(result["roles"]["limit_state"]) == result["counts"]["runouts"]
    assert not result["roles"]["design"] and not result["roles"]["censored"]


def test_failures_without_a_design_set_are_ambiguous(tmp_path):
    rows = [r for r in basquin_rows() if r[8] == 0 or r[1] == "Yes"]  # no uniaxial tension failure
    result = _classify(write_workbook(tmp_path / "NoDesign.xlsx", rows=rows))
    assert result["workflow"] == "ambiguous" and "no master curve" in result["reasons"][0]
    assert not any(result["roles"].values())


def test_runouts_without_limits_are_ambiguous(tmp_path):
    path = write_workbook(tmp_path / "NoLim.xlsx", rows=_limit_rows(), **_without_limits())
    result = _classify(path)
    assert result["workflow"] == "ambiguous" and "Sigm1" in result["reasons"][0]


def test_notched_only_dataset_is_out_of_scope(tmp_path):
    rows = [r for r in basquin_rows() if r[2] == "Yes"]
    assert _classify(write_workbook(tmp_path / "Notch.xlsx", rows=rows))["workflow"] == (
        "out_of_scope"
    )


def test_pipeline_writes_the_classification_and_routes_on_it(tmp_path):
    from multifat.cli.pipeline import run_workbook

    res = tmp_path / "res"
    run_workbook(write_workbook(tmp_path / "W.xlsx"), models=["CROSSLAND"], results_dir=res)
    saved = json.loads((res / "W_classification.json").read_text(encoding="utf-8"))
    assert saved["workflow"] == "fatigue_life"
    run_workbook(
        write_workbook(tmp_path / "L.xlsx", rows=_limit_rows()), models=["CROSSLAND"],
        results_dir=res,
    )  # fmt: skip
    assert json.loads((res / "L_classification.json").read_text())["workflow"] == "fatigue_limit"
    assert (res / "FatigueLimit" / "CROSSLAND_L_FI.csv").exists()
    rows = [r for r in basquin_rows() if r[8] == 0 or r[1] == "Yes"]
    with pytest.raises(WorkbookError):
        run_workbook(write_workbook(tmp_path / "A.xlsx", rows=rows), results_dir=res)


def test_runouts_below_the_reference_life_are_reported():
    """25CrMo4: every run-out stopped below the 2e6-cycle reference of its Summary limits."""
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "Experimental_Data" / "3_25CrMo4.xlsx"
    result = _classify(path)
    assert result["workflow"] == "fatigue_limit"
    assert any("finite-cycle survivals" in d for d in result["diagnostics"])
