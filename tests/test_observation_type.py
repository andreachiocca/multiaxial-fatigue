"""Observation type column of the Dataset sheet (decision d113, issue #18).

The column declares a row a failure, an interrupted test (run-out, censored) or an established
fatigue limit (limit state). Declared types drive the classification; undeclared rows keep the rule
before d113; a mixed life workbook also gets the fatigue index of its declared limits.
"""

from __future__ import annotations

import sys
from pathlib import Path

import openpyxl
import pandas as pd
import pytest

from multifat.cli.pipeline import run_workbook
from multifat.imp.classification import classify_dataset, limit_state_runouts
from multifat.imp.experimental_data import load_material_workbook
from multifat.imp.schema import OBSERVATION_TYPES, WorkbookError, parse_observation_type
from workbook_factory import DATASET_HEADERS, basquin_rows, write_workbook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from add_observation_type import HEADER, add_column, observation_types  # noqa: E402

WORKBOOKS = sorted((ROOT / "Experimental_Data").glob("*.xlsx"))


def _synthetic(path: Path, runout_type: str | None, extra_limits=()) -> Path:
    """Synthetic life workbook with an Observation type column; its run-out typed `runout_type`."""
    rows = []
    for row in basquin_rows() + [list(r) for r in extra_limits]:
        kind = "failure" if row[1] == "No" else runout_type
        rows.append(list(row) + [kind])
    return write_workbook(path, dataset_headers=DATASET_HEADERS + [HEADER], rows=rows)


@pytest.mark.parametrize("path", WORKBOOKS, ids=lambda p: p.stem)
def test_every_supplied_row_is_declared_by_the_d113_rule(path):
    """Workbooks with failures hold run-outs; workbooks without failures hold fatigue limits."""
    ds = pd.read_excel(path, sheet_name="Dataset")
    run = next(c for c in ds.columns if "run" in str(c).lower() and "out" in str(c).lower())
    rows = ds[ds[run].notna()]
    kinds = rows[HEADER].map(parse_observation_type)
    assert kinds.notna().all()
    is_runout = rows[run].astype(str).str.lower().str.startswith("yes")
    assert ((kinds == "failure") == ~is_runout).all()
    expected = "run-out" if (~is_runout).any() else "fatigue limit"
    assert set(kinds[is_runout]) <= {expected}


@pytest.mark.parametrize("path", WORKBOOKS, ids=lambda p: p.stem)
def test_declared_types_reproduce_the_routing_before_d113(path, tmp_path):
    """Blanking the column gives the same workflow and roles as the declared types."""
    declared = load_material_workbook(path)
    for material, block in declared.items():
        undeclared = {
            key: {
                k: {**c, "meta": {**c["meta"], "observation_type": None}} for k, c in cases.items()
            }
            for key, cases in (("cases", block["cases"]), ("runouts", block["runouts"]))
        }
        legacy = classify_dataset({**block, **undeclared}, material)
        current = classify_dataset(block, material)
        assert (current["workflow"], current["roles"]) == (legacy["workflow"], legacy["roles"])


def test_parse_observation_type():
    assert parse_observation_type(" Fatigue  Limit ") == "fatigue limit"
    assert parse_observation_type("Run out") == parse_observation_type("runout") == "run-out"
    assert parse_observation_type("Null") is None and parse_observation_type(None) is None
    with pytest.raises(ValueError, match="not one of"):
        parse_observation_type("staircase")
    assert OBSERVATION_TYPES == ("failure", "run-out", "fatigue limit")


def test_a_type_that_contradicts_run_out_is_refused(tmp_path):
    wb = _synthetic(tmp_path / "Synth.xlsx", "run-out")
    book = openpyxl.load_workbook(wb)
    ds = book["Dataset"]
    ds.cell(2, len(DATASET_HEADERS) + 1).value = "fatigue limit"  # row 2 is a failure
    book.save(wb)
    with pytest.raises(WorkbookError, match="contradicts Run out No"):
        load_material_workbook(wb)


def test_declared_run_outs_without_failures_are_not_limit_states(tmp_path):
    rows = [[5e6, "Yes", "No", 0, 0, 0, 0, 0, 300, 0, 0, "run-out"]]
    wb = write_workbook(
        tmp_path / "Lim.xlsx", dataset_headers=DATASET_HEADERS + [HEADER], rows=rows
    )
    (block,) = load_material_workbook(wb).values()
    result = classify_dataset(block, "Lim")
    assert result["workflow"] == "ambiguous" and not result["roles"]["limit_state"]
    assert "declared as interrupted tests" in result["reasons"][0]


def test_mixed_workbook_gets_both_assessments(tmp_path):
    limits = [
        [1e7, "Yes", "No", 0, 0, 0, 0, 0, 0, 0, 290],
        [1e7, "Yes", "No", 0, 90, 0, 0, 0, 250, 0, 150],
    ]
    wb = _synthetic(tmp_path / "Mixed.xlsx", "fatigue limit", extra_limits=limits)
    (block,) = load_material_workbook(wb).values()
    result = classify_dataset(block, "Mixed")
    assert result["workflow"] == "fatigue_life" and "mixed" in result["reasons"][1]
    assert len(result["roles"]["limit_state"]) == 3 and not result["roles"]["censored"]
    assert len(limit_state_runouts(block, "Mixed")) == 3

    out = tmp_path / "res"
    outcomes = run_workbook(wb, models=["CROSSLAND"], results_dir=out)
    assert {o.model: o.status for o in outcomes} == {"CROSSLAND": "ok"}  # the life result
    assert (out / "CROSSLAND_Mixed.csv").exists()
    fi = pd.read_csv(out / "FatigueLimit" / "CROSSLAND_Mixed_FI.csv")
    assert len(fi) == 3 and fi["FI_load"].notna().all()
    assert (out / "FatigueLimit" / "Mixed_run_report.tsv").exists()


def test_censored_run_outs_of_a_life_workbook_get_no_fatigue_index(tmp_path):
    wb = _synthetic(tmp_path / "Life.xlsx", "run-out")
    out = tmp_path / "res"
    run_workbook(wb, models=["CROSSLAND"], results_dir=out)
    assert not (out / "FatigueLimit").exists()


def test_rule_of_the_migration_tool():
    assert observation_types({2: "No", 3: "Yes"}) == {2: "failure", 3: "run-out"}
    assert observation_types({2: "Yes", 3: "yes"}) == {2: "fatigue limit", 3: "fatigue limit"}
    with pytest.raises(ValueError, match="neither Yes nor No"):
        observation_types({2: "maybe"})


def test_migration_tool_keeps_values_and_refuses_a_second_pass(tmp_path):
    """On a workbook without the column the tool adds it and leaves every other value as it is."""
    plain = write_workbook(tmp_path / "Plain.xlsx")
    before = pd.read_excel(plain, sheet_name="Dataset")
    info = add_column(plain)
    after = pd.read_excel(plain, sheet_name="Dataset")
    assert after.drop(columns=[HEADER]).equals(before)
    assert (info["failure"], info["run-out"], info["fatigue limit"]) == (len(before) - 1, 1, 0)
    with pytest.raises(ValueError, match="already exists"):
        add_column(plain)
