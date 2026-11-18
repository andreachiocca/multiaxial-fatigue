"""Observation identities, dataset fingerprints and common-subset comparison.

Ported from upstream tests/test_eligibility_subsets.py (A. Chiocca, 904c619) to pytest and our API.
"""

import csv

import numpy as np
import pandas as pd
import pytest

from multifat.imp.experimental_data import load_material_workbook
from multifat.report.export_csv import save_fatigue_csv
from multifat.report.subsets import (
    IDENTITY_COLUMNS,
    common_subset,
    dataset_fingerprint,
    observation_index,
    write_subset_reports,
)
from workbook_factory import basquin_rows, write_workbook


def _block(path, rows):
    (block,) = load_material_workbook(write_workbook(path, rows=rows)).values()
    return block


def _frame(rows):
    cols = [
        "Material",
        "Method",
        "Dataset_fingerprint",
        "Observation_id",
        "Case_id",
        "Test_status",
        "Nf_exp",
        "Nf_expected",
        "Error_ln",
        "Error_log10_dp",
        "Flag",
    ]
    return pd.DataFrame(rows, columns=cols)


def test_csv_appends_identity_columns(tmp_path):
    kw = dict(
        method_name="FS",
        material_name="X",
        N_exp=[1e5, 2e5],
        CP=[1.0, 0.9],
        N_pred=[1.1e5, 1.8e5],
        N_runout=[5e6],
        CP_runout=[0.5],
        N_pred_runout=[9e6],
        DP_fit_runout=[0.6],
    )
    ids = [
        {
            "Test_status": "failure",
            "Specimen_id": f"S{i}",
            "Case_id": "c",
            "Observation_id": f"o{i}",
            "Dataset_fingerprint": "f",
        }
        for i in range(2)
    ]
    run = [
        {
            "Test_status": "runout",
            "Specimen_id": "R",
            "Case_id": "r",
            "Observation_id": "or",
            "Dataset_fingerprint": "f",
        }
    ]
    a = save_fatigue_csv(out_dir=str(tmp_path / "a"), **kw)
    b = save_fatigue_csv(
        out_dir=str(tmp_path / "b"), observation_ids=ids, observation_ids_runout=run, **kw
    )
    ra, rb = (list(csv.reader(open(p, newline="", encoding="utf-8"))) for p in (a, b))
    assert rb[0][-len(IDENTITY_COLUMNS) :] == IDENTITY_COLUMNS
    assert [r[: len(ra[0])] for r in rb] == ra
    assert rb[-1][-5:] == ["runout", "R", "r", "or", "f"]


def test_fingerprint_ignores_runouts_and_row_order(tmp_path):
    rows = basquin_rows()
    base = dataset_fingerprint(_block(tmp_path / "A.xlsx", rows))
    swapped = rows[1:2] + rows[:1] + rows[2:]
    assert dataset_fingerprint(_block(tmp_path / "B.xlsx", swapped)) == base
    more = rows + [[7e6, "Yes", "No", 0, 0, 0, 0, 0, 250, 0, 0]]
    assert dataset_fingerprint(_block(tmp_path / "C.xlsx", more)) == base
    changed = [list(r) for r in rows]
    changed[0][0] *= 2
    assert dataset_fingerprint(_block(tmp_path / "D.xlsx", changed)) != base


def test_observation_index_separates_replicates(tmp_path):
    from workbook_factory import revised_row, write_revised_workbook

    rows = [revised_row({"Number of cycles to failure": n}) for n in (1e5, 3e5, 1e5)]
    rows.append(revised_row({"Run out (Yes/No)": "Yes", "Number of cycles to failure": 5e6}))
    (block,) = load_material_workbook(write_revised_workbook(tmp_path / "A.xlsx", rows)).values()
    idx = observation_index(block)
    first, twin = idx[("case", 0)], idx[("case", 2)]
    assert first["Case_id"] == twin["Case_id"] == idx[("case", 1)]["Case_id"]
    assert first["Observation_id"] != twin["Observation_id"]
    assert idx[("runout", 0)]["Test_status"] == "runout"


def test_runouts_excluded_from_common_subset():
    df = _frame(
        [
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
            ["X", "B", "f", "o1", "c1", "failure", 1e5, 2e5, -0.2, 0.02, 1],
            ["X", "A", "f", "o2", "c2", "runout", 5e6, 9e6, np.nan, np.nan, 2],
            ["X", "B", "f", "o2", "c2", "runout", 5e6, 8e6, np.nan, np.nan, 2],
        ]
    )
    sub = common_subset(df, ["A", "B"], "life")
    assert set(sub.Observation_id) == {"o1"} and len(sub) == 2


def test_common_subset_life_vs_damage():
    df = _frame(
        [
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
            ["X", "B", "f", "o1", "c1", "failure", 1e5, np.nan, np.nan, 0.02, 1],
        ]
    )
    assert common_subset(df, ["A", "B"], "life").empty
    assert len(common_subset(df, ["A", "B"], "damage")) == 2


def test_duplicate_method_observation_raises():
    df = _frame(
        [
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
            ["X", "B", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
        ]
    )
    with pytest.raises(ValueError, match="Duplicate"):
        common_subset(df, ["A", "B"], "life")


def test_subset_reports_written(tmp_path):
    df = _frame(
        [
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 1],
            ["X", "B", "f", "o1", "c1", "failure", 1e5, 2e5, -0.2, 0.02, 1],
        ]
    )
    write_subset_reports(df, tmp_path, ["A", "B"])
    assert (tmp_path / "individual_availability.csv").exists()
    assert (tmp_path / "common_subsets.json").exists()


def test_pipeline_csv_has_identities(tmp_path):
    from multifat.cli.pipeline import run_workbook

    path = write_workbook(tmp_path / "W.xlsx")
    run_workbook(path, models=["CROSSLAND"], results_dir=tmp_path / "res", calibrate=False)
    df = pd.read_csv(tmp_path / "res" / "CROSSLAND_W.csv")
    columns = list(df.columns)
    start = columns.index(IDENTITY_COLUMNS[0])  # later columns are appended after them
    assert columns[start : start + len(IDENTITY_COLUMNS)] == IDENTITY_COLUMNS
    assert columns[-3:] == ["Survival_probability", "Error_load", "SN_model"]
    assert df.Observation_id.is_unique and df.Dataset_fingerprint.nunique() == 1
    assert set(df.loc[df.Flag == 2, "Test_status"]) == {"runout"}
    assert set(df.loc[df.Flag != 2, "Test_status"]) == {"failure"}


def test_compare_common_methods(tmp_path):
    from multifat.cli import compare
    from multifat.cli.pipeline import run_workbook

    path = write_workbook(tmp_path / "W.xlsx")
    res = tmp_path / "res"
    run_workbook(path, models=["CROSSLAND", "SINES"], results_dir=res, calibrate=False)
    out = tmp_path / "cmp"
    compare.main(
        [
            "--results_dir",
            str(res),
            "--out_dir",
            str(out),
            "--flag_filter",
            "other",
            "--common-methods",
            "CROSSLAND",
            "SINES",
            "--comparison-metric",
            "damage",
        ]
    )
    sel = pd.read_csv(out / "selected_common_observations.csv")
    assert set(sel.Method) == {"CROSSLAND", "SINES"} and (sel.Flag == 1).all()
    assert sel.groupby("Method").size().nunique() == 1
    for name in ("subset_definition.txt", "individual_availability.csv", "common_subsets.json"):
        assert (out / name).exists()
    stats = pd.read_csv(out / "summary_error_stats.csv")
    assert set(stats["Method"]) == {"CROSSLAND", "SINES"}


def test_review_subsets_shim_help():
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    out = subprocess.run(
        [sys.executable, str(root / "scripts" / "3_review_subsets.py"), "--help"],
        capture_output=True,
        text=True,
        cwd=root,
    )
    assert out.returncode == 0 and "--common-methods" in out.stdout


def test_conflicting_flags_message_suggests_filter():
    """Minor: methods with different design sets under --flag_filter all get a clear hint."""
    df = _frame(
        [
            ["X", "A", "f", "o1", "c1", "failure", 1e5, 1e5, 0.1, 0.01, 0],
            ["X", "B", "f", "o1", "c1", "failure", 1e5, 2e5, -0.2, 0.02, 1],
        ]
    )
    with pytest.raises(ValueError, match="--flag_filter other"):
        common_subset(df, ["A", "B"], "life")


def test_design_identities_use_the_series_of_the_design_rows():
    """Minor: an _ext series filtered differently must not blank the design-row identities."""
    from types import SimpleNamespace

    from multifat.cli.run import design_identities

    obs = {("case", i): {"Test_status": "failure", "Observation_id": f"o{i}"} for i in range(4)}
    ev = SimpleNamespace(observations=obs, extra={"ids": [3]},
                         cal={"design_case_ids": [0, 1, 2], "design_case_ids_ext": [0, 1]})  # fmt: skip
    ids = design_identities(ev, "_ext", n_rows=4)
    assert [r["Observation_id"] for r in ids] == ["o0", "o1", "o2", "o3"]
