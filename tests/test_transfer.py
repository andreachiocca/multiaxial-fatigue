"""Cross-batch transfer (issue #71): predict the tests of one batch from another one.

The design side supplies the Summary properties, the constants, the design set and the master
curve; only the target's tests are assessed. Checked on the 42CrMo4 batches.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multifat.cli import pipeline
from multifat.cli import run as runmod

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Experimental_Data"
DESIGN, TARGET = "42CrMo4_Santus", "42CrMo4_Reis"
SSF_REIS = ROOT / "Calibrations" / "42CrMo4_Reis" / "SSF_published.json"


@pytest.fixture(scope="module")
def transfer():
    return runmod.load_transfer(TARGET, DESIGN)


def test_transfer_dataset_layout(transfer):
    data, material, design_data = transfer
    assert material == "42CrMo4_Reis_from_42CrMo4_Santus"
    d_data, d_mat = runmod.load_material(DESIGN)
    t_data, t_mat = runmod.load_material(TARGET)
    block = data[material]
    # resistance from the design batch
    assert block["material"] == d_data[d_mat]["material"]
    assert design_data[material]["cases"].keys() == d_data[d_mat]["cases"].keys()
    assert design_data[material]["runouts"] == {}
    # loads: every failure and run-out of the target, ids shifted past the design ids
    off = runmod.TRANSFER_ID_OFFSET
    target_ids = {i for i in block["cases"] if i >= off}
    assert target_ids == {int(i) + off for i in t_data[t_mat]["cases"]}
    assert set(block["runouts"]) == {int(i) + off for i in t_data[t_mat]["runouts"]}
    assert block["transfer"] == {"target": TARGET, "design": DESIGN}


def test_design_side_must_be_a_subset(transfer):
    data, material, design_data = transfer
    wrong = {material: {**design_data[material], "cases": {-1: {}}}}
    with pytest.raises(ValueError, match="subset"):
        runmod.evaluate_model(
            data, material, "CROSSLAND", _params(data, material), design_data=wrong
        )


def _params(data, material):
    return runmod.build_params(runmod.exp.get_material_properties(data, material))


def test_transfer_reuses_the_design_batch_master_curve(transfer):
    data, material, design_data = transfer
    own_data, own_mat = runmod.load_material(DESIGN)
    own = runmod.evaluate_model(own_data, own_mat, "CROSSLAND", _params(own_data, own_mat))
    moved = runmod.evaluate_model(
        data, material, "CROSSLAND", _params(data, material), design_data=design_data
    )
    assert moved.cal["design_case_ids"] == own.cal["design_case_ids"]
    for key in ("A", "b"):
        assert moved.cal["fit"][key] == own.cal["fit"][key]
    # only the target's failures are scored
    off = runmod.TRANSFER_ID_OFFSET
    t_data, t_mat = runmod.load_material(TARGET)
    assert sorted(moved.case_ids) == sorted(int(i) + off for i in t_data[t_mat]["cases"])
    assert moved.extra["ids"] == []


def test_default_evaluation_is_unchanged_by_the_option():
    data, material = runmod.load_material(TARGET)
    a = runmod.evaluate_model(data, material, "CROSSLAND", _params(data, material))
    b = runmod.evaluate_model(
        data, material, "CROSSLAND", _params(data, material), design_data=None
    )
    assert a.case_ids == b.case_ids and a.values == b.values


def test_transfer_run_writes_design_rows_of_the_design_batch(tmp_path):
    outcomes = pipeline.run_workbook(
        DATA / "14_42CrMo4_Reis.xlsx",
        models=["CROSSLAND"],
        results_dir=tmp_path,
        design_workbook=DATA / "6_42CrMo4_Santus.xlsx",
    )
    assert [o.status for o in outcomes] == ["ok"]
    frame = pd.read_csv(tmp_path / "CROSSLAND_42CrMo4_Reis_from_42CrMo4_Santus.csv")
    own = pd.read_csv(_own_run(tmp_path / "own", "6_42CrMo4_Santus.xlsx"))
    design, own_design = frame[frame["Flag"] == 0], own[own["Flag"] == 0]
    assert sorted(design["Observation_id"]) == sorted(own_design["Observation_id"])
    assert np.array_equal(
        design.sort_values("Observation_id")["DP_fit"].to_numpy(),
        own_design.sort_values("Observation_id")["DP_fit"].to_numpy(),
    )
    assessed = frame[frame["Flag"] == 1]
    assert assessed["Specimen_id"].str.startswith("Plain_multiaxial").all()
    assert len(assessed) == 25  # every Reis failure, the axial ones included


def _own_run(dest: Path, workbook: str) -> Path:
    pipeline.run_workbook(DATA / workbook, models=["CROSSLAND"], results_dir=dest)
    return dest / "CROSSLAND_42CrMo4_Santus.csv"


def test_transfer_refuses_constants_of_another_material(tmp_path):
    # a published-constants record of Reis cannot serve a transfer designed on QT
    with pytest.raises(ValueError, match="published-constants record of the design material"):
        pipeline.run_workbook(
            DATA / "14_42CrMo4_Reis.xlsx",
            models=["CROSSLAND"],
            results_dir=tmp_path,
            constants=SSF_REIS,
            design_workbook=DATA / "6_42CrMo4_Santus.xlsx",
        )


def test_transfer_refuses_calibration_records(tmp_path):
    sidecar = tmp_path / "sidecar.json"
    sidecar.write_text(json.dumps({"params": {"k_FS": 1.0}}), encoding="utf-8")
    assert not pipeline.is_published_constants(sidecar, DESIGN)
    with pytest.raises(ValueError, match="published-constants record"):
        pipeline.run_workbook(
            DATA / "14_42CrMo4_Reis.xlsx",
            models=["CROSSLAND"],
            results_dir=tmp_path,
            constants=sidecar,
            design_workbook=DATA / "6_42CrMo4_Santus.xlsx",
        )


def test_transfer_takes_published_constants_of_the_design_batch(tmp_path):
    """SSF: the surface published for the Reis batch, assessed on FABEST (out of sample)."""
    outcomes = pipeline.run_workbook(
        DATA / "16_42CrMo4_FABEST.xlsx",
        models=["SSF"],
        results_dir=tmp_path,
        constants=SSF_REIS,
        design_workbook=DATA / "14_42CrMo4_Reis.xlsx",
    )
    assert [o.status for o in outcomes] == ["ok"]
    sidecar = json.loads(
        (tmp_path / "42CrMo4_FABEST_from_42CrMo4_Reis_calibration.json").read_text("utf-8")
    )
    assert sidecar["constants_source"]["kind"] == "published constants"
    assert sidecar["params"]["a_SSF"] == 2.69
    frame = pd.read_csv(tmp_path / "SSF_42CrMo4_FABEST_from_42CrMo4_Reis.csv")
    assert frame[frame["Flag"] == 1]["Specimen_id"].str.startswith("FAB").all()
    assert frame[frame["Flag"] == 0]["Specimen_id"].str.startswith("Plain_multiaxial").all()


def test_transfer_native_life_model(tmp_path):
    """#80: MOREL (native life law) in a transfer keeps the design batch's own native lives."""
    outcomes = pipeline.run_workbook(
        DATA / "14_42CrMo4_Reis.xlsx",
        models=["MOREL"],
        results_dir=tmp_path,
        design_workbook=DATA / "6_42CrMo4_Santus.xlsx",
    )
    assert [o.status for o in outcomes] == ["ok"]
    moved = pd.read_csv(tmp_path / "MOREL_42CrMo4_Reis_from_42CrMo4_Santus.csv")
    pipeline.run_workbook(
        DATA / "6_42CrMo4_Santus.xlsx", models=["MOREL"], results_dir=tmp_path / "own"
    )
    own = pd.read_csv(tmp_path / "own" / "MOREL_42CrMo4_Santus.csv")
    design = moved[moved["Flag"] == 0].sort_values("Observation_id")
    own_design = own[own["Flag"] == 0].sort_values("Observation_id")
    assert list(design["Observation_id"]) == list(own_design["Observation_id"])
    assert np.array_equal(design["Nf_native"].to_numpy(), own_design["Nf_native"].to_numpy())
    assessed = moved[moved["Flag"] == 1]
    assert len(assessed) == 25 and np.isfinite(assessed["Error_load"]).all()
