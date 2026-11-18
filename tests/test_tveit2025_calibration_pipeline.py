"""TVEIT2025 calibration in the pipeline: automatic, archived, reused, skipped with the reason.

A synthetic workbook (Summary Sigm1 = 420 MPa) carries axial R = -1 failures on the closed-form curve of
known constants (exact fit), one axial test with a mean stress and one torsion test. Its inputs file,
in a temporary Calibrations folder, sets A and coarse time-integration settings, and the C grid is
coarsened, so that a calibration takes seconds.
"""

import json

import pandas as pd
import pytest

from multifat.cli.pipeline import run_workbook
from multifat.eval.calibration_confidence import LEVELS
from multifat.eval import tveit2025_calibration as tc
from workbook_factory import write_workbook

TRUE = dict(C=0.5, K=2e-5, L=6.0, A=0.25, f_m1=420.0)


def _N(sigma_m, sigma_a):
    p = TRUE
    return tc.uniaxial_life(p["C"], p["K"], p["L"], sigma_m, sigma_a, p["A"], p["f_m1"]).Nf


def _row(N, sigma_m, sigma_a, tau_a=0.0):
    return [N, "No", "No", 0, 0, 0, 0, sigma_m, sigma_a, 0, tau_a]


@pytest.fixture
def synthetic(tmp_path, monkeypatch):
    """The workbook Syn.xlsx and its inputs in a temporary Calibrations folder."""
    monkeypatch.setattr(tc, "C_PER_DECADE", 1)
    monkeypatch.setattr(tc, "CALIBRATIONS_DIR", tmp_path / "Calibrations")
    rows = [_row(_N(0, s), 0, s) for s in (450.0, 480.0, 520.0, 600.0)]
    rows += [_row(_N(55.0, 450.0), 55.0, 450.0), _row(2e5, 0, 0, 280.0)]
    folder = tmp_path / "Calibrations" / "Syn"
    folder.mkdir(parents=True)
    inputs = {
        "schema": tc.INPUTS_SCHEMA,
        "material_id": "Syn",
        "A": {"value": TRUE["A"], "source": "the constants that generated the lives"},
        "numerical_settings": {"npts_TVEIT": 401, "max_blocks_TVEIT": 50},
    }
    (folder / tc.INPUTS_FILE).write_text(json.dumps(inputs), encoding="utf-8")
    return write_workbook(tmp_path / "Syn.xlsx", rows=rows)


def _run(workbook, results, **options):
    outcomes = run_workbook(workbook, models=["TVEIT2025"], results_dir=results, **options)
    sidecar = json.loads((results / "Syn_calibration.json").read_text(encoding="utf-8"))
    return outcomes, sidecar


def test_calibrated_without_supplied_constants_archived_and_flagged(tmp_path, synthetic):
    archive = tmp_path / "archive"
    (outcome,), sidecar = _run(synthetic, tmp_path / "res", calibration_archive=archive)
    entry = sidecar["model_calibrations"]["TVEIT2025"]
    record = entry["record"]
    assert entry["status"] == "calibrated" and record["mode"] == "author"
    assert [p["sigma_a"] for p in record["calibration_points"]] == [450.0, 520.0, 600.0]
    assert record["E_min"] < tc.EXACT_FIT and record["chosen"]["G1"] and record["chosen"]["G2"]
    assert record["inputs_sha256"] and record["mean_stress_slope"]["route"] == "supplied"
    params = sidecar["params"]
    assert params["A_TVEIT"] == TRUE["A"] and params["Sigm1_TVEIT"] == 420.0
    assert params["npts_TVEIT"] == 401 and params["max_blocks_TVEIT"] == 50
    stored = json.loads(
        (archive / "Syn" / f"{sidecar['calibration_id']}.json").read_text(encoding="utf-8")
    )
    assert stored["model_calibrations"]["TVEIT2025"]["params"]["C_TVEIT"] == params["C_TVEIT"]
    # the calibration tests are design rows (Flag = 0) of TVEIT2025
    table = pd.read_csv(tmp_path / "res" / "TVEIT2025_Syn.csv")
    design = table.loc[table.Flag == 0, "Nf_exp"].tolist()
    for point in record["calibration_points"]:
        assert any(abs(N / point["N_exp"] - 1.0) < 1e-9 for N in design)
    # the run report carries the confidence; a hand-entered A is a minor finding
    assert outcome.status == "ok" and "calibration confidence" in outcome.notes
    assert "A_supplied" in outcome.notes and entry["confidence"]["level"] in LEVELS


def test_a_second_run_reuses_the_archived_calibration(tmp_path, synthetic, monkeypatch):
    archive = tmp_path / "archive"
    _, first = _run(synthetic, tmp_path / "a", calibration_archive=archive)

    def fail(*args, **kwargs):
        raise AssertionError("calibrated again")

    monkeypatch.setattr(tc, "fit_ckl", fail)
    _, second = _run(synthetic, tmp_path / "b", calibration_archive=archive)
    entry = second["model_calibrations"]["TVEIT2025"]
    assert entry["reused_from"]["calibration_id"] == first["calibration_id"]
    assert second["params"] == first["params"]


def test_supplied_constants_are_not_calibrated(tmp_path, synthetic):
    constants = tmp_path / "c.json"
    hand = {k + "_TVEIT": TRUE[k] for k in ("A", "C", "K", "L")}
    record = {"hand_entered_constants": hand, "material_id": "Syn", "provenance": {"by": "test"}}
    constants.write_text(json.dumps(record), encoding="utf-8")
    (outcome,), sidecar = _run(synthetic, tmp_path / "res", constants=constants)
    assert "model_calibrations" not in sidecar
    assert outcome.status == "ok" and "calibration confidence" not in outcome.notes


def _no_fit(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("calibrated")

    monkeypatch.setattr(tc, "fit_ckl", fail)


@pytest.mark.parametrize(
    "by_hand",
    [{"Sigm1_TVEIT": 400.0}, {"npts_TVEIT": 801}, {"A_TVEIT": TRUE["A"]}],
    ids=["limit", "setting", "A"],
)
def test_partly_supplied_constants_are_not_calibrated(tmp_path, synthetic, monkeypatch, by_hand):
    # The calibration would fit C, K and L with other inputs than those of the evaluation
    _no_fit(monkeypatch)
    constants = tmp_path / "c.json"
    record = {"hand_entered_constants": by_hand, "material_id": "Syn", "provenance": {"by": "t"}}
    constants.write_text(json.dumps(record), encoding="utf-8")
    (outcome,), sidecar = _run(synthetic, tmp_path / "res", constants=constants)
    entry = sidecar["model_calibrations"]["TVEIT2025"]
    (key,) = by_hand
    assert entry["status"] == "partly_supplied" and key in entry["reason"]
    assert outcome.status == "skipped" and "entered by hand" in outcome.reason


def test_a_summary_constant_is_entered_by_hand(tmp_path, monkeypatch):
    # A Summary column named like a declared parameter is passed to the model, so it counts too
    monkeypatch.setattr(tc, "CALIBRATIONS_DIR", tmp_path / "Calibrations")
    _no_fit(monkeypatch)
    rows = [_row(_N(0, s), 0, s) for s in (450.0, 480.0, 520.0, 600.0)]
    workbook = write_workbook(tmp_path / "Syn.xlsx", rows=rows, extra_summary={"A_TVEIT": 0.3})
    (outcome,), sidecar = _run(workbook, tmp_path / "res")
    assert sidecar["model_calibrations"]["TVEIT2025"]["status"] == "partly_supplied"
    assert outcome.status == "skipped" and "A_TVEIT" in outcome.reason


def test_a_repeated_transfer_reuses_its_calibration(tmp_path, synthetic, monkeypatch):
    # A transfer archives its record under <target>_from_<design>; the next one must find it
    target = write_workbook(tmp_path / "Tgt.xlsx", rows=[_row(_N(0, 500.0), 0, 500.0)] * 2)
    archive = tmp_path / "archive"

    def transfer(results):
        run_workbook(
            target,
            models=["TVEIT2025"],
            results_dir=results,
            design_workbook=synthetic,
            calibration_archive=archive,
        )
        path = results / "Tgt_from_Syn_calibration.json"
        return json.loads(path.read_text(encoding="utf-8"))

    first = transfer(tmp_path / "a")
    assert first["model_calibrations"]["TVEIT2025"]["status"] == "calibrated"
    _no_fit(monkeypatch)
    second = transfer(tmp_path / "b")
    entry = second["model_calibrations"]["TVEIT2025"]
    assert entry["reused_from"]["calibration_id"] == first["calibration_id"]
    assert second["params"] == first["params"]


def test_not_calibrated_with_no_calibrate(tmp_path, synthetic):
    (outcome,), sidecar = _run(synthetic, tmp_path / "res", calibrate=False)
    assert "model_calibrations" not in sidecar
    assert outcome.status == "skipped" and "A_TVEIT" in outcome.reason


def test_without_three_levels_tveit_is_skipped_with_the_reason(tmp_path, synthetic):
    rows = [_row(_N(0, s), 0, s) for s in (450.0, 600.0)] + [_row(2e5, 0, 0, 280.0)]
    workbook = write_workbook(tmp_path / "Syn.xlsx", rows=rows)
    (outcome,), sidecar = _run(workbook, tmp_path / "res")
    entry = sidecar["model_calibrations"]["TVEIT2025"]
    assert entry["status"] == "insufficient_calibration_data"
    assert outcome.status == "skipped" and "three stress levels" in outcome.reason
