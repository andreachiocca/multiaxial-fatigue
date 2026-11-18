"""Every calibration record kept under Calibrations/ loads through the pipeline loader.

The records come from the upstream calibration history and are kept unchanged. Each must give
constants for the pipeline, and a Findley amplitude setting, when it has one, must be the production
MCE path (the only one the pipeline accepts).
"""

from pathlib import Path

import pytest

from multifat.cli.pipeline import load_constants

ROOT = Path(__file__).resolve().parents[1]
RECORDS = sorted(p for p in (ROOT / "Calibrations").glob("*/*.json"))


def test_records_are_present():
    assert RECORDS, "Calibrations/ holds no records"


@pytest.mark.parametrize("record", RECORDS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_record_loads_through_the_pipeline_loader(record):
    constants, origin = load_constants(record, record.parent.name)
    assert constants, f"{record.name}: no constants"
    assert all(isinstance(v, float) for v in constants.values())
    assert origin["kind"]
    if origin["kind"] == "upstream calibration record":
        assert "implementation revision" in origin.get("warning", "")


@pytest.mark.parametrize("record", RECORDS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_findley_amplitude_metadata_is_read(record):
    import json

    payload = json.loads(record.read_text(encoding="utf-8"))
    settings = (payload.get("amplitude_settings") or {}).get("FIN")
    if settings is None:
        # Legacy records predate the amplitude metadata; they are the MCE production path.
        return
    assert settings["algorithm"] == "mce_harmonic"
    # The loader accepts the record only for the production path
    load_constants(record, record.parent.name)


def test_tveit_stage1_file_runs_through_the_pipeline_constants(tmp_path, monkeypatch):
    """#60: a hand-entered TVEIT2025 file is accepted by --constants and labelled as such.

    The model run itself is stubbed (a TVEIT2025 integration takes minutes); the constants that
    reach it, Sigm1_TVEIT and Taum1_TVEIT included, and the sidecar label are checked.
    """
    import json

    from multifat.cli import pipeline

    seen = {}

    def fake_base(job):
        seen.update(job[3])
        return [pipeline.Outcome(job[1], "ok")]

    monkeypatch.setattr(pipeline, "_evaluate_base", fake_base)
    record = (
        Path(__file__).resolve().parents[1] / "Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json"
    )
    outcomes = pipeline.run_workbook(
        Path(__file__).resolve().parents[1] / "Experimental_Data/6_42CrMo4_Santus.xlsx",
        models=["TVEIT2025"],
        results_dir=tmp_path,
        constants=record,
        calibrate=False,
    )
    assert [o.status for o in outcomes] == ["ok"]
    assert (seen["Sigm1_TVEIT"], seen["Taum1_TVEIT"]) == (433.0, 344.0)
    assert {"A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT", "npts_TVEIT"} <= set(seen)
    sidecar = json.loads((tmp_path / "42CrMo4_Santus_calibration.json").read_text("utf-8"))
    source = sidecar["constants_source"]
    assert source["kind"] == "hand-entered constants"
    assert "entered_by" in source["provenance"]


def test_hand_entered_constants_must_name_their_material(tmp_path):
    import json

    from multifat.cli.pipeline import load_constants

    path = tmp_path / "c.json"
    path.write_text(json.dumps({"hand_entered_constants": {"A_TVEIT": 0.2},
                                "material_id": "Other"}), encoding="utf-8")  # fmt: skip
    with pytest.raises(ValueError, match="holds constants for material 'Other'"):
        load_constants(path, "42CrMo4_Santus")
