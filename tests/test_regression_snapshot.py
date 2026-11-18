"""tools/regression_snapshot.py: compare on the base snapshot's columns, data-folder option."""

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "regression_snapshot", ROOT / "tools/regression_snapshot.py"
)
snap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(snap)


def _snapshot(path: Path, csv_text: str) -> Path:
    (path / "runs").mkdir(parents=True)
    f = path / "runs" / "FS_X.csv"
    f.write_text(csv_text, encoding="utf-8")
    files = {"runs/FS_X.csv": hashlib.sha256(f.read_bytes()).hexdigest()}
    (path / "manifest.json").write_text(json.dumps({"codes": {}, "files": files}))
    return path


def test_compare_ignores_appended_columns(tmp_path):
    base = _snapshot(tmp_path / "a", "Material,Nf_exp\nX,100\nX,200\n")
    new = _snapshot(tmp_path / "b", "Material,Nf_exp,Test_status\nX,100,failure\nX,200,runout\n")
    assert snap.compare(base, new) == 0


def test_compare_flags_changed_base_column(tmp_path):
    base = _snapshot(tmp_path / "a", "Material,Nf_exp\nX,100\n")
    new = _snapshot(tmp_path / "b", "Material,Nf_exp,Test_status\nX,101,failure\n")
    assert snap.compare(base, new) == 1


def test_compare_flags_added_rows(tmp_path):
    base = _snapshot(tmp_path / "a", "Material,Nf_exp\nX,100\n")
    new = _snapshot(tmp_path / "b", "Material,Nf_exp\nX,100\nX,300\n")
    assert snap.compare(base, new) == 1


def test_material_args_from_data_folder(tmp_path):
    (tmp_path / "6_A.xlsx").write_bytes(b"")
    (tmp_path / "B.xlsx").write_bytes(b"")
    assert snap.material_args(tmp_path) == [str(tmp_path / "6_A.xlsx"), str(tmp_path / "B.xlsx")]


def test_calibrate_main_takes_data_folder(tmp_path):
    import pytest

    from multifat.cli import calibrate

    with pytest.raises(FileNotFoundError, match=tmp_path.name):
        calibrate.main(["--data-folder", str(tmp_path)])


def test_life_workbooks_come_from_the_classification():
    names = {p.stem for p in snap.life_workbooks(None)}
    assert {"6_42CrMo4_Santus", "14_42CrMo4_Reis", "16_42CrMo4_FABEST", "Al7075_T6"} <= names
    # fatigue-limit and notched-only workbooks go through the pipeline stage only
    assert not {"2_CX18", "4_30NCD16", "17_Ti6Al4V_Delahay", "10_En3B"} & names
    assert len(snap.workbooks(None)) > len(names)


def test_pipeline_hashes_ignore_volatile_fields(tmp_path):
    out = tmp_path / "snap"
    sidecar = out / "pipeline" / "6_X" / "X_calibration.json"
    sidecar.parent.mkdir(parents=True)
    report = sidecar.with_name("X_run_report.tsv")

    def write(stamp, seconds, k_fs):
        sidecar.write_text(json.dumps({"params": {"k_FS": k_fs}, "calibration_id": stamp,
                                       "provenance": {"created_utc": stamp, "git": {"commit": stamp}}}),
                           encoding="utf-8")  # fmt: skip
        report.write_text(f"model\tstatus\tseconds\tcalibration_id\nFS\tok\t{seconds}\t{stamp}\n",
                          encoding="utf-8")  # fmt: skip
        return snap.stable_bytes(out, sidecar), snap.stable_bytes(out, report)

    first = write("2026-11-12T10:00", 1.5, 1.0)
    assert write("2026-11-13T11:00", 9.9, 1.0) == first  # another run, same numbers
    assert write("2026-11-13T11:00", 9.9, 1.1)[0] != first[0]  # a changed k_FS is seen
    archived = out / "pipeline" / "6_X" / "calibrations" / "X" / "20261112_abc.json"
    archived.parent.mkdir(parents=True)
    archived.write_text("{}", encoding="utf-8")
    assert snap.stable_bytes(out, archived) is None
    script_csv = out / "runs" / "FS_X.csv"
    script_csv.parent.mkdir(parents=True)
    script_csv.write_bytes(b"a,b\n1,2\n")  # bytes: write_text would give \r\n on Windows
    assert snap.stable_bytes(out, script_csv) == b"a,b\n1,2\n"
