"""Calibration history (issue #15): persistent ids, no overwriting, selection and compatibility."""

from __future__ import annotations

import json

import pytest

from multifat.cli.pipeline import main as pipeline_main
from multifat.cli.pipeline import run_workbook
from multifat.cli.records import (
    archive_record,
    check_record_compatibility,
    content_signature,
    find_record,
)
from workbook_factory import SUMMARY_HEADERS, SUMMARY_VALUES, write_workbook


def _no_k_fi(path):
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if not h.startswith("k FI")]
    return write_workbook(
        path,
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )


def test_signature_ignores_time_git_and_paths_only():
    record = {
        "params": {"k_FI": 0.3},
        "provenance": {"created_utc": "a", "git": {"commit": "x"}, "workbook": "p/W.xlsx", "n": 1},
    }
    moved = json.loads(json.dumps(record))
    moved["provenance"].update(created_utc="b", git={"commit": "y"}, workbook="q/W.xlsx")
    assert content_signature(record) == content_signature(moved)
    changed = json.loads(json.dumps(record))
    changed["params"]["k_FI"] = 0.31
    assert content_signature(record) != content_signature(changed)


def test_archive_never_overwrites_and_reuses_equivalent_records(tmp_path):
    first = archive_record({"params": {"k_FI": 0.3}}, tmp_path, "M")
    again = archive_record({"params": {"k_FI": 0.3}}, tmp_path, "M")
    other = archive_record({"params": {"k_FI": 0.4}}, tmp_path, "M")
    assert again == (first[0], first[1], True) and first[2] is False
    assert other[0] != first[0] and len(list((tmp_path / "M").glob("*.json"))) == 2
    stored = json.loads(first[1].read_text(encoding="utf-8"))
    assert stored["calibration_id"] == first[0]
    assert find_record(tmp_path, "M", first[0]) == first[1]
    with pytest.raises(FileNotFoundError, match="available"):
        find_record(tmp_path, "M", "nope")


def test_pipeline_archives_and_results_refer_to_the_record(tmp_path):
    wb = _no_k_fi(tmp_path / "H.xlsx")
    res = tmp_path / "res"
    run_workbook(wb, models=["FIN"], results_dir=res)
    side = json.loads((res / "H_calibration.json").read_text(encoding="utf-8"))
    calibration_id = side["calibration_id"]
    archived = res / "calibrations" / "H" / f"{calibration_id}.json"
    assert archived.is_file() and side["calibration_procedure"] == 2
    report = (res / "H_run_report.tsv").read_text(encoding="utf-8").splitlines()
    assert report[0].split("\t")[-1] == "calibration_id"
    assert report[1].split("\t")[-1] == calibration_id
    # An equivalent rerun keeps the same record (no duplicate)
    run_workbook(wb, models=["FIN"], results_dir=res)
    assert len(list((res / "calibrations" / "H").glob("*.json"))) == 1


def test_selection_by_id_reuses_constants(tmp_path):
    wb = _no_k_fi(tmp_path / "S.xlsx")
    first = tmp_path / "first"
    run_workbook(wb, models=["FIN"], results_dir=first)
    side = json.loads((first / "S_calibration.json").read_text(encoding="utf-8"))
    out = tmp_path / "second"
    with pytest.raises(SystemExit) as done:
        pipeline_main(
            [
                str(wb), "--models", "FIN", "--results-dir", str(out), "--no-compare",
                "--calibration-archive", str(first / "calibrations"),
                "--calibration-id", side["calibration_id"],
            ]
        )  # fmt: skip
    assert done.value.code in (0, None)
    reused = json.loads((out / "S_calibration.json").read_text(encoding="utf-8"))
    assert reused["params"]["k_FI"] == side["params"]["k_FI"]
    assert reused["constants_source"]["calibration_id"] == side["calibration_id"]
    assert "warning" not in reused["constants_source"]


def test_incompatible_records_are_rejected(tmp_path):
    record = {"params": {"k_FI": 0.3}, "provenance": {"failure_fingerprint": "x"}}
    assert "insufficiently documented" in check_record_compatibility(record, {})[0]
    for old_revision in ("1999-01-01-old", "2026-10-07-audit1"):
        record["provenance"]["implementation_revision"] = old_revision
        assert "other formulations" in check_record_compatibility(record, {})[0]
    wb = _no_k_fi(tmp_path / "Old.xlsx")
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="other formulations"):
        run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res", constants=stale)
