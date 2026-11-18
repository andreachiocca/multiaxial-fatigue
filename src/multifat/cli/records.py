"""Calibration history of the pipeline: persistent, never-overwritten calibration records (#15).

Every calibration sidecar that the pipeline writes gets an identity::

    calibration_id = <UTC creation time>_<first 12 hex digits of its content signature>

and a copy is archived as ``<archive>/<MATERIAL>/<calibration_id>.json``. The content signature
is the SHA-256 of the record without the fields that do not change the calibration (creation
time, Git state, file paths), so an equivalent rerun finds its earlier record and reuses that
identity instead of creating a duplicate. Archived files are created exclusively and never
replaced. ``multifat-pipeline --calibration-id ID`` selects a record of the archive; its
compatibility is checked like any ``--constants`` file (``check_record_compatibility``).
"""

from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path

# Version of the calibration procedure of multifat.cli.calibrate (selection rules, search
# intervals, acceptance of optima). Increment it when these change. 2: design flags of the
# scatter-objective tests, k_MGSE from Zhu Eq. (12), L_LI in [0, 20], bound-limited constants
# unidentified (PR #34).
CALIBRATION_PROCEDURE = 2
VOLATILE_PROVENANCE = ("created_utc", "git", "workbook")


def content_signature(record: dict) -> str:
    """SHA-256 of a record without its volatile fields (time, Git state, paths, its own id)."""
    stable = {k: v for k, v in record.items() if k not in ("calibration_id", "constants_source")}
    if "provenance" in stable:
        stable["provenance"] = {
            k: v for k, v in stable["provenance"].items() if k not in VOLATILE_PROVENANCE
        }
    source = record.get("constants_source")
    if source:  # where supplied constants came from: by content, not by path
        stable["constants_source"] = {k: v for k, v in source.items() if k != "file"}
    blob = json.dumps(stable, sort_keys=True, default=float, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def archive_record(record: dict, archive: Path, material: str) -> tuple[str, Path, bool]:
    """Store `record` under a persistent id; (id, path, reused) where reused marks an equivalent
    earlier record whose id is returned instead of writing a duplicate."""
    folder = Path(archive) / material
    folder.mkdir(parents=True, exist_ok=True)
    digest = content_signature(record)[:12]
    earlier = sorted(folder.glob(f"*_{digest}.json"))
    if earlier:
        return earlier[0].stem, earlier[0], True
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    calibration_id = f"{stamp}_{digest}"
    path = folder / f"{calibration_id}.json"
    with open(path, "x", encoding="utf-8") as handle:  # exclusive: never replaces a record
        json.dump({**record, "calibration_id": calibration_id}, handle, indent=1, default=float)
    return calibration_id, path, False


def find_record(archive: Path, material: str, calibration_id: str) -> Path:
    """Path of the archived record `calibration_id` of `material`, or a clear error."""
    path = Path(archive) / material / f"{calibration_id}.json"
    if not path.is_file():
        available = sorted(p.stem for p in (Path(archive) / material).glob("*.json"))
        raise FileNotFoundError(
            f"no calibration {calibration_id!r} for {material} in {archive}; "
            f"available: {available or 'none'}"
        )
    return path


def check_record_compatibility(record: dict, data_block: dict) -> list[str]:
    """Reasons a pipeline sidecar cannot be reused here (empty list: compatible).

    Rejected: a record without the provenance needed to judge it, or made with another
    implementation revision of the methods (the constants belong to other formulations).
    Data and procedure differences are not rejections: they are warnings returned by
    ``record_warnings``, since transferring constants to new data can be intended.
    """
    from multifat.models.references import IMPLEMENTATION_REVISION

    prov = record.get("provenance") or {}
    missing = [k for k in ("implementation_revision", "failure_fingerprint") if k not in prov]
    if missing:
        return [f"insufficiently documented: provenance lacks {', '.join(missing)}"]
    if prov["implementation_revision"] != IMPLEMENTATION_REVISION:
        return [
            f"made with implementation revision {prov['implementation_revision']}, the methods "
            f"are now {IMPLEMENTATION_REVISION}: its constants belong to other formulations"
        ]
    return []


def record_warnings(record: dict, data_block: dict) -> list[str]:
    """Differences that do not prevent reuse but must be recorded with the results."""
    from multifat.report.subsets import dataset_fingerprint

    prov = record.get("provenance") or {}
    notes = []
    if prov.get("failure_fingerprint") != dataset_fingerprint(data_block):
        notes.append("made on other failure data: the constants are transferred, not calibrated")
    if record.get("calibration_procedure", 1) != CALIBRATION_PROCEDURE:
        notes.append(
            f"calibration procedure {record.get('calibration_procedure', 1)}, current "
            f"{CALIBRATION_PROCEDURE}: a recalibration may give different constants"
        )
    return notes
