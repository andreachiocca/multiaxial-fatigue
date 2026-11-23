"""Provenance manifests and staleness checks for committed homepage data (issue #80)."""

from __future__ import annotations

import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import DATA_DIR  # noqa: E402

RANKING_CSV = DATA_DIR / "load_error_ranking.csv"

RANKING_MATERIALS: tuple[str, ...] = (
    "42CrMo4_FABEST",
    "42CrMo4_FABEST_bar",
    "Al7075_T651_Wang",
    "42CrMo4_Santus",
    "42CrMo4_Reis",
    "AISI316L",
    "Al7075_T6",
    "Ti6Al4V_Meneghetti",
    "Ti6Al4V_Berto",
    # fatigue-limit datasets (d110; CX18 and GGG60 left out until #24 resolves their data)
    "25CrMo4",
    "30NCD16",
    "34Cr4",
    "982FA",
    "S65A",
    "Ti6Al4V_Delahay",
)


def manifest_path_for(csv_path: Path) -> Path:
    """Path of the provenance manifest paired with ``csv_path``."""
    return csv_path.with_suffix(".manifest.json")


def _csv_sha256(csv_path: Path) -> str:
    """Canonical LF SHA-256 so Windows autocrlf and Linux checkouts agree."""
    return hashlib.sha256(csv_path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def _git_head() -> str | None:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True, text=True
    )
    return proc.stdout.strip() or None


def _implementation_revision() -> str:
    from multifat.models.references import IMPLEMENTATION_REVISION

    return IMPLEMENTATION_REVISION


def registered_models() -> list[str]:
    """Sorted base names of all enabled registered models."""
    from multifat.models.references import DISABLED_METHODS
    from multifat.models.registry import param_table

    return sorted(name for name, _, _ in param_table() if name not in DISABLED_METHODS)


def dataset_fingerprints(materials: Sequence[str]) -> dict[str, str]:
    """Post-overlay dataset fingerprints keyed by material name.

    A dataset without failures (fatigue-limit states only) is hashed with its run-outs, which are
    its assessed observations; the fingerprint of a life dataset is the failure fingerprint.
    """
    from multifat.cli import run as runmod
    from multifat.report.subsets import dataset_fingerprint

    out: dict[str, str] = {}
    for mat in materials:
        data, resolved = runmod.load_material(mat)
        block = data[resolved]
        out[mat] = dataset_fingerprint(block, runouts=not block["cases"])
    return out


def write_manifest(
    csv_path: Path,
    materials: Sequence[str],
    methods_with_results: Sequence[str],
    *,
    manifest_path: Path | None = None,
) -> Path:
    """Write a provenance manifest alongside a generated homepage data CSV."""
    target = manifest_path or manifest_path_for(csv_path)
    ordered_materials = sorted(materials)
    payload = {
        "csv": csv_path.name,
        "csv_sha256": _csv_sha256(csv_path),
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git_head(),
        "implementation_revision": _implementation_revision(),
        "registered_models": registered_models(),
        "materials": ordered_materials,
        "dataset_fingerprints": dataset_fingerprints(ordered_materials),
        "methods_with_results": sorted(set(str(m) for m in methods_with_results)),
    }
    target.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return target


def check_manifest(
    csv_path: Path,
    expected_materials: Sequence[str],
    *,
    manifest_path: Path | None = None,
    current_fingerprints: dict[str, str] | None = None,
    current_models: Sequence[str] | None = None,
    current_revision: str | None = None,
) -> list[str]:
    """Reasons a homepage data CSV and its manifest are missing or stale."""
    label = csv_path.name
    mpath = manifest_path or manifest_path_for(csv_path)
    if not csv_path.is_file():
        return [f"{label}: missing CSV ({csv_path})"]
    if not mpath.is_file():
        return [f"{label}: missing manifest ({mpath})"]

    manifest = json.loads(mpath.read_text(encoding="utf-8"))
    problems: list[str] = []

    actual_sha = _csv_sha256(csv_path)
    if manifest.get("csv_sha256") != actual_sha:
        problems.append(f"{label}: CSV sha256 differs from manifest")

    expected = sorted(expected_materials)
    if sorted(manifest.get("materials", [])) != expected:
        problems.append(f"{label}: materials changed")

    models = sorted(registered_models() if current_models is None else current_models)
    recorded_models = list(manifest.get("registered_models", []))
    added = sorted(set(models) - set(recorded_models))
    removed = sorted(set(recorded_models) - set(models))
    if added or removed:
        problems.append(f"{label}: models added {added}, removed {removed}")

    revision = _implementation_revision() if current_revision is None else current_revision
    if manifest.get("implementation_revision") != revision:
        problems.append(f"{label}: implementation revision changed")

    fps = dataset_fingerprints(expected) if current_fingerprints is None else current_fingerprints
    recorded_fps = manifest.get("dataset_fingerprints") or {}
    changed = [m for m in expected if recorded_fps.get(m) != fps.get(m)]
    if changed:
        problems.append(f"{label}: datasets changed {changed}")
    return problems


def check() -> list[str]:
    """Reasons committed homepage data CSVs are out of date (empty list: current)."""
    fps = dataset_fingerprints(RANKING_MATERIALS)
    models = registered_models()
    revision = _implementation_revision()
    return [
        *check_manifest(
            RANKING_CSV,
            RANKING_MATERIALS,
            current_fingerprints=fps,
            current_models=models,
            current_revision=revision,
        ),
    ]
