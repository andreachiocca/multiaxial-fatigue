"""Declarative corrections applied to as-supplied experimental workbooks.

An overlay file ``imp/overlay_files/<stem>.overlay.json`` (shipped with the package) lists the
cell-level changes that turn a workbook as supplied by its source into the revised workbook used by
this project. Overlays are found by content: a workbook gets the overlay whose ``base_sha256`` equals
its SHA-256, wherever the file lives and whatever it is called, so a copy keeps its corrections and an
edited file gets none. Every operation checks the value it replaces before writing. The importer
applies the overlay to the loaded worksheet, so everything downstream sees the same data as it does
from the revised workbook.
"""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any
import warnings

if TYPE_CHECKING:
    from multifat.imp.schema import ValidationReport

OVERLAY_DIR = Path(__file__).with_name("overlay_files")
SCHEMA_VERSION = 1


class OverlayError(ValueError):
    """An overlay operation does not match the workbook it is applied to."""


def sha256_file(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def _overlays_by_base_hash() -> dict[str, Path]:
    index: dict[str, Path] = {}
    for path in sorted(OVERLAY_DIR.glob("*.overlay.json")):
        base = json.loads(path.read_text(encoding="utf-8"))["base_sha256"]
        if base in index:
            raise OverlayError(f"{path.name} and {index[base].name} name the same base file")
        index[base] = path
    return index


@lru_cache(maxsize=None)
def _overlay_target_workbooks(overlay_dir: Path) -> set[str]:
    names: set[str] = set()
    for path in sorted(overlay_dir.glob("*.overlay.json")):
        wb_name = json.loads(path.read_text(encoding="utf-8")).get("workbook")
        if wb_name:
            names.add(str(wb_name))
    return names


def overlay_path(workbook: str | Path) -> Path | None:
    """The overlay written against this exact file content, or None."""
    return _overlays_by_base_hash().get(sha256_file(workbook))


def load_overlay(path: str | Path) -> dict[str, Any]:
    overlay = json.loads(Path(path).read_text(encoding="utf-8"))
    if overlay.get("schema_version") != SCHEMA_VERSION:
        raise OverlayError(f"{path}: unsupported schema_version {overlay.get('schema_version')!r}")
    return overlay


def _same(current: Any, expected: Any) -> bool:
    if isinstance(current, bool) or isinstance(expected, bool):
        return current is expected
    if isinstance(current, (int, float)) and isinstance(expected, (int, float)):
        return float(current) == float(expected)
    return current == expected


def _column_index(ws, header: str) -> int:
    for cell in ws[1]:
        if cell.value == header:
            return cell.column
    raise OverlayError(f"{ws.title}: no column headed {header!r}")


def _dataset_row(ws, authors: str, specimen: str) -> int:
    rows = [
        r
        for r in range(2, ws.max_row + 1)
        if ws.cell(r, 1).value == authors and ws.cell(r, 4).value == specimen
    ]
    if len(rows) != 1:
        raise OverlayError(f"Dataset: {len(rows)} rows match ({authors!r}, {specimen!r}), need 1")
    return rows[0]


def _apply_one(wb, op: dict[str, Any]) -> None:
    kind = op["op"]
    if kind == "set_summary":
        ws = wb["Summary"]
        cell = ws[f"{op['column']}{op['row']}"]
        if not _same(cell.value, op["from"]):
            raise OverlayError(
                f"Summary!{op['column']}{op['row']} is {cell.value!r}, expected {op['from']!r}"
            )
        cell.value = op["to"]
    elif kind == "set_case":
        ws = wb["Dataset"]
        row = _dataset_row(ws, op["authors"], op["specimen"])
        cell = ws.cell(row, _column_index(ws, op["field"]))
        if not _same(cell.value, op["from"]):
            raise OverlayError(
                f"Dataset row {row} ({op['specimen']}) {op['field']!r} is {cell.value!r}, "
                f"expected {op['from']!r}"
            )
        cell.value = op["to"]
    elif kind == "add_case":
        ws = wb["Dataset"]
        row = op["row"]
        if any(ws.cell(row, c).value is not None for c in range(1, ws.max_column + 1)):
            raise OverlayError(f"Dataset row {row} is not empty; add_case needs a free row")
        for header, value in op["values"].items():
            ws.cell(row, _column_index(ws, header)).value = value
    elif kind == "drop_case":
        # Clear one test row; row numbers of the other rows are kept, and the importer reads an
        # empty row as no test. The row is identified by number (specimen names need not be
        # unique) and checked against its authors, specimen name and the 'check' values.
        ws = wb["Dataset"]
        row = op["row"]
        expected = {
            "Authors": op["authors"],
            "Specimen name": op["specimen"],
            **op.get("check", {}),
        }
        for header, value in expected.items():
            current = ws.cell(row, _column_index(ws, header)).value
            if not _same(current, value):
                raise OverlayError(
                    f"Dataset row {row} {header!r} is {current!r}, expected {value!r}"
                )
        for c in range(1, ws.max_column + 1):
            ws.cell(row, c).value = None
    else:
        raise OverlayError(f"unknown overlay operation {kind!r}")


def apply_overlay_to_workbook(
    wb, workbook: str | Path, report: ValidationReport | None = None
) -> dict[str, Any] | None:
    """Apply the overlay of ``workbook`` to the loaded ``wb`` in memory.

    Returns the overlay identity recorded in provenance, or None when there is no overlay for
    this file or the file is not the base it was written against. Emits a ``RuntimeWarning`` and
    records a warning on ``report`` when an overlay exists for the file name but the SHA-256
    does not match.
    """
    path = overlay_path(workbook)
    if path is None:
        wb_name = Path(workbook).name
        if wb_name in _overlay_target_workbooks(OVERLAY_DIR):
            msg = f"{wb_name}: SHA-256 does not match any overlay; overlay was not applied"
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            if report is not None:
                report.warnings.append(msg)
        return None
    overlay = load_overlay(path)
    for op in overlay["operations"]:
        _apply_one(wb, op)
    return {
        "name": path.name,
        "sha256": sha256_file(path),
        "base_sha256": overlay["base_sha256"],
        "operations": len(overlay["operations"]),
        "expected_fingerprint": overlay["expected_fingerprint"],
    }
