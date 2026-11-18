"""Al7075_T651_Wang: fatigue-life dataset built from Wang, Shang, Sun & Chen (2018), Tables 1-5."""

from __future__ import annotations

import sys
from pathlib import Path

import openpyxl
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from build_wang_al7075_t651 import SOURCE, build, read_transcription  # noqa: E402

from multifat.imp.classification import classify_dataset  # noqa: E402
from multifat.imp.experimental_data import load_material_workbook  # noqa: E402

WORKBOOK = ROOT / "Experimental_Data" / "18_Al7075_T651_Wang.xlsx"


def _rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return [tuple(r) for r in wb[sheet].iter_rows(values_only=True)]


def test_committed_workbook_is_the_build(tmp_path):
    info = build(SOURCE, tmp_path / "wang.xlsx")
    for sheet in ("Summary", "Dataset"):
        assert _rows(tmp_path / "wang.xlsx", sheet) == _rows(WORKBOOK, sheet)
    assert (info["Sigm1"], info["Taum1"], info["nu"], info["rows"]) == (221.99, 127.1, 0.3036, 52)


def test_transcription_agrees_with_the_paper():
    """Table 4 prints tau_max = tau_a + tau_m; the table sizes are 13, 16 and 23 tests."""
    tests = read_transcription()
    assert [sum(t["table"] == k for t in tests) for k in "345"] == [13, 16, 23]
    torsion = {int(t["specimen"]): t["tau_a"] + t["tau_m"] for t in tests if t["table"] == "4"}
    assert torsion[14] == pytest.approx(225.17) and torsion[22] == pytest.approx(296.0)
    assert torsion[18] == pytest.approx(321.97) and torsion[23] == pytest.approx(285.62)


def test_imported_as_a_fatigue_life_dataset():
    (block,) = load_material_workbook(WORKBOOK).values()
    assert len(block["cases"]) == 52
    assert classify_dataset(block, "Al7075_T651_Wang")["workflow"] == "fatigue_life"


def test_fingerprint_ignores_last_bit_trigonometric_differences():
    """30 and 45 deg failures hold cos/sin values that differ by one ulp between platforms
    (Windows CI); their histories are hashed rounded to 1e-6 MPa, so the fingerprint agrees."""
    import copy

    import numpy as np

    from multifat.report.subsets import dataset_fingerprint

    (block,) = load_material_workbook(WORKBOOK).values()
    perturbed = copy.deepcopy(block)
    for case in perturbed["cases"].values():
        for key, value in case.get("harmonics", {}).items():
            if key.startswith("S_"):
                case["harmonics"][key] = np.nextafter(np.asarray(value, float), np.inf)
    assert dataset_fingerprint(perturbed) == dataset_fingerprint(block)
