"""42CrMo4_FABEST: the committed workbook is the conversion of the FABEST source file.

Source: FABEST database FAB_v20261014, doi:10.5281/zenodo.23359038 (CC BY 4.0), kept unchanged in
Experimental_Data/sources/FABEST/. tools/convert_fabest.py documents the scope and reading rules.
"""

from pathlib import Path

import openpyxl
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Experimental_Data/sources/FABEST/FAB_v20261014.xlsm"
WORKBOOK = ROOT / "Experimental_Data/16_42CrMo4_FABEST.xlsx"
pytestmark = pytest.mark.skipif(not SOURCE.exists(), reason="FABEST source not present")


def _rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return [tuple(r) for r in wb[sheet].iter_rows(values_only=True)]


def test_committed_workbook_is_the_conversion(tmp_path):
    import sys

    sys.path.insert(0, str(ROOT / "tools"))
    from convert_fabest import convert

    info = convert(SOURCE, tmp_path / "FABEST.xlsx")
    for sheet in ("Summary", "Dataset"):
        assert _rows(tmp_path / "FABEST.xlsx", sheet) == _rows(WORKBOOK, sheet)
    assert info["rows"] == 89
    assert sum(f for f, _ in info["counts"].values()) == 79
    assert sum(r for _, r in info["counts"].values()) == 10
    # Authors' Basquin strengths at 1e6 cycles (sheets FAB060/070 and FAB062/072, SN_FL)
    assert info["Sigm1"] == pytest.approx((514.92456 + 519.64165) / 2, abs=1e-4)
    assert info["Taum1"] == pytest.approx((368.11606 + 364.89692) / 2, abs=1e-4)


def test_imported_dataset():
    from multifat.imp.experimental_data import load_material_workbook

    block = load_material_workbook(WORKBOOK)["42CrMo4_FABEST"]
    assert (len(block["cases"]), len(block["runouts"])) == (79, 10)
    assert not block["excluded_observations"]
    phases = {c["meta"].get("phase_tor_deg") for c in block["cases"].values()}
    assert {0.0, 90.0} <= phases
    assert block["material"]["E_MPa"] == pytest.approx(206362.0)
