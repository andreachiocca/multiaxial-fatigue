"""Ti6Al4V_Delahay: fatigue-limit dataset built from Delahay & Palin-Luc 2006 (decision d109).

Source: Int J Fatigue 28 (2006) 474-484, doi:10.1016/j.ijfatigue.2005.06.048, Tables 1 and 2,
transcribed in Experimental_Data/sources/Delahay2006/ti6al4v_tables.csv; tools/
build_delahay_ti6al4v.py documents the scope and builds the workbook.
"""

import math
import sys
from pathlib import Path

import openpyxl
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from build_delahay_ti6al4v import SOURCE, build, read_transcription, threshold_stress  # noqa: E402

from multifat.imp.classification import classify_dataset  # noqa: E402
from multifat.imp.experimental_data import load_material_workbook  # noqa: E402

WORKBOOK = ROOT / "Experimental_Data" / "17_Ti6Al4V_Delahay.xlsx"


def _rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return [tuple(r) for r in wb[sheet].iter_rows(values_only=True)]


def test_committed_workbook_is_the_build(tmp_path):
    info = build(SOURCE, tmp_path / "delahay.xlsx")
    for sheet in ("Summary", "Dataset"):
        assert _rows(tmp_path / "delahay.xlsx", sheet) == _rows(WORKBOOK, sheet)
    assert (info["Sigm1"], info["Taum1"], info["rows"]) == (652.0, 411.0, 4)


def test_transcription_agrees_with_the_paper_equations():
    """Checks of the transcribed limits against quantities the paper derives from them."""
    props, limits = read_transcription()
    # Eq. (6): sigma* = sqrt(2 sigma_Ten^2 - sigma_RotBen^2) = 563 MPa (Table 1)
    assert threshold_stress(limits) == pytest.approx(props["sigma_star_MPa"], abs=0.5)
    # Section 4.2, p. 483: tau*/tau_ref = 0.934 with tau_ref = 411 and tau* = 384 MPa
    torsion = next(lim for lim in limits if lim["load"] == "torsion")
    assert torsion["ta"] == 411.0 and 384.0 / torsion["ta"] == pytest.approx(0.934, abs=5e-4)
    # Fig. 2: the combined loadings have sigma_a / tau_a = 1.732
    for lim in limits:
        if lim["load"] == "plane bending + torsion":
            assert lim["sa"] / lim["ta"] == pytest.approx(math.sqrt(3.0), abs=5e-3)


def test_imported_as_a_fatigue_limit_dataset():
    block = load_material_workbook(WORKBOOK)["Ti6Al4V_Delahay"]
    m = block["material"]
    assert (m["E_MPa"], m["nu"], m["Tensile_strength_MPa"], m["Yield_strength_MPa"]) == (
        110000.0,
        0.3,
        1090.0,
        980.0,
    )
    assert (m["Tensile_fatigue_limit_MPa"], m["Torsional_fatigue_limit_MPa"]) == (652.0, 411.0)
    (sigm1,) = block["fatigue_limit_sources"]["Sigm1"]
    assert sigm1["loading"] == "plane bending" and sigm1["reference_cycles"] == 1e6
    assert not block["identity_issues"]
    assert (len(block["cases"]), len(block["runouts"])) == (0, 4)
    phases = sorted(r["meta"].get("phase_tor_deg") or 0.0 for r in block["runouts"].values())
    assert phases == [0.0, 0.0, 0.0, 90.0]
    assert classify_dataset(block, "Ti6Al4V_Delahay")["workflow"] == "fatigue_limit"
