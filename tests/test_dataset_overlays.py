"""Overlay layer: as-supplied workbooks plus declared corrections reproduce the revised dataset.

The workbooks edited after supply are kept byte-identical to their source; their corrections
are shipped with the package. Each overlay is found by the SHA-256 of the file it was written against,
so it follows a copy of that file and never applies to an edited one.
"""

import hashlib
import json
from pathlib import Path
import warnings

import pytest
from openpyxl import load_workbook

from multifat.imp.experimental_data import (
    import_experimental_data_grouped_by_material,
    load_material_workbook,
)
from multifat.imp import overlays as overlays_mod
from multifat.imp.overlays import (
    OverlayError,
    _apply_one,
    apply_overlay_to_workbook,
    load_overlay,
    overlay_path,
)
from multifat.imp.schema import ValidationReport, workbook_files
from multifat.report.subsets import dataset_fingerprint

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "Experimental_Data"
OVERLAID = [
    "Al7075_T6",
    "2_CX18",
    "3_25CrMo4",
    "4_30NCD16",
    "7_982FA",
    "8_AISI316L",
    "14_42CrMo4_Reis",
    "1_Ti6Al4V_Meneghetti",
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("stem", OVERLAID)
def test_overlay_names_the_as_supplied_bytes(stem):
    path = overlay_path(FOLDER / f"{stem}.xlsx")
    assert path is not None
    overlay = json.loads(path.read_text(encoding="utf-8"))
    assert overlay["base_sha256"] == _sha256(FOLDER / f"{stem}.xlsx")
    assert overlay["operations"], "an overlay with no operations should not exist"


@pytest.mark.parametrize("stem", OVERLAID)
def test_overlay_reproduces_expected_fingerprint(stem):
    overlay = json.loads(overlay_path(FOLDER / f"{stem}.xlsx").read_text(encoding="utf-8"))
    data = load_material_workbook(FOLDER / f"{stem}.xlsx")
    (block,) = data.values()
    assert block["overlay"]["sha256"] is not None
    assert block["overlay"]["operations"] == len(overlay["operations"])
    assert dataset_fingerprint(block) == overlay["expected_fingerprint"]


def test_aisi316l_counts_match_provenance():
    block = load_material_workbook(FOLDER / "8_AISI316L.xlsx")["AISI316L"]
    assert len(block["cases"]) == 39
    assert len(block["runouts"]) == 13


def test_overlay_off_gives_the_as_supplied_dataset():
    on = load_material_workbook(FOLDER / "8_AISI316L.xlsx")["AISI316L"]
    off = load_material_workbook(FOLDER / "8_AISI316L.xlsx", overlay=False)["AISI316L"]
    assert off["overlay"] is None
    assert len(off["cases"]) == 31
    assert on["material"]["E_MPa"] == 193000.0
    assert off["material"]["E_MPa"] != on["material"]["E_MPa"]


def test_no_overlay_for_workbooks_without_one():
    block = load_material_workbook(FOLDER / "15_Ti6Al4V_Berto.xlsx")["Ti6Al4V_Berto"]
    assert block["overlay"] is None


def test_folder_import_applies_overlays_by_default():
    data = import_experimental_data_grouped_by_material(FOLDER)
    assert data["CX18"]["material"]["nu"] == 0.3
    assert data["AISI316L"]["overlay"]["name"] == "8_AISI316L.overlay.json"


def test_tampered_from_value_aborts(tmp_path):
    wb = load_workbook(FOLDER / "8_AISI316L.xlsx", data_only=True)
    wb["Summary"]["E2"].value = 1.0
    with pytest.raises(OverlayError, match="expected"):
        apply_overlay_to_workbook(wb, FOLDER / "8_AISI316L.xlsx")


@pytest.mark.parametrize(
    ("op", "match"),
    [
        (
            {
                "op": "set_case",
                "authors": "Janssens",
                "specimen": "Uniaxial_1",
                "field": "Number of cycles to failure",
                "from": -1,
                "to": 100,
            },
            "expected",
        ),
        (
            {
                "op": "set_case",
                "authors": "Janssens",
                "specimen": "NoSuchSpecimen",
                "field": "Number of cycles to failure",
                "from": 1,
                "to": 2,
            },
            "need 1",
        ),
        (
            {
                "op": "set_case",
                "authors": "Janssens",
                "specimen": "Uniaxial_1",
                "field": "NoSuchColumn",
                "from": 1,
                "to": 2,
            },
            "no column headed",
        ),
        (
            {"op": "add_case", "row": 2, "values": {"Specimen name": "Dup"}},
            "is not empty",
        ),
        ({"op": "delete_case"}, "unknown overlay operation"),
    ],
)
def test_apply_one_rejects_invalid_operations(op, match):
    wb = load_workbook(FOLDER / "8_AISI316L.xlsx", data_only=True)
    with pytest.raises(OverlayError, match=match):
        _apply_one(wb, op)


def test_load_overlay_rejects_unsupported_schema_version(tmp_path):
    bad = tmp_path / "bad.overlay.json"
    bad.write_text(json.dumps({"schema_version": 99, "operations": []}), encoding="utf-8")
    with pytest.raises(OverlayError, match="unsupported schema_version"):
        load_overlay(bad)


def test_duplicate_base_hash_overlays_abort(tmp_path, monkeypatch):
    same_hash = "0" * 64
    for name in ("a.overlay.json", "b.overlay.json"):
        (tmp_path / name).write_text(
            json.dumps({"schema_version": 1, "base_sha256": same_hash, "operations": []}),
            encoding="utf-8",
        )
    monkeypatch.setattr(overlays_mod, "OVERLAY_DIR", tmp_path)
    overlays_mod._overlays_by_base_hash.cache_clear()
    try:
        with pytest.raises(OverlayError, match="name the same base file"):
            overlays_mod._overlays_by_base_hash()
    finally:
        overlays_mod._overlays_by_base_hash.cache_clear()


def test_renamed_copy_keeps_its_overlay(tmp_path):
    copy = tmp_path / "some_other_name.xlsx"
    copy.write_bytes((FOLDER / "8_AISI316L.xlsx").read_bytes())
    (block,) = load_material_workbook(copy).values()
    assert block["overlay"]["name"] == "8_AISI316L.overlay.json"
    assert len(block["cases"]) == 39


def test_modified_workbook_is_not_overlaid(tmp_path):
    wb = load_workbook(FOLDER / "8_AISI316L.xlsx", data_only=True)
    copy = tmp_path / "8_AISI316L.xlsx"
    wb.save(copy)  # same content, different bytes: no overlay is written against it
    assert overlay_path(copy) is None
    with pytest.warns(RuntimeWarning, match="overlay was not applied"):
        assert apply_overlay_to_workbook(wb, copy) is None

    report = ValidationReport(file=copy.name)
    with pytest.warns(RuntimeWarning, match="overlay was not applied"):
        load_material_workbook(copy, report=report)
    assert any("overlay was not applied" in warning for warning in report.warnings)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        load_material_workbook(copy, overlay=False)
    assert not [w for w in caught if issubclass(w.category, RuntimeWarning)]

    renamed = tmp_path / "other_name.xlsx"
    wb.save(renamed)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert apply_overlay_to_workbook(wb, renamed) is None
    assert not [w for w in caught if issubclass(w.category, RuntimeWarning)]


def test_experimental_data_manifest_is_current():
    manifest_path = FOLDER / "MANIFEST.json"
    assert manifest_path.is_file(), "Experimental_Data/MANIFEST.json is missing"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    files = workbook_files(FOLDER)
    assert set(manifest) == {p.name for p in files}
    for path in files:
        entry = manifest[path.name]
        assert entry["sha256"] == _sha256(path), (
            f"{path.name}: SHA-256 differs from MANIFEST.json; "
            "supplied workbooks must not be edited in place"
        )
        (block,) = load_material_workbook(path).values()
        assert entry["dataset_fingerprint"] == dataset_fingerprint(
            block
        ), f"{path.name}: dataset_fingerprint differs from MANIFEST.json"


# The five uniaxial and one non-proportional PSI rows were added from the
# publication using nominal (uncorrected) stresses. The 2026-10-29 correction
# applies only to their axial amplitudes, never to torsional amplitudes.
@pytest.mark.parametrize(
    ("specimen", "original", "corrected"),
    [
        ("Uniaxial_10", 240, 277.4),
        ("Uniaxial_11", 260, 300.5),
        ("Uniaxial_12", 290, 335.2),
        ("Uniaxial_13", 300, 346.7),
        ("Uniaxial_14", 330, 381.4),
        ("Non_prop_12", 233, 269.7),
    ],
)
def test_aisi316l_psi_axial_corrections_are_loaded(specimen, original, corrected):
    raw = json.loads(overlay_path(FOLDER / "8_AISI316L.xlsx").read_text(encoding="utf-8"))
    op = next(
        item
        for item in raw["operations"]
        if item["op"] == "add_case" and item["values"]["Specimen name"] == specimen
    )
    metadata = op["psi_axial_stress_correction"]
    assert (
        metadata.get("uncorrected_amplitude_MPa", metadata.get("published_rounded_amplitude_MPa"))
        == original
    )
    assert metadata["multiplier"] == pytest.approx(1.1558)
    assert op["values"]["Amplitute AXIAL stress  (MPa)"] == pytest.approx(corrected)
    source_amplitude = metadata.get("original_record_amplitude_MPa", original)
    assert corrected == pytest.approx(round(source_amplitude * 1.1558, 1))

    block = load_material_workbook(FOLDER / "8_AISI316L.xlsx")["AISI316L"]
    cases = [case for case in block["cases"].values() if case["meta"]["specimen_id"] == specimen]
    assert len(cases) == 1
    case = cases[0]
    assert case["harmonics"]["S_sin"][1, 1] == pytest.approx(corrected)
    assert case["harmonics"]["S_cos"][1, 1] == pytest.approx(0.0)
    assert case[1][1]["stress"][1, 1] == pytest.approx(corrected)
    assert len(block["cases"]) == 39
    assert len(block["runouts"]) == 13


@pytest.mark.parametrize(
    ("specimen", "shear"),
    [
        ("Torsional_8", 150.1),
        ("Torsional_9", 155.9),
        ("Non_prop_12", 134.7),
    ],
)
def test_aisi316l_psi_torsional_amplitudes_are_unchanged(specimen, shear):
    raw = json.loads(overlay_path(FOLDER / "8_AISI316L.xlsx").read_text(encoding="utf-8"))
    op = next(
        item
        for item in raw["operations"]
        if item["op"] == "add_case" and item["values"]["Specimen name"] == specimen
    )
    assert op["values"]["Amplitute SHEAR stress  (MPa)"] == pytest.approx(shear)
    block = load_material_workbook(FOLDER / "8_AISI316L.xlsx")["AISI316L"]
    case = next(case for case in block["cases"].values() if case["meta"]["specimen_id"] == specimen)
    sin_shear = case["harmonics"]["S_sin"][0, 1]
    cos_shear = case["harmonics"]["S_cos"][0, 1]
    assert (sin_shear**2 + cos_shear**2) ** 0.5 == pytest.approx(shear)
    if specimen == "Non_prop_12":
        assert case["meta"]["phase_tor_deg"] == pytest.approx(90.0)


def test_aisi316l_psi_correction_scope_is_six_of_eight_added_cases():
    overlay = json.loads(overlay_path(FOLDER / "8_AISI316L.xlsx").read_text(encoding="utf-8"))
    added = [op for op in overlay["operations"] if op["op"] == "add_case"]
    corrected = [op for op in added if "psi_axial_stress_correction" in op]
    assert len(added) == 8
    assert len(corrected) == 6
    assert {op["values"]["Specimen name"] for op in added} - {
        op["values"]["Specimen name"] for op in corrected
    } == {"Torsional_8", "Torsional_9"}
    nonprop = next(op for op in corrected if op["values"]["Specimen name"] == "Non_prop_12")
    assert nonprop["psi_axial_stress_correction"]["original_record_amplitude_MPa"] == 233.35
    assert nonprop["psi_axial_stress_correction"]["applied_amplitude_MPa"] == 269.7
    assert nonprop["psi_shear_source_precision"]["original_record_amplitude_MPa"] == 134.7
    for spec in ("Torsional_8", "Torsional_9"):
        op = next(item for item in added if item["values"]["Specimen name"] == spec)
        assert "psi_axial_stress_correction" not in op
        assert op["psi_shear_source_precision"]["correction_type"].startswith(
            "published figure rounding"
        )


def test_reis_fatigue_limits_are_the_authors_trend_lines_at_1e6():
    """d100: Sigm1, Taum1 of 42CrMo4_Reis = Anes et al. IJF 62 (2014) Table 5 at N = 1e6."""
    block = load_material_workbook(FOLDER / "14_42CrMo4_Reis.xlsx")["42CrMo4_Reis"]
    assert block["material"]["Tensile_fatigue_limit_MPa"] == 452.0
    assert block["material"]["Torsional_fatigue_limit_MPa"] == 372.0
    (src,) = block["fatigue_limit_sources"]["Sigm1"]
    assert src["reference"] == "10.1016/j.ijfatigue.2013.05.004"
    assert src["reference_cycles"] == 1e6 and src["loading"] == "push-pull"
    # the Table 4 trend lines behind Table 5 (sigma_a = 1204.3 N^-0.071, tau_a = 864.78 N^-0.061)
    assert abs(1204.3 * 1e6**-0.071 - 452.0) < 1.0
    assert abs(864.78 * 1e6**-0.061 - 372.0) < 1.0
    off = load_material_workbook(FOLDER / "14_42CrMo4_Reis.xlsx", overlay=False)["42CrMo4_Reis"]
    assert off["material"]["Tensile_fatigue_limit_MPa"] is None


def test_ti6al4v_summary_keeps_only_the_meneghetti_batch():
    """d108: every failure is Meneghetti's; the Delahay 2006 row no longer enters the values."""
    block = load_material_workbook(FOLDER / "1_Ti6Al4V_Meneghetti.xlsx")["Ti6Al4V_Meneghetti"]
    m = block["material"]
    assert (m["Tensile_strength_MPa"], m["Yield_strength_MPa"]) == (980.0, 915.0)
    assert (m["Tensile_fatigue_limit_MPa"], m["Torsional_fatigue_limit_MPa"]) == (500.0, 422.5)
    for key in ("Sigm1", "Taum1"):
        assert [s["authors"] for s in block["fatigue_limit_sources"][key]] == ["Meneghetti et al."]
    # d109: the 4 Delahay rows (endurance limits entered as run-outs) leave the dataset
    assert (len(block["cases"]), len(block["runouts"])) == (10, 0)
    off = load_material_workbook(FOLDER / "1_Ti6Al4V_Meneghetti.xlsx", overlay=False)[
        "Ti6Al4V_Meneghetti"
    ]
    assert off["material"]["Tensile_strength_MPa"] == 1035.0  # supplied average of both rows
    assert (len(off["cases"]), len(off["runouts"])) == (10, 4)


def _drop(row=58, **check):
    return {"op": "drop_case", "row": row, "authors": "Delahay T", "specimen": "SPECIMEN_5",
            "check": check}  # fmt: skip


def test_drop_case_clears_the_checked_row_only():
    wb = load_workbook(FOLDER / "1_Ti6Al4V_Meneghetti.xlsx", data_only=True)
    ws = wb["Dataset"]
    before = [[c.value for c in r] for r in ws.iter_rows(min_row=57, max_row=59)]
    _apply_one(wb, _drop(**{"Amplitute AXIAL stress  (MPa)": 652}))
    assert all(c.value is None for c in ws[58])
    assert [c.value for c in ws[57]] == before[0] and [c.value for c in ws[59]] == before[2]


@pytest.mark.parametrize(
    ("op", "match"),
    [
        (_drop(row=57), "expected 'Delahay T'"),
        (_drop(**{"Amplitute AXIAL stress  (MPa)": 583}), "expected 583"),
        (_drop(**{"No such column": 1}), "no column headed"),
    ],
)
def test_drop_case_rejects_a_row_that_does_not_match(op, match):
    wb = load_workbook(FOLDER / "1_Ti6Al4V_Meneghetti.xlsx", data_only=True)
    with pytest.raises(OverlayError, match=match):
        _apply_one(wb, op)


def test_al7075_limits_at_the_pedranz_basis_and_measured_plastic_terms():
    """d111: Pedranz 2025 Appendix B Table 5 limits at 5e6; Benedetti 2021 Table 4 eps'_f, c."""
    block = load_material_workbook(FOLDER / "Al7075_T6.xlsx")["Al7075_T6"]
    m = block["material"]
    assert (m["Tensile_fatigue_limit_MPa"], m["Torsional_fatigue_limit_MPa"]) == (143.8, 100.9)
    (sigm1,) = block["fatigue_limit_sources"]["Sigm1"]
    assert sigm1["reference_cycles"] == 5e6 and sigm1["loading"] == "push-pull"
    assert sigm1["reference"] == "10.1016/j.ijmecsci.2025.109961"
    assert (m["epsilon_f"], m["c_fat"], m["c0"]) == (0.131, -0.674, -0.674)
    assert m["gamma_f"] == pytest.approx(3**0.5 * 0.131, abs=1e-6)
    # sigma'_f and b stay the Basquin fit of the HCF tests (the LCF values extrapolate 2x high)
    assert m["sigma_f"] is None and m["b_fat"] is None
    assert not block["identity_issues"]
    off = load_material_workbook(FOLDER / "Al7075_T6.xlsx", overlay=False)["Al7075_T6"]
    assert off["material"]["Tensile_fatigue_limit_MPa"] == 150.0
