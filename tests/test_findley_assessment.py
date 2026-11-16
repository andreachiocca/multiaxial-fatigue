"""Explicit Findley shear-amplitude assessment in multifat-run (FIN and FIN_ext).

Rewritten against this repository's architecture: the Findley constant is calibrated inside the
run with the same criterion it is assessed with, so the calibration and assessment options cannot
disagree. The upstream calibration store is not carried over.
"""

import hashlib
import json
import warnings
from pathlib import Path

import pandas as pd
import pytest

import multifat.imp.experimental_data as exp
from multifat.cli.compare import collect_all_results
from multifat.cli.run import (
    harmonic_case_or_skip,
    shear_amplitude_for,
    write_assessment_metadata,
)
from multifat.eval.findley_calibration import FindleyCalibrationUnavailable, findley_method
from multifat.models.cp_methods.fin import Findley
from multifat.report.export_csv import save_fatigue_csv

ROOT = Path(__file__).resolve().parents[1]


def test_default_selection_is_the_production_path():
    assert shear_amplitude_for("FIN", "method_default") is None
    assert type(findley_method("method_default")) is Findley


def test_explicit_selection_resolves_algorithm_and_path():
    rss = shear_amplitude_for("FIN", "rss")
    assert rss["algorithm"] == "mce_harmonic"
    assert rss["evaluation_path"] == "harmonic"
    assert rss["selected_option"] == "rss"
    mcc = shear_amplitude_for("FIN", "mcc_harmonic")
    assert mcc["algorithm"] == "mcc_harmonic"


def test_explicit_mcc_uses_the_mcc_criterion():
    from multifat.eval.findley_calibration import FindleyMCC

    assert isinstance(findley_method("mcc_harmonic"), FindleyMCC)


@pytest.mark.parametrize("model", ["SWT", "BP", "MATAKE"])
def test_every_model_accepts_a_valid_selection_and_rejects_an_unknown_one(model):
    """The selection is no longer Findley-only (multifat.eval.shear_amplitude_option); a model
    without a selectable amplitude (BP) runs unchanged, with a warning, in evaluate_model."""
    assert shear_amplitude_for(model, "mcc_harmonic")["algorithm"] == "mcc_harmonic"
    assert shear_amplitude_for(model, "longest_chord")["algorithm"] == "mcc_harmonic"
    with pytest.raises(ValueError):
        shear_amplitude_for(model, "mcc")


def test_cases_without_stored_harmonics_are_skipped_with_a_warning():
    case = {"harmonics": {}}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert harmonic_case_or_skip(case, "AISI316L", "FIN", 3) is False
    assert any(issubclass(w.category, RuntimeWarning) for w in caught)
    assert "case 3" in str(caught[0].message)


def test_case_with_stored_harmonics_is_kept():
    data = exp.import_experimental_data_grouped_by_material(ROOT / "Experimental_Data")
    case_id, case = next(iter(data["AISI316L"]["cases"].items()))
    assert harmonic_case_or_skip(case, "AISI316L", "FIN", case_id) is True


def test_missing_harmonics_raise_in_the_strict_check():
    from multifat.eval.findley_calibration import require_stress_harmonics

    with pytest.raises(FindleyCalibrationUnavailable, match="S_sin and S_cos"):
        require_stress_harmonics({"harmonics": {}}, "AISI316L", 7)


def _write_csv(tmp_path: Path, result_label: str | None) -> Path:
    return Path(
        save_fatigue_csv(
            out_dir=str(tmp_path),
            method_name="FIN",
            material_name="AISI316L",
            N_exp=[1.0e5, 2.0e5],
            CP=[1.0, 1.1],
            N_pred=[1.0e5, 2.1e5],
            loading_type=["UNIAXIAL", "UNIAXIAL"],
            result_label=result_label,
        )
    )


def test_result_label_names_the_file_and_the_method_column(tmp_path):
    path = _write_csv(tmp_path, "FIN.mcc_harmonic")
    assert path.name == "FIN.mcc_harmonic_AISI316L.csv"
    frame = pd.read_csv(path)
    assert set(frame["Method"]) == {"FIN.mcc_harmonic"}


def test_default_csv_keeps_the_method_name(tmp_path):
    path = _write_csv(tmp_path, None)
    assert path.name == "FIN_AISI316L.csv"
    assert set(pd.read_csv(path)["Method"]) == {"FIN"}


def test_metadata_sidecar_records_the_selection(tmp_path):
    path = _write_csv(tmp_path, "FIN.mcc_harmonic")
    amplitude = shear_amplitude_for("FIN", "mcc_harmonic")
    sidecar = write_assessment_metadata(
        path, method="FIN", result_label="FIN.mcc_harmonic", amplitude=amplitude
    )
    assert sidecar.name == "FIN.mcc_harmonic_AISI316L.meta.json"
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    assert meta["selected_option"] == "mcc_harmonic"
    assert meta["resolved_algorithm"] == "mcc_harmonic"
    assert meta["evaluation_path"] == "harmonic"
    assert meta["csv_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_compare_reader_finds_the_labelled_results(tmp_path):
    _write_csv(tmp_path, "FIN.mcc_harmonic")
    frame = collect_all_results(tmp_path)
    assert "FIN.mcc_harmonic" in set(frame["Method"])
