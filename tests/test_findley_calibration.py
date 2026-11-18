"""Configurable Findley calibration amplitude (MCC or MCE) in multifat-calibrate.

Port of the idea of upstream tests/test_findley_calibration.py (A. Chiocca, 36ff0a1, 0105b69).
Not ported: the versioned calibration store (JSON history, signatures, fingerprints) and the
Script 0 integration; the MCC choice is the existing ``longest_chord`` plane shear rule, so the
calibration uses the established plane search and is checked against the isolated harmonic
evaluator of ``fin_harmonic``.
"""

import json

import numpy as np
import pytest

from test_findley_harmonic import ranking_case
from workbook_factory import write_workbook

from multifat.cli import calibrate as cal
from multifat.cli.pipeline import load_constants
from multifat.eval.evaluator import best_cp_for_case, shear_rule_of
from multifat.eval.findley_calibration import (
    FindleyCalibrationUnavailable,
    amplitude_identity,
    findley_method,
    require_stress_harmonics,
    resolve_findley_algorithm,
)
from multifat.imp import experimental_data as exp
from multifat.models.cp_methods.fin import Findley
from multifat.models.cp_methods.fin_harmonic import evaluate_findley
from multifat.models.utils.shear_amplitude import DEFINITION_VERSION


@pytest.mark.parametrize(
    "selection,expected",
    [
        ("method_default", "mce_harmonic"),
        ("rss", "mce_harmonic"),
        ("mce_harmonic", "mce_harmonic"),
        ("mcc_harmonic", "mcc_harmonic"),
    ],
)
def test_resolution(selection, expected):
    assert resolve_findley_algorithm(selection) == expected


@pytest.mark.parametrize("bad", [None, [], "mcc", ""])
def test_unknown_selection_rejected(bad):
    with pytest.raises(ValueError):
        resolve_findley_algorithm(bad)


def test_identity_records_resolved_algorithm_and_path():
    assert amplitude_identity("method_default") == {
        "algorithm": "mce_harmonic",
        "definition_version": DEFINITION_VERSION,
        "evaluation_path": "legacy",
    }
    assert amplitude_identity("rss")["algorithm"] == "mce_harmonic"
    assert amplitude_identity("rss")["evaluation_path"] == "harmonic"
    assert amplitude_identity("mcc_harmonic")["algorithm"] == "mcc_harmonic"


def test_default_and_mce_use_the_production_findley_rule():
    for selection in ("method_default", "mce_harmonic", "rss"):
        method = findley_method(selection)
        assert type(method) is Findley and shear_rule_of(method) == "rss"
    mcc = findley_method("mcc_harmonic")
    assert mcc.name == "FIN" and shear_rule_of(mcc) == "longest_chord"
    assert shear_rule_of(Findley()) == "rss"  # the production class is untouched


def _best(case, method, planes):
    harmonics = case["harmonics"]
    return best_cp_for_case(
        S0=case["S0"],
        E0=case["E0"],
        S1=case["S1"],
        E1=case["E1"],
        R_list=planes,
        method=method,
        params=case["params"],
        S_sin=harmonics["S_sin"],
        S_cos=harmonics["S_cos"],
        refine=False,
    )


@pytest.mark.parametrize("proportional", [True, False])
@pytest.mark.parametrize(
    "selection",
    ["mcc_harmonic", "mce_harmonic", "rss"],
)
def test_plane_search_variant_matches_isolated_harmonic_evaluator(proportional, selection):
    case = ranking_case(proportional)
    best = _best(case, findley_method(selection), case["R_list"])
    expected = evaluate_findley(**case, shear_amplitude_algorithm=selection).values[""]
    np.testing.assert_allclose(best.ext, expected, rtol=3e-15)


def test_mcc_and_mce_rank_planes_differently():
    case = ranking_case(False)
    winners = {}
    for selection in ("mcc_harmonic", "mce_harmonic"):
        method = findley_method(selection)
        values = [_best(case, method, case["R_list"][i : i + 1]).ext for i in range(2)]
        winners[selection] = int(np.argmax(values))
    assert winners == {"mcc_harmonic": 1, "mce_harmonic": 0}


def test_require_stress_harmonics_rejects_missing_or_invalid():
    good = {"harmonics": {"S_sin": np.eye(3), "S_cos": np.eye(3)}}
    require_stress_harmonics(good, "Mat", 1)
    for case in (
        {},
        {"harmonics": None},
        {"harmonics": {"S_sin": np.eye(3)}},
        {"harmonics": {"S_sin": np.eye(3), "S_cos": np.full((3, 3), np.nan)}},
        {"harmonics": {"S_sin": np.eye(3), "S_cos": np.zeros((2, 2))}},
    ):
        with pytest.raises(FindleyCalibrationUnavailable, match="Mat/FIN case 7"):
            require_stress_harmonics(case, "Mat", 7)


# ---- calibration -----------------------------------------------------------------------


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    path = write_workbook(tmp_path_factory.mktemp("fin") / "Synth.xlsx")
    data = exp.load_material_workbook(path)
    return data, exp.material_id(path.stem)


def _calibrate_fin(data, material, **kwargs):
    return cal.calibrate_material_params(
        data=data,
        material_name=material,
        R_list=cal._build_R_list(n1=3, n2=6),
        blocks={"FI"},
        n_grid_1d=5,
        n_refine_1d=1,
        **kwargs,
    )


def test_default_calibration_is_unchanged_and_records_its_amplitude(synthetic):
    data, material = synthetic
    default = _calibrate_fin(data, material)
    explicit = _calibrate_fin(data, material, shear_amplitude_algorithm="method_default")
    assert default["FI"]["k_FI"] == explicit["FI"]["k_FI"]
    assert default["FI"]["scatter_STD"] == explicit["FI"]["scatter_STD"]
    assert default["amplitude_settings"] == {"FIN": amplitude_identity("method_default")}


@pytest.mark.parametrize("selection", ["mce_harmonic", "rss", "mcc_harmonic"])
def test_explicit_amplitude_on_proportional_calibration_failures(synthetic, selection):
    # Axial and torsion failures are proportional: MCC = MCE, so the fit coincides.
    data, material = synthetic
    default = _calibrate_fin(data, material)
    chosen = _calibrate_fin(data, material, shear_amplitude_algorithm=selection)
    assert chosen["FI"]["k_FI"] == pytest.approx(default["FI"]["k_FI"], abs=1e-9)
    assert chosen["FI"]["scatter_STD"] == pytest.approx(default["FI"]["scatter_STD"], rel=1e-9)
    assert chosen["amplitude_settings"]["FIN"] == amplitude_identity(selection)


def test_mcc_selection_searches_planes_with_the_longest_chord_rule(synthetic, monkeypatch):
    data, material = synthetic
    seen = []
    real = cal.shear_rule_of

    def spy(method):
        rule = real(method)
        seen.append((method.name, rule))
        return rule

    monkeypatch.setattr(cal, "shear_rule_of", spy)
    _calibrate_fin(data, material, shear_amplitude_algorithm="mcc_harmonic")
    assert set(seen) == {("FIN", "longest_chord")}
    seen.clear()
    _calibrate_fin(data, material)
    assert set(seen) == {("FIN", "rss")}


def test_explicit_amplitude_requires_stored_harmonics(synthetic):
    data, material = synthetic
    broken = {
        material: {
            **data[material],
            "cases": {
                cid: {**case, "harmonics": None} for cid, case in data[material]["cases"].items()
            },
        }
    }
    with pytest.raises(FindleyCalibrationUnavailable, match="S_sin"):
        _calibrate_fin(broken, material, shear_amplitude_algorithm="mcc_harmonic")
    # the default keeps the proportional endpoint fallback
    assert _calibrate_fin(broken, material)["FI"]["k_FI"] > 0


def test_unknown_algorithm_rejected(synthetic):
    data, material = synthetic
    with pytest.raises(ValueError, match="Unknown Findley"):
        _calibrate_fin(data, material, shear_amplitude_algorithm="mcc")


def test_explicit_amplitude_applies_to_every_calibrated_plane_criterion(synthetic):
    data, material = synthetic
    out = cal.calibrate_material_params(
        data=data,
        material_name=material,
        R_list=cal._build_R_list(n1=3, n2=6),
        blocks={"FS"},
        n_grid_1d=5,
        n_refine_1d=1,
        shear_amplitude_algorithm="mcc_harmonic",
    )
    # Since the option covers every critical-plane criterion, k_FS records its amplitude too.
    assert out["FI"] is None
    assert out["amplitude_settings"]["FS"]["algorithm"] == "mcc_harmonic"
    assert set(out["amplitude_settings"]) == {"FS"}


def test_cli_accepts_the_option(tmp_path, monkeypatch, capsys):
    path = write_workbook(tmp_path / "Cli.xlsx")
    calls = []
    real = cal.calibrate_material_params

    def fake(**kwargs):
        calls.append(kwargs["shear_amplitude_algorithm"])
        return real(**{**kwargs, "n_grid_1d": 5, "n_refine_1d": 1, "plane_stride": 20})

    monkeypatch.setattr(cal, "calibrate_material_params", fake)
    cal.main(["--data-folder", str(path.parent), "--shear-amplitude-algorithm", "mcc_harmonic"])
    assert calls == ["mcc_harmonic"]
    assert "mcc_harmonic" in capsys.readouterr().out


# ---- pipeline compatibility ------------------------------------------------------------


def _record(tmp_path, amplitude):
    rec = {
        "schema_version": 1,
        "material_id": "Synth",
        "calibration_id": "20261020T000000Z_abc",
        "parameters": {"k_FI": 0.4},
    }
    if amplitude is not None:
        rec["amplitude_settings"] = {"FIN": {"algorithm": amplitude}}
    path = tmp_path / "r.json"
    path.write_text(json.dumps(rec), encoding="utf-8")
    return path


@pytest.mark.parametrize("amplitude", [None, "mce_harmonic"])
def test_record_with_mce_findley_is_accepted(tmp_path, amplitude):
    params, _ = load_constants(_record(tmp_path, amplitude), "Synth")
    assert params == {"k_FI": 0.4}


def test_record_with_mcc_findley_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="mcc_harmonic"):
        load_constants(_record(tmp_path, "mcc_harmonic"), "Synth")
