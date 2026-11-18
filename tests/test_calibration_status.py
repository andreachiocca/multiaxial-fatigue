"""Calibration failure handling: the named infeasibility, batch isolation and statuses."""

from __future__ import annotations

import json

import pytest

import multifat.cli.calibrate as calibrate_module
from multifat.cli.calibrate import (
    CalibrationInfeasible,
    _make_plane_frames,
    calibrate_many,
    calibrate_material_params,
    calibrate_workbook,
    main,
)
from multifat.cli.pipeline import run_workbook
from multifat.imp import experimental_data as exp
from workbook_factory import (
    SUMMARY_HEADERS,
    SUMMARY_VALUES,
    basquin_rows,
    write_workbook,
)

LIVES = (1e5, 2e5, 3e5, 4e5)
TENSION = (400.0, 350.0, 300.0, 260.0)
TORSION = (200.0, 320.0, 440.0, 560.0)


def rising_damage_rows():
    """Tension and torsion rows on which the LI damage rises with life (wrong ranking)."""
    rows = [[N, "No", "No", 0, 0, 0, 0, 0, s, 0, 0] for N, s in zip(LIVES, TENSION)]
    rows += [[N, "No", "No", 0, 0, 0, 0, 0, 0, 0, t] for N, t in zip(LIVES, TORSION)]
    return rows


def few_case_rows():
    """Two tension and two torsion tests on Basquin lines: enough to calibrate, cheap to run."""
    rows = [[N, "No", "No", 0, 0, 0, 0, 0, 1500.0 * (2 * N) ** -0.09, 0, 0] for N in (1e4, 1e6)]
    rows += [[N, "No", "No", 0, 0, 0, 0, 0, 0, 0, 860.0 * (2 * N) ** -0.07] for N in (2e4, 1.5e6)]
    return rows


def write_rising_damage_workbook(path):
    """The synthetic infeasible workbook, without the k_FS and k_FI Summary columns."""
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if not h.startswith(("k FS", "k FI"))]
    return write_workbook(
        path,
        rows=rising_damage_rows(),
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )


def _report(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    return {r.split("\t")[0]: dict(zip(head, r.split("\t"))) for r in lines[1:]}


def test_infeasible_li_is_named_with_its_slope_range(tmp_path):
    wb = write_workbook(tmp_path / "Bad.xlsx", rows=rising_damage_rows())
    data = exp.load_material_workbook(wb)
    with pytest.raises(CalibrationInfeasible) as caught:
        calibrate_material_params(
            data=data,
            material_name=exp.material_id(wb.stem),
            R_list=_make_plane_frames(),
            plane_stride=10,
            blocks={"LI"},
        )
    exc = caught.value
    assert exc.block == "LI"
    assert exc.parameter == "L_LI"
    assert exc.bounds == (0.0, 20.0)
    assert str(exc).startswith("No finite calibration objective")
    assert "DP increases with life for all L_LI in [0, 20]" in str(exc)
    assert "slope +" in exc.diagnostic
    assert str(exc).endswith(exc.diagnostic)


def test_calibrate_many_isolates_an_infeasible_material(tmp_path):
    write_workbook(tmp_path / "1_Good1.xlsx", rows=few_case_rows())
    write_rising_damage_workbook(tmp_path / "2_Bad.xlsx")
    write_workbook(tmp_path / "3_Good2.xlsx", rows=few_case_rows())

    records = calibrate_many(exp_folder=tmp_path, plane_stride=10)

    assert list(records) == ["Good1", "Bad", "Good2"]
    assert records["Good1"]["status"] == "calibrated"
    assert records["Good2"]["status"] == "calibrated"
    assert records["Good1"]["FS"]["k_FS"] > 0.0 and records["Good1"]["reason"] == ""
    assert records["Bad"]["status"] == "infeasible"
    assert records["Bad"]["infeasible"]["block"] == "FS"
    assert str(records["Bad"]["reason"]).startswith("No finite calibration objective")


def test_calibrate_main_reports_every_material_and_exits_one(tmp_path, capsys):
    write_workbook(tmp_path / "1_Good1.xlsx", rows=few_case_rows())
    write_rising_damage_workbook(tmp_path / "2_Bad.xlsx")
    write_workbook(tmp_path / "3_Good2.xlsx", rows=few_case_rows())

    with pytest.raises(SystemExit) as caught:
        main(["--data-folder", str(tmp_path)])
    assert caught.value.code == 1

    out = capsys.readouterr().out
    assert "Calibration summary" in out
    assert "Bad: infeasible at FS" in out
    summary = out.split("Calibration summary", 1)[1]
    assert "Good1" in summary and "Good2" in summary
    assert summary.count("calibrated") == 2
    assert summary.count("infeasible") == 1


def test_one_infeasible_block_makes_the_material_infeasible_without_partial_constants(
    tmp_path, monkeypatch
):
    """Pins the deferred block isolation: an infeasible LI block must not leave the other four.

    No natural workbook isolates an LI-only failure, so the LI fit is forced to fail while the
    other blocks run for real. Enabling isolation would keep k_FS, k_FI, k_MGSE and a_MSWT and
    fail here, which must then be a deliberate change.
    """
    wb = write_workbook(tmp_path / "1_Synth.xlsx", rows=few_case_rows())
    others = {"FS", "FI", "MGSE_ZHU", "MSWT"}
    alone = calibrate_material_params(
        data=exp.load_material_workbook(wb),
        material_name="Synth",
        R_list=_make_plane_frames(),
        plane_stride=10,
        blocks=others,
    )
    assert all(alone[b] is not None for b in others)

    real_fit = calibrate_module._fit_single_param

    def fit_with_li_infeasible(**kwargs):
        if kwargs["block"] == "LI":
            raise CalibrationInfeasible(
                "No finite calibration objective",
                block="LI",
                parameter="L_LI",
                bounds=(0.0, 5.0),
                diagnostic="forced for the test",
            )
        return real_fit(**kwargs)

    monkeypatch.setattr(calibrate_module, "_fit_single_param", fit_with_li_infeasible)

    record = calibrate_workbook(wb)

    assert record["status"] == "infeasible"
    assert record["infeasible"]["block"] == "LI"
    for name in ("k_FS", "k_FI", "k_MGSE", "a_MSWT"):
        assert name not in record["params"]


def test_fatigue_limit_and_notched_only_are_not_applicable(tmp_path, capsys):
    limits = [list(r) for r in basquin_rows()]
    for row in limits:
        row[1] = "Yes"
    write_workbook(tmp_path / "1_Limits.xlsx", rows=limits)
    notched = [list(r) for r in basquin_rows()]
    for row in notched:
        row[2] = "Yes"
    write_workbook(tmp_path / "2_Notched.xlsx", rows=notched)

    records = calibrate_many(exp_folder=tmp_path)
    assert records["Limits"]["status"] == "not_applicable"
    assert "fatigue-limit dataset" in records["Limits"]["reason"]
    assert records["Notched"]["status"] == "not_applicable"
    assert "notched-only" in records["Notched"]["reason"]

    record = calibrate_workbook(tmp_path / "1_Limits.xlsx")
    assert record["status"] == "not_applicable"
    assert record["errors"] == []

    main(["--data-folder", str(tmp_path)])  # no error status: exits 0
    out = capsys.readouterr().out
    assert out.count("not_applicable") >= 2


def test_pipeline_marks_the_default_after_an_infeasible_calibration(tmp_path):
    wb = write_rising_damage_workbook(tmp_path / "Synth.xlsx")
    out = tmp_path / "res"
    outcomes = run_workbook(
        wb,
        models=["FS", "FIN", "LI", "MGSE_ZHU", "MSWT", "MDC"],
        results_dir=out,
    )
    by_model = {outcome.model: outcome for outcome in outcomes}

    assert by_model["FS"].status == "skipped" and "calibration infeasible at FS" in (
        by_model["FS"].reason
    )
    assert by_model["FIN"].status == "skipped" and "k_FI: calibration infeasible at FS" in (
        by_model["FIN"].reason
    )
    # The defaults would run, but the 4 tension tests over 0.6 decades of this synthetic workbook
    # do not identify a master curve, so the models are skipped (issue #74); the annotation of a
    # default constant is tested on its own below.
    for name in ("LI", "MGSE_ZHU", "MSWT"):
        assert by_model[name].status == "skipped"
        assert "master curve not identifiable" in by_model[name].reason

    record = json.loads((out / "Synth_calibration.json").read_text(encoding="utf-8"))
    assert record["status"] == "infeasible"
    assert record["infeasible"]["block"] == "FS"
    # Pinned: the constants the all-or-nothing calibration still produced before this change.
    assert record["params"] == pytest.approx(
        {
            "epsilon_f": 0.4890412621359223,
            "c_fat": -0.58,
            "gamma_f": 0.8470443130170272,
            "c0": -0.58,
        },
        rel=1e-12,
    )

    report = _report(out / "Synth_run_report.tsv")
    assert "master curve not identifiable" in report["LI"]["reason"]
    assert "calibration infeasible at FS" in report["FS"]["reason"]


def test_a_model_running_on_a_lost_constant_default_is_marked():
    from multifat.cli.pipeline import Outcome, _annotate_calibration
    from multifat.models.registry import get_models

    state = {"status": "infeasible", "block": "FS", "lost": ["L_LI"], "at_bound": {}}
    model = get_models(["LI"])[0]
    outcome = _annotate_calibration(Outcome("LI", "ok"), model, state)
    assert outcome.status == "ok (default constant)"
    assert "default 1.0 used" in outcome.notes and "calibration infeasible at FS" in outcome.notes


def test_constant_on_an_arbitrary_search_bound_is_not_used():
    """An upper search limit is not an optimum; a physical limit (a_MSWT = 1) is."""
    from multifat.cli.calibrate import stopped_by_search_bound

    assert stopped_by_search_bound("k_FI", 3.0, (0.0, 3.0))
    assert stopped_by_search_bound("L_LI", 20.0, (0.0, 20.0))
    assert not stopped_by_search_bound("a_MSWT", 1.0, (0.0, 1.0))
    assert not stopped_by_search_bound("k_FI", 0.0, (0.0, 3.0))  # k = 0 is admissible
    assert not stopped_by_search_bound("k_FI", 0.9, (0.0, 3.0))
