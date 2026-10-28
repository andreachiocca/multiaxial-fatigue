"""End-to-end pipeline on synthetic workbooks: skipped models, calibration fill-in, reports."""

from __future__ import annotations

import json

import numpy as np
import pytest

from multifat.cli.pipeline import run_workbook
from multifat.imp.schema import WorkbookError
from workbook_factory import SUMMARY_HEADERS, SUMMARY_VALUES, basquin_rows, write_workbook


def _report(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    return {r.split("\t")[0]: dict(zip(head, r.split("\t"))) for r in lines[1:]}


def test_missing_constants_are_skipped_and_others_run(tmp_path):
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if not h.startswith("Tensile strength")]
    wb = write_workbook(
        tmp_path / "Synth.xlsx",
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )
    out = tmp_path / "res"
    outcomes = run_workbook(wb, models=["FIN", "CROSSLAND", "SINES"], results_dir=out)
    status = {o.model: o.status for o in outcomes}
    assert status == {
        "FIN": "ok",
        "CROSSLAND": "ok",
        "SINES": "skipped",
    }
    rep = _report(out / "Synth_run_report.tsv")
    assert "['Su']" in rep["SINES"]["reason"]
    assert "FIN_ext" not in rep and "CROSSLAND_ext" not in rep
    assert (out / "FIN_Synth.csv").exists() and (out / "CROSSLAND_Synth.csv").exists()
    assert "kept: 15" in (out / "Synth_validation.txt").read_text(encoding="utf-8")
    cal = json.loads((out / "Synth_calibration.json").read_text(encoding="utf-8"))
    assert cal["params"] == {}  # nothing calibratable is missing
    prov = cal["provenance"]
    assert prov["failures"] == 15 and prov["runouts"] == 1 and len(prov["workbook_sha256"]) == 64


def test_missing_workbook_constant_is_calibrated(tmp_path):
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if not h.startswith("k FI")]
    wb = write_workbook(
        tmp_path / "NoKFI.xlsx",
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )
    out = tmp_path / "res"
    outcomes = run_workbook(wb, models=["FIN"], results_dir=out)
    assert {o.model: (o.status, o.calibrated) for o in outcomes} == {
        "FIN": ("ok", "k_FI"),
    }
    cal = json.loads((out / "NoKFI_calibration.json").read_text(encoding="utf-8"))
    assert 0.0 <= cal["params"]["k_FI"] <= 3.0


def test_explicit_ext_and_all_variants_output(tmp_path):
    """A requested extension is not accompanied by an unwanted plain result."""
    wb = write_workbook(tmp_path / "Variant.xlsx")
    ext_dir = tmp_path / "ext"
    result = run_workbook(wb, models=["FIN_ext"], results_dir=ext_dir)
    assert [(row.model, row.status) for row in result] == [("FIN_ext", "ok")]
    assert (ext_dir / "FIN_ext_Variant.csv").exists()
    assert not (ext_dir / "FIN_Variant.csv").exists()

    both_dir = tmp_path / "both"
    result = run_workbook(wb, models=["FIN", "CROSSLAND"], all_variants=True, results_dir=both_dir)
    assert {row.model: row.status for row in result} == {
        "FIN": "ok",
        "FIN_ext": "ok",
        "CROSSLAND": "ok",
    }
    assert (both_dir / "FIN_Variant.csv").exists()
    assert (both_dir / "FIN_ext_Variant.csv").exists()
    assert not (both_dir / "CROSSLAND_ext_Variant.csv").exists()


def test_invalid_explicit_ext_is_rejected(tmp_path):
    wb = write_workbook(tmp_path / "InvalidVariant.xlsx")
    with pytest.raises(ValueError, match="no _ext variant"):
        run_workbook(wb, models=["CROSSLAND_ext"], results_dir=tmp_path / "res")


def test_workbook_without_design_set_is_rejected(tmp_path):
    rows = [r for r in basquin_rows() if r[8] == 0 or r[10] != 0]
    wb = write_workbook(tmp_path / "NoUni.xlsx", rows=rows)
    with pytest.raises(WorkbookError, match="design set"):
        run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res")


def test_fitted_basquin_constants_feed_mdc_and_flag_their_tests(tmp_path):
    import pandas as pd

    wb = write_workbook(tmp_path / "Bq.xlsx", rows=basquin_rows(sf=1500, b=-0.09, tf=860, b0=-0.07))
    out = tmp_path / "res"
    outcomes = {o.model: o for o in run_workbook(wb, models=["MDC"], results_dir=out)}
    assert outcomes["MDC"].status == "ok"
    assert outcomes["MDC"].calibrated == "b0,b_fat,sigma_f,tau_f"
    df = pd.read_csv(out / "MDC_Bq.csv")
    tors = df[(df.Flag == 0) & (df.Loading_type == "Torsion")]
    assert len(tors) == 5  # the torsion tests behind tau_f, b0 are design points, not predictions
    assert (df[df.Flag == 1].Loading_type != "Torsion").all()
    # On its own calibration tests MDC's native life is exact (tau_r = tau_a in torsion).
    assert np.allclose(tors.Error_log10_native, 0.0, atol=1e-9)
    cal = json.loads((out / "Bq_calibration.json").read_text(encoding="utf-8"))
    assert cal["source"]["tau_f"].startswith("Basquin fit")


def test_flat_sn_data_skip_mdc_with_reason(tmp_path):
    rows = [r[:8] + [250.0] + r[9:] if r[8] and not r[10] else r for r in basquin_rows()]
    wb = write_workbook(tmp_path / "Flat.xlsx", rows=rows)
    o = {o.model: o for o in run_workbook(wb, models=["MDC"], results_dir=tmp_path / "res")}
    assert o["MDC"].status == "skipped" and "not identifiable" in o["MDC"].reason


def test_all_writes_eligibility_matrix(tmp_path):
    """Upstream 904c619 eligibility matrix, from the run reports of every workbook."""
    import pandas as pd

    from multifat.cli.pipeline import main

    good = write_workbook(tmp_path / "Good.xlsx")
    notched = [list(r) for r in basquin_rows()]
    for r in notched:
        r[2] = "Yes"
    write_workbook(tmp_path / "Notched.xlsx", rows=notched)
    runouts = [list(r) for r in basquin_rows()]
    for r in runouts:
        r[1] = "Yes"
    i = SUMMARY_HEADERS.index("Tensile fatigue limit* (MPa)")
    no_lim = SUMMARY_VALUES[:i] + [None] + SUMMARY_VALUES[i + 1 :]
    # run-outs only and no fatigue limit: neither a master curve nor a fatigue index
    write_workbook(tmp_path / "Runouts.xlsx", rows=runouts, summary_values=no_lim)
    nok = SUMMARY_HEADERS.index("k FI (-)")
    vals = SUMMARY_VALUES[:nok] + [None] + SUMMARY_VALUES[nok + 1 :]
    write_workbook(tmp_path / "NoK.xlsx", summary_values=vals)
    res = tmp_path / "res"
    books = [str(tmp_path / f"{m}.xlsx") for m in ("Good", "Notched", "Runouts", "NoK")]
    with pytest.raises(SystemExit):
        main(books + ["--results-dir", str(res), "--no-compare", "--no-calibrate",
                      "--models", "CROSSLAND,FIN"])  # fmt: skip
    m = pd.read_csv(res / "Comparison" / "eligibility_matrix.csv", index_col=0)
    assert m.loc["CROSSLAND", "Good"] == "available"
    assert m.loc["FIN", "NoK"] == "calibration_required"
    assert m.loc["CROSSLAND", "Notched"] == "out_of_scope"
    assert m.loc["CROSSLAND", "Runouts"] == "insufficient_inputs"
    avail = pd.read_csv(res / "Comparison" / "method_availability.csv", index_col=0)
    assert avail.loc["CROSSLAND", "materials"] == 2  # Good and NoK (no k_FI needed)
    assert avail.loc["FIN", "materials"] == 1
    detail = json.loads((res / "Comparison" / "eligibility.json").read_text(encoding="utf-8"))
    assert {d["material"] for d in detail} == {"Good", "Notched", "Runouts", "NoK"}
    assert good.exists()


def _upstream_record(path, material, params):
    """Minimal record in the upstream Calibrations/ schema (A. Chiocca, dc4fe6d)."""
    rec = {
        "schema_version": 1,
        "status": "success",
        "calibration_id": "20261016T000000Z_abc",
        "material_id": material,
        "parameters": params,
    }
    path.write_text(json.dumps(rec), encoding="utf-8")
    return path


def _no_plane_constants(path):
    keep = [i for i, h in enumerate(SUMMARY_HEADERS) if not h.startswith(("k FS", "k FI"))]
    return write_workbook(
        path,
        summary_headers=[SUMMARY_HEADERS[i] for i in keep],
        summary_values=[SUMMARY_VALUES[i] for i in keep],
    )


def test_constants_from_upstream_record(tmp_path):
    wb = _no_plane_constants(tmp_path / "Rec.xlsx")
    rec = _upstream_record(
        tmp_path / "r.json",
        "Rec",
        {"k_FS": 2.0, "k_FI": 0.4, "L_LI": 3.0, "k_MGSE": 1.0, "a_MSWT": 0.53},
    )
    out = run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res", constants=rec)
    assert [o.status for o in out if o.model == "FIN"] == ["ok"]
    side = json.loads((tmp_path / "res" / "Rec_calibration.json").read_text(encoding="utf-8"))
    assert side["params"]["k_FI"] == 0.4
    assert side["params"]["a_MSWT"] == pytest.approx(0.53 / (2 - 0.53))
    src = side["constants_source"]
    assert src["calibration_id"] == "20261016T000000Z_abc" and len(src["sha256"]) == 64
    assert "a/(2 - a)" in src["note"]
    assert "FI" not in (side.get("selection") or {})  # nothing recalibrated


def test_constants_from_our_sidecar(tmp_path):
    wb = _no_plane_constants(tmp_path / "Side.xlsx")
    side = tmp_path / "Side_calibration.json"
    side.write_text(
        json.dumps(
            {"params": {"k_FI": 0.35, "k_FS": 1.5}, "provenance": {"workbook": "x/Side.xlsx"}}
        ),
        encoding="utf-8",
    )
    run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res", constants=side)
    out = json.loads((tmp_path / "res" / "Side_calibration.json").read_text(encoding="utf-8"))
    assert out["params"]["k_FI"] == 0.35


def test_record_for_other_material_rejected(tmp_path):
    wb = _no_plane_constants(tmp_path / "Mine.xlsx")
    rec = _upstream_record(tmp_path / "r.json", "Other", {"k_FI": 0.4})
    with pytest.raises(ValueError, match="Other"):
        run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res", constants=rec)


def test_sidecar_has_failure_fingerprint(tmp_path):
    from multifat.imp.experimental_data import load_material_workbook
    from multifat.report.subsets import dataset_fingerprint

    wb = write_workbook(tmp_path / "F.xlsx")
    run_workbook(wb, models=["CROSSLAND"], results_dir=tmp_path / "res", calibrate=False)
    side = json.loads((tmp_path / "res" / "F_calibration.json").read_text(encoding="utf-8"))
    (block,) = load_material_workbook(wb).values()
    assert side["provenance"]["failure_fingerprint"] == dataset_fingerprint(block)


def test_empty_data_mapping_is_used_as_given(tmp_path):
    """An empty dict is data, not a request to reload the workbook (deliberate ``is None`` test)."""
    from multifat.cli.pipeline import _absent_suppliable_constants

    wb = write_workbook(tmp_path / "F.xlsx")
    with pytest.raises(KeyError):
        _absent_suppliable_constants(wb, ["FIN"], {})


def test_calibrate_main_skips_materials_without_failures(tmp_path, capsys):
    from multifat.cli import calibrate

    rows = [list(r) for r in basquin_rows()]
    for r in rows:
        r[1] = "Yes"
    write_workbook(tmp_path / "OnlyRunouts.xlsx", rows=rows)
    calibrate.main(["--data-folder", str(tmp_path)])
    assert "OnlyRunouts: no plain failures, skipped" in capsys.readouterr().out


def test_numbered_workbook_is_calibrated(tmp_path):
    """The material ID drops the ordering number in calibration too (review focus 1)."""
    wb = _no_plane_constants(tmp_path / "6_Pre.xlsx")
    out = run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res")
    assert [o.status for o in out if o.model == "FIN"] == ["ok"]
    side = json.loads((tmp_path / "res" / "Pre_calibration.json").read_text(encoding="utf-8"))
    assert side["material"] == "Pre" and side["params"]["k_FI"] > 0


def test_duplicate_material_ids_rejected_in_main(tmp_path):
    """Review I1: 1_X.xlsx and 2_X.xlsx must not silently overwrite each other's results."""
    from multifat.cli.pipeline import main

    write_workbook(tmp_path / "1_X.xlsx")
    write_workbook(tmp_path / "2_X.xlsx")
    with pytest.raises(SystemExit, match="X"):
        main([str(tmp_path / "1_X.xlsx"), str(tmp_path / "2_X.xlsx"), "--results-dir",
              str(tmp_path / "res"), "--no-compare", "--models", "CROSSLAND"])  # fmt: skip
    assert not (tmp_path / "res" / "CROSSLAND_X.csv").exists()


def test_missing_nu_with_calibration(tmp_path):
    """Review I2: the default (calibrating) pipeline must skip strain models, not crash."""
    i = SUMMARY_HEADERS.index("Poisson's ratio (-)")
    keep = [j for j, h in enumerate(SUMMARY_HEADERS) if not h.startswith(("k FS", "k FI"))]
    vals = list(SUMMARY_VALUES)
    vals[i] = None
    wb = write_workbook(tmp_path / "NoNu.xlsx", summary_headers=[SUMMARY_HEADERS[j] for j in keep],
                        summary_values=[vals[j] for j in keep])  # fmt: skip
    out = {o.model: o for o in run_workbook(wb, models=["FIN", "SWT", "MSWT"],
                                             results_dir=tmp_path / "res")}  # fmt: skip
    assert out["FIN"].status == "ok"
    assert out["SWT"].status == "skipped" and "E and nu" in out["SWT"].reason
    assert out["MSWT"].status == "skipped" and "E and nu" in out["MSWT"].reason


def test_fatigue_limit_status_in_eligibility(tmp_path):
    """Review I3: limit-state sets get their own status and count no failures."""
    import pandas as pd

    from multifat.cli.pipeline import main

    rows = [list(r) for r in basquin_rows()]
    for r in rows:
        r[0], r[1] = 1e6, "Yes"
    write_workbook(tmp_path / "Lim.xlsx", rows=rows)
    res = tmp_path / "res"
    with pytest.raises(SystemExit):
        main([str(tmp_path / "Lim.xlsx"), str(write_workbook(tmp_path / "Good.xlsx")), "--results-dir", str(res), "--no-compare",
              "--models", "CROSSLAND"])  # fmt: skip
    m = pd.read_csv(res / "Comparison" / "eligibility_matrix.csv", index_col=0)
    assert m.loc["CROSSLAND", "Lim"] == "fatigue_index"
    detail = json.loads((res / "Comparison" / "eligibility.json").read_text(encoding="utf-8"))
    (lim,) = [d for d in detail if d["material"] == "Lim" and d["method"] == "CROSSLAND"]
    assert lim["n_cases"] == 0  # limit states are not counted as failures


def test_record_for_other_dataset_is_flagged(tmp_path):
    """Review I4: a record made on different failures is used but flagged in the sidecar."""
    wb = _no_plane_constants(tmp_path / "Rec.xlsx")
    rec = tmp_path / "r.json"
    rec.write_text(json.dumps({"schema_version": 1, "status": "success", "calibration_id": "c1",
                               "material_id": "Rec", "parameters": {"k_FI": 0.4},
                               "counts": {"failure": 70}}), encoding="utf-8")  # fmt: skip
    run_workbook(wb, models=["FIN"], results_dir=tmp_path / "res", constants=rec)
    side = json.loads((tmp_path / "res" / "Rec_calibration.json").read_text(encoding="utf-8"))
    assert "70 failures" in side["constants_source"]["warning"]
    assert "failure_fingerprint" in side["provenance"]
    assert "input_fingerprint" not in side["provenance"]


def test_constants_needs_exactly_one_workbook(tmp_path):
    """Minor: --constants with several workbooks is a usage error, not a traceback."""
    from multifat.cli.pipeline import main

    a, b = write_workbook(tmp_path / "A.xlsx"), write_workbook(tmp_path / "B.xlsx")
    rec = tmp_path / "r.json"
    rec.write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        main([str(a), str(b), "--constants", str(rec), "--no-compare"])
    assert e.value.code == 2  # argparse usage error


def test_material_accepts_full_stem_and_reports_unknown(tmp_path, monkeypatch, capsys):
    """Minor: --material 6_X works; an unknown material is reported without a traceback."""
    from multifat.cli import pipeline, run

    write_workbook(tmp_path / "6_X.xlsx")
    assert run.resolve_workbook("6_X", tmp_path) == tmp_path / "6_X.xlsx"
    monkeypatch.setattr(pipeline, "DATA_FOLDER", str(tmp_path))
    with pytest.raises(SystemExit) as e:
        pipeline.main(["--material", "Nope", "--no-compare", "--results-dir", str(tmp_path / "r")])
    assert e.value.code == 2 and "Nope" in capsys.readouterr().err


def test_single_workbook_keeps_existing_eligibility(tmp_path):
    """Minor: a one-material run does not overwrite the matrix of an earlier --all."""
    from multifat.cli.pipeline import main

    wb = write_workbook(tmp_path / "A.xlsx")
    comp = tmp_path / "res" / "Comparison"
    comp.mkdir(parents=True)
    (comp / "eligibility_matrix.csv").write_text("earlier\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        main([str(wb), "--results-dir", str(tmp_path / "res"), "--no-compare", "--no-calibrate",
              "--models", "CROSSLAND"])  # fmt: skip
    assert (comp / "eligibility_matrix.csv").read_text(encoding="utf-8") == "earlier\n"


def test_run_workbook_loads_each_workbook_twice_at_most(tmp_path, monkeypatch):
    """Minor: validation plus one import, not four imports per workbook."""
    from multifat.imp import experimental_data as exp

    calls = []
    real = exp.load_material_workbook
    monkeypatch.setattr(
        exp, "load_material_workbook", lambda *a, **k: calls.append(1) or real(*a, **k)
    )
    wb = write_workbook(tmp_path / "A.xlsx")
    run_workbook(wb, models=["CROSSLAND"], results_dir=tmp_path / "res", calibrate=False)
    assert len(calls) <= 2 + 1  # validate, load once, and the one model job (sequential run)
