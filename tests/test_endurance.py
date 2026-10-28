"""Fatigue index error on fatigue-limit datasets (decision d057, Boehme & Papuga 2024 Eqs. 8-9)."""

import numpy as np
import pandas as pd
import pytest

from multifat.eval.endurance import (
    calibrate_endurance_constant,
    case_from_amplitudes,
    fatigue_index,
    is_fatigue_limit_dataset,
    reference_cases,
)
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list
from workbook_factory import basquin_rows, write_workbook

R_LIST = build_R_list(n1=10, n2=20)
PARAMS = {"Sigm1": 420.0, "Taum1": 300.0, "E": 206000.0, "nu": 0.3, "Su": 900.0, "Sy": 750.0}


def _refs():
    return reference_cases(PARAMS["Sigm1"], PARAMS["Taum1"], PARAMS["E"], PARAMS["nu"])


def test_fi_zero_at_reference():
    refs = _refs()
    (cross,) = get_models(["CROSSLAND"])
    (row,) = fatigue_index(cross, PARAMS, [refs["axial"]], refs, R_LIST)
    assert row["FI_axial"] == 0.0
    (fin,) = get_models(["FIN"])
    k = calibrate_endurance_constant(fin, "k_FI", PARAMS, refs, (0.0, 3.0), R_LIST)
    (row,) = fatigue_index(fin, {**PARAMS, "k_FI": k}, [refs["axial"]], refs, R_LIST)
    assert row["FI_axial"] == 0.0


def test_fi_linear_criterion_scaling():
    refs = _refs()
    (cross,) = get_models(["CROSSLAND"])
    case = case_from_amplitudes(axial=1.1 * PARAMS["Sigm1"], E=PARAMS["E"], nu=PARAMS["nu"])
    (row,) = fatigue_index(cross, PARAMS, [case], refs, R_LIST)
    assert row["FI_axial"] == pytest.approx(0.1, rel=1e-9)


def test_endurance_constant_equalises_refs():
    refs = _refs()
    (fin,) = get_models(["FIN"])
    k = calibrate_endurance_constant(fin, "k_FI", PARAMS, refs, (0.0, 3.0), R_LIST)
    assert k is not None and 0.0 < k < 3.0
    rows = fatigue_index(fin, {**PARAMS, "k_FI": k}, [refs["axial"], refs["torsion"]], refs, R_LIST)
    dp_ax, dp_tor = rows[0]["DP"], rows[1]["DP"]
    assert abs(dp_ax - dp_tor) / dp_ax < 1e-8


def _limit_rows():
    rows = [list(r) for r in basquin_rows()]
    for r in rows:
        r[0], r[1] = 1e6, "Yes"  # every plain row a fatigue-limit state at 1e6 cycles
    return rows


def test_pipeline_routes_fatigue_limit_workbook(tmp_path):
    from multifat.cli.pipeline import run_workbook
    from multifat.cli.run import validate_workbook
    from multifat.imp.experimental_data import load_material_workbook

    path = write_workbook(tmp_path / "Lim.xlsx", rows=_limit_rows())
    (block,) = load_material_workbook(path).values()
    assert is_fatigue_limit_dataset(block)
    rep = validate_workbook(path)
    assert rep.ok and any("fatigue-limit dataset" in s for s in rep.info)
    res = tmp_path / "res"
    out = {
        o.model: o for o in run_workbook(path, models=["CROSSLAND", "FIN", "MDC"], results_dir=res)
    }
    assert out["CROSSLAND"].status == "ok" and out["FIN"].status == "ok"
    assert out["MDC"].status == "skipped"
    df = pd.read_csv(res / "FatigueLimit" / "CROSSLAND_Lim_FI.csv")
    assert len(df) == len(block["runouts"]) and np.isfinite(df.FI_axial).all()
    assert not list(res.glob("CROSSLAND_Lim.csv"))  # no master-curve results
    summary = pd.read_csv(res / "FatigueLimit" / "Lim_fatigue_index_summary.csv")
    assert set(summary.Method) == {"CROSSLAND", "FIN"}


def test_fatigue_limit_pipeline_honours_explicit_ext(tmp_path):
    from multifat.cli.pipeline import run_workbook

    path = write_workbook(tmp_path / "LimitVariant.xlsx", rows=_limit_rows())
    outdir = tmp_path / "res"
    outcomes = run_workbook(path, models=["FIN_ext"], results_dir=outdir)
    assert [(row.model, row.status) for row in outcomes] == [("FIN_ext", "ok")]
    assert (outdir / "FatigueLimit" / "FIN_ext_LimitVariant_FI.csv").exists()
    assert not (outdir / "FatigueLimit" / "FIN_LimitVariant_FI.csv").exists()
    summary = pd.read_csv(outdir / "FatigueLimit" / "LimitVariant_fatigue_index_summary.csv")
    assert list(summary.Method) == ["FIN_ext"]


def test_runout_workbook_without_limits_rejected(tmp_path):
    from multifat.cli.run import validate_workbook
    from workbook_factory import SUMMARY_HEADERS, SUMMARY_VALUES

    i = SUMMARY_HEADERS.index("Tensile fatigue limit* (MPa)")
    vals = SUMMARY_VALUES[:i] + [None] + SUMMARY_VALUES[i + 1 :]
    path = write_workbook(tmp_path / "NoLim.xlsx", rows=_limit_rows(), summary_values=vals)
    rep = validate_workbook(path)
    assert not rep.ok and any("fatigue limit" in e for e in rep.errors)


def test_compare_ignores_fatigue_limit_outputs(tmp_path):
    from multifat.cli.compare import collect_all_results

    (tmp_path / "FatigueLimit").mkdir()
    pd.DataFrame({"Material": ["X"], "Method": ["A"], "FI_axial": [0.1]}).to_csv(
        tmp_path / "FatigueLimit" / "A_X_FI.csv", index=False
    )
    pd.DataFrame({"Material": ["X"], "Method": ["A"], "Nf_exp": [1e5], "CP_value": [1.0]}).to_csv(
        tmp_path / "A_X.csv", index=False
    )
    df = collect_all_results(tmp_path)
    assert "FI_axial" not in df.columns and len(df) == 1


def test_summary_counts_nonfinite_separately():
    """MWCM is infinite past its rho cut-off: such rows are counted, not averaged."""
    from multifat.eval.endurance import fatigue_index_summary

    df = pd.DataFrame({"Method": ["A", "A", "A"], "FI_axial": [0.1, -0.1, np.inf],
                       "FI_torsion": [0.2, 0.0, np.inf]})  # fmt: skip
    (row,) = fatigue_index_summary(df, "X").to_dict("records")
    assert row["n"] == 3 and row["n_nonfinite"] == 1
    assert row["FI_axial_mean"] == pytest.approx(0.0) and row["FI_axial_max"] == pytest.approx(0.1)


def test_fi_load_equals_fi_for_linear_criterion():
    """d063: FI_load = 1/s - 1 with DP(s * case) = DP_ref; Crossland is linear in the loads."""
    refs = _refs()
    (cross,) = get_models(["CROSSLAND"])
    case = case_from_amplitudes(axial=1.1 * PARAMS["Sigm1"], E=PARAMS["E"], nu=PARAMS["nu"])
    (row,) = fatigue_index(cross, PARAMS, [case], refs, R_LIST)
    assert row["FI_load"] == pytest.approx(0.1, rel=1e-8)
    (ref_row,) = fatigue_index(cross, PARAMS, [refs["axial"]], refs, R_LIST)
    assert abs(ref_row["FI_load"]) < 1e-10


def test_fi_load_is_stress_scaled_for_swt():
    """SWT ~ sigma^2: FI_axial = 1.1^2 - 1 but FI_load = 0.1."""
    refs = _refs()
    (swt,) = get_models(["SWT"])
    case = case_from_amplitudes(axial=1.1 * PARAMS["Sigm1"], E=PARAMS["E"], nu=PARAMS["nu"])
    (row,) = fatigue_index(swt, PARAMS, [case], refs, R_LIST)
    assert row["FI_axial"] == pytest.approx(0.21, rel=1e-6)
    assert row["FI_load"] == pytest.approx(0.1, rel=1e-6)
