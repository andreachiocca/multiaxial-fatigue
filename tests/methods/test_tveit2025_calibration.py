"""Calibration of TVEIT2025: Tveit et al. (2025) Sect. 4, Eqs. (22)-(26); Ottosen et al. (2008) Sect. 5.

Paper numbers: Ottosen, Stenstroem & Ristinmaa (2008), Int. J. Fatigue 30, 996-1006, Eq. (31): the SAE
4340 fatigue limits sigma_a = 490 MPa at sigma_m = 0 and sigma_a = sigma_m = 400 MPa give A = 0.225.

Cross-implementation check (computed values, not paper numbers): the closed form of Eqs. (23), (24) and
(29) with the SAE 4340 constants of Eqs. (31)-(32) against the reference table of an independent
implementation in the FABER WG4.3 model appendix (alpha_2, alpha_4 to 0.01 MPa, N to its three printed
digits). Model consistency: the closed form against the model's own time integration.

Reproducibility: with the stage-1 inputs, the stage-1 hand-entered constants
(tests/data/tveit2025_stage1_calibrations.csv) reproduce their Eq. (25) value, the three-point rule
picks their three tests from the workbooks, and the profile of Eq. (25) over C reaches their E_min.
C, K and L themselves are not reproducible: uniaxial data do not identify C, so equivalent sets span
three to four decades of C.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multifat.eval import tveit2025_calibration as tc
from multifat.imp import experimental_data as exp
from multifat.models.invariant_methods.tveit2025 import tveit2025_life

ROOT = Path(__file__).resolve().parents[2]
STAGE1 = pd.read_csv(ROOT / "tests" / "data" / "tveit2025_stage1_calibrations.csv", comment="#")

# Ottosen et al. (2008) Eqs. (31)-(32), SAE 4340
SAE4340 = dict(C=1.25, K=2.65e-5, L=14.4, A=0.225, f_m1=490.0)
# FABER WG4.3 model appendix, reference results for code verification (computed, closed form)
APPENDIX = [
    # sigma_m, sigma_a, alpha_2, alpha_4, N
    (0.0, 670.0, 133.51, 38.07, 9.54e3),
    (0.0, 600.0, 109.02, 52.05, 3.96e4),
    (0.0, 560.0, 95.65, 59.92, 1.01e5),
    (0.0, 510.0, 79.66, 69.64, 6.06e5),
    (100.0, 650.0, 196.38, 101.88, 8.59e3),
    (200.0, 620.0, 255.72, 167.65, 9.41e3),
    (400.0, 560.0, 374.63, 299.05, 1.14e4),
]


def _points(row):
    out = []
    for item in row.points.split("|"):
        specimen, sigma_a, N = item.split(":")
        out.append((specimen, float(sigma_a), float(N)))
    return out


def _axial(sigma_m, sigma_a):
    S_mean, S_sin = np.zeros((3, 3)), np.zeros((3, 3))
    S_mean[1, 1], S_sin[1, 1] = sigma_m, sigma_a
    return S_mean, S_sin, np.zeros((3, 3))


def _test(row, sigma_a, N, sigma_m=0.0, source=("a", "doi-a")):
    return dict(
        case_id=row,
        specimen_id=f"S{row}",
        excel_row=row,
        sigma_m=sigma_m,
        sigma_a=sigma_a,
        N_exp=N,
        source=list(source),
    )


# Routes to A (Tveit et al. 2025, Sect. 4)


def test_eq31_mean_stress_slope_from_two_fatigue_limits():
    spec = {"mean_stress_fatigue_limit": {"sigma_a": 400, "sigma_m": 400}}
    slope = tc.mean_stress_slope(spec, 490.0, {})
    assert slope["A"] == pytest.approx(0.225, abs=1e-12)
    assert slope["route"] == "mean_stress_fatigue_limit"


@pytest.mark.parametrize(
    "inputs, f_m1, Su, route, expected",
    [
        ({"mean_stress_fatigue_limit": {"sigma_a": 100.0, "R": 0.1}}, 150.0, 595, "mean_stress_fatigue_limit", 9.0 / 22.0),
        ({"A": {"value": 0.3, "source": "test"}}, 500.0, 980, "supplied", 0.3),
        ({"material_type": "steel"}, 433.0, 875.0, "fkm", 0.20625),
        ({"material_type": "titanium alloy"}, 500.0, 980.0, "goodman", 500.0 / 980.0),
        ({}, 500.0, 980.0, "goodman", 500.0 / 980.0),
    ],
    ids=["R=0.1 limit", "hand-entered", "FKM steel", "Goodman, no FKM values", "Goodman, no type"],
)  # fmt: skip
def test_routes_to_A_in_the_order_of_the_authors(inputs, f_m1, Su, route, expected):
    slope = tc.mean_stress_slope(inputs, f_m1, {"Tensile_strength_MPa": Su})
    assert slope["route"] == route and slope["A"] == pytest.approx(expected)


def test_routes_to_A_reject_what_the_closed_form_cannot_use():
    with pytest.raises(tc.CalibrationInputError, match="outside"):
        tc.mean_stress_slope({}, 500.0, {"Tensile_strength_MPa": 400.0})
    with pytest.raises(tc.CalibrationInputError, match="tensile mean"):
        spec = {"mean_stress_fatigue_limit": {"sigma_a": 100.0, "sigma_m": -20.0}}
        tc.mean_stress_slope(spec, 150.0, {})
    with pytest.raises(tc.CalibrationInputError, match="not both"):
        spec = {"mean_stress_fatigue_limit": {"sigma_a": 100.0, "R": 0.1}, "A": {"value": 0.3}}
        tc.mean_stress_slope(spec, 150.0, {})
    with pytest.raises(tc.CalibrationInputError, match="Su"):
        tc.mean_stress_slope({}, 150.0, {"Tensile_strength_MPa": None})


def test_inputs_file_is_optional_and_checked(tmp_path):
    assert tc.read_inputs(tmp_path / "absent.json", "M") == ({}, None)
    path = tmp_path / "in.json"
    path.write_text(json.dumps({"schema": tc.INPUTS_SCHEMA, "material_id": "M", "Sigm1": 1}))
    with pytest.raises(tc.CalibrationInputError, match="unknown keys Sigm1"):
        tc.read_inputs(path, "M")
    path.write_text(json.dumps({"schema": tc.INPUTS_SCHEMA, "material_id": "N"}))
    with pytest.raises(tc.CalibrationInputError, match="not 'M'"):
        tc.read_inputs(path, "M")


def test_committed_inputs_files_are_valid():
    for path in sorted((ROOT / "Calibrations").glob(f"*/{tc.INPUTS_FILE}")):
        inputs, digest = tc.read_inputs(path, path.parent.name)
        assert digest and inputs["material_id"] == path.parent.name


# Closed form (Ottosen et al. 2008)


@pytest.mark.parametrize("row", APPENDIX, ids=lambda r: f"sm{r[0]:g}-sa{r[1]:g}")
def test_closed_form_matches_independent_implementation(row):
    sigma_m, sigma_a, alpha_2, alpha_4, N = row
    p = SAE4340
    life = tc.uniaxial_life(
        p["C"], p["K"], p["L"], sigma_m, sigma_a, p["A"], p["f_m1"], check_unique=True
    )
    assert life.status == tc.VALID
    assert life.alpha_2 == pytest.approx(alpha_2, abs=0.006)
    assert life.alpha_4 == pytest.approx(alpha_4, abs=0.006)
    assert float(f"{life.Nf:.3g}") == pytest.approx(N)


@pytest.mark.parametrize("sigma_m, sigma_a", [(0.0, 600.0), (200.0, 620.0)])
def test_closed_form_agrees_with_the_time_integration(sigma_m, sigma_a):
    """Ottosen et al. (2008): the periodic state is reached after a few cycles."""
    p = SAE4340
    closed = tc.uniaxial_life(p["C"], p["K"], p["L"], sigma_m, sigma_a, p["A"], p["f_m1"]).Nf
    t_m1 = p["f_m1"] / math.sqrt(3.0)  # irrelevant for a uniaxial cycle
    integrated = tveit2025_life(
        *_axial(sigma_m, sigma_a), p["f_m1"], t_m1, p["A"], p["C"], p["K"], p["L"], max_blocks=300
    ).Nf  # extrapolated from the damage of block 300, long after the transient
    assert abs(math.log10(integrated / closed)) < 0.03


def test_endurance_limit_of_the_closed_form():
    p = SAE4340
    for sigma_m in (0.0, 200.0):
        sigma_a_lim = p["f_m1"] - p["A"] * sigma_m  # Haigh line
        below = tc.uniaxial_life(
            p["C"], p["K"], p["L"], sigma_m, 0.999 * sigma_a_lim, p["A"], p["f_m1"]
        )
        above = tc.uniaxial_life(
            p["C"], p["K"], p["L"], sigma_m, 1.01 * sigma_a_lim, p["A"], p["f_m1"]
        )
        assert below.status == tc.BELOW and math.isinf(below.Nf)
        assert above.status == tc.VALID and 1e5 < above.Nf < math.inf


# Calibration tests: one source, three points


def test_three_point_rule():
    pool = [
        _test(2, 300.0, 1e7),  # at the limit: not above it
        _test(3, 320.0, 2e6),
        _test(4, 320.0, 8e5),
        _test(5, 320.0, 1e6),  # median replicate of LOW
        _test(6, 400.0, 1e5),
        _test(7, 400.0, 3e5),  # MID: two replicates equally far from their median -> lower row 6
        _test(8, 450.0, 9e4),
        _test(9, 500.0, 1e4),  # HIGH; sqrt(320 * 500) = 400 is the MID target
    ]
    chosen = tc.three_point_rule(pool, 300.0)
    assert [(c["role"], c["excel_row"]) for c in chosen] == [("LOW", 5), ("MID", 6), ("HIGH", 9)]
    assert [c["replicates"] for c in chosen] == [3, 2, 1]
    with pytest.raises(tc.CalibrationInputError, match="three stress levels"):
        tc.three_point_rule(pool[:6], 300.0)


def test_calibration_tests_come_from_the_source_of_the_limits():
    tests = [
        _test(2, 450.0, 1e5, source=("santus et al.", "doi-s")),
        _test(3, 500.0, 1e5, source=("benedetti et al.", "doi-b")),
        _test(4, 520.0, 1e5, source=("anes et al.", "doi-a2")),
    ]

    def limits(authors, reference):
        return {"sources": [{"authors": authors, "reference": reference}]}

    def rows(pool):
        return [t["excel_row"] for t in pool]

    pool, how = tc.calibration_source(tests, limits("Benedetti et al.", "DOI-B"), None)
    assert rows(pool) == [3] and how["match"] == "same reference"
    pool, how = tc.calibration_source(tests, limits("Anes et al.", "doi-a1"), None)
    assert rows(pool) == [4] and how["match"] == "same authors"
    pool, how = tc.calibration_source(tests, limits("Other", "doi-x"), None)
    assert rows(pool) == [2, 3, 4] and how["match"] == "no match"
    pool, how = tc.calibration_source(tests, {"sources": []}, None)
    assert rows(pool) == [2, 3, 4] and how["match"] == "unknown"
    pool, how = tc.calibration_source(tests, limits("Other", "doi-x"), "Santus")
    assert rows(pool) == [2] and how["match"] == "declared"


# Confidence of the calibration


def _codes(findings):
    return {f.code: f.severity for f in findings}


def test_input_findings():
    failures = [_test(2, 380.0, 1e6), _test(3, 450.0, 1e5)]
    runouts = [_test(4, 440.0, 1e7)]
    slope = tc.mean_stress_slope({}, 433.0, {"Tensile_strength_MPa": 875.0})
    found = tc._input_findings(
        433.0, failures, runouts, {"match": "no match", "source": "x"}, slope
    )
    assert _codes(found) == {
        "failures_below_limit": "major",
        "runouts_above_limit": "minor",
        "limit_source_differs": "major",
        "A_estimated": "minor",  # the authors' own route without a mean-stress limit
        "material_type_for_FKM": "minor",
    }
    slope = tc.mean_stress_slope({"mean_stress_fatigue_limit": {"sigma_a": 100, "R": 0.1}}, 150, {})
    assert tc._input_findings(150.0, [], [], {"match": "same reference", "source": ""}, slope) == []


def test_point_findings():
    close = [
        {**_test(2, 523.2, 9e4), "role": "LOW", "replicates": 1},
        {**_test(3, 534.7, 7e4), "role": "MID", "replicates": 2},
        {**_test(4, 550.5, 4e4, source=("b", "doi-b")), "role": "HIGH", "replicates": 1},
    ]
    assert _codes(tc._point_findings(close, 500.0)) == {
        "narrow_stress_range": "major",
        "low_point_near_limit": "minor",
        "single_specimen_level": "minor",
        "points_mixed_sources": "major",
    }
    wide = [{**p, "sigma_a": s, "replicates": 3} for p, s in zip(close, (600.0, 700.0, 800.0))]
    wide[2]["source"] = wide[0]["source"]
    assert tc._point_findings(wide, 500.0) == []


def test_fit_findings():
    point = {"residual_closed": 0.35, "residual_transient": 0.1}
    chosen = {"C": 1e-4, "K": 1e-5, "L": math.exp(math.log(0.05)), "points": [point]}
    # C free over decades: its bound is the flat E falling on a grid end, not a finding of its own
    found = tc._fit_findings(chosen, {"C_range": [1e-4, 3.16]})
    assert _codes(found) == {
        "poor_fit": "minor",
        "constant_at_bound": "major",
        "C_not_identified": "minor",
    }
    bound = next(f.what for f in found if f.code == "constant_at_bound")
    assert bound.startswith("L of the chosen set")
    # C pinned by the equivalent sets: on its bound it is a finding
    found = tc._fit_findings(chosen, {"C_range": [1e-4, 5e-4]})
    assert _codes(found) == {"poor_fit": "minor", "constant_at_bound": "major"}
    assert "C, L of the chosen set" in next(f.what for f in found if f.code == "constant_at_bound")


# Reuse of an identical calibration


def test_an_archived_calibration_of_the_same_identity_is_reused(tmp_path):
    identity = {"procedure": tc.PROCEDURE, "workbook_sha256": "w", "inputs_sha256": None}
    entry = {"status": "calibrated", "identity": identity, "params": {"C_TVEIT": 0.1}}
    record = {"calibration_id": "id1", "model_calibrations": {"TVEIT2025": entry}}
    (tmp_path / "id1.json").write_text(json.dumps(record))
    (tmp_path / "TVEIT2025_inputs.json").write_text(json.dumps({"schema": tc.INPUTS_SCHEMA}))
    reused = tc.reusable_calibration([tmp_path / "absent", tmp_path], identity)
    assert reused["params"] == {"C_TVEIT": 0.1} and reused["reused_from"]["calibration_id"] == "id1"
    assert tc.reusable_calibration([tmp_path], {**identity, "inputs_sha256": "other"}) is None


# Reproducibility of stage 1 (stage-1 inputs)


@pytest.mark.parametrize("row", list(STAGE1.itertuples()), ids=lambda r: r.material)
def test_stage1_constants_reproduce_their_objective(row):
    points = [(0.0, sigma_a, N) for _, sigma_a, N in _points(row)]
    E = tc.eq25_objective(row.C, row.K, row.L, points, row.A, row.f_m1)
    assert E == pytest.approx(row.E_stage1, rel=1e-4)
    assert E <= tc.EQUIVALENT * row.E_min_stage1  # the stage-1 set is in the equivalent family


@pytest.mark.parametrize("row", list(STAGE1.itertuples()), ids=lambda r: r.material)
def test_three_point_rule_on_the_workbooks_and_profile_reaches_stage1_minimum(row):
    workbook = ROOT / "Experimental_Data" / row.workbook
    block = exp.load_material_workbook(workbook)[row.material]
    failures = tc.axial_tests(block["cases"], tc.dataset_sources(workbook), failures=True)
    pool, _ = tc.calibration_source(failures, {"sources": []}, row.source)
    chosen = tc.three_point_rule(pool, row.f_m1)
    expected = _points(row)
    assert [c["specimen_id"] for c in chosen] == [s for s, _, _ in expected]
    assert [c["sigma_a"] for c in chosen] == pytest.approx([a for _, a, _ in expected])
    assert [c["N_exp"] for c in chosen] == pytest.approx([n for _, _, n in expected])
    triples = [(c["sigma_m"], c["sigma_a"], c["N_exp"]) for c in chosen]
    fit = tc.fit_ckl(triples, row.A, row.f_m1)
    assert fit.E_min <= row.E_min_stage1 * (1 + 1e-3)  # the profile reaches the stage-1 minimum
    assert min(p["E"] for p in fit.profile if p["C"] <= 0.01) > 0.0  # not an exact fit
    # C is not identified: equivalent sets span more than two decades of C
    assert max(p["C"] for p in fit.family) / min(p["C"] for p in fit.family) > 100.0
