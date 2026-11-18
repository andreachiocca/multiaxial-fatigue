"""Tests for relative RMS load error ranking and complexity KPI (issue #57, decision d116)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from multifat.cli.rank import format_ranking_summary, main as rank_main
from multifat.models.complexity import get_model_complexity
from multifat.models.registry import complexity_table
from multifat.report.ranking import (
    compute_cell_stats,
    evaluate_ranking,
    find_equivalence_classes,
    select_family_candidates,
)


def test_model_complexity_definitions():
    """Verify p_in and p_fit match the review article complexity proposal."""
    expected = {
        # Critical plane
        "MSWT": (1, 3),
        "FIN": (0, 3),
        "FS": (2, 3),
        "SWT": (1, 2),
        "MGSE_ZHU": (1, 3),
        "MDC": (1, 2),
        "MMCCM": (3, 4),
        # Invariant
        "MMMH": (2, 2),
        "GAM": (2, 2),
        "VHN": (3, 2),
        "TVEIT2025": (0, 4),
        # Integral
        "CAIM": (2, 2),
        "PCR": (2, 2),
        "BP": (2, 2),
        "MOREL": (2, 3),
        # Energy
        "EK": (1, 2),
        "SWTD": (1, 2),
    }
    for model_name, (exp_pin, exp_pfit) in expected.items():
        comp = get_model_complexity(model_name)
        assert comp.p_in == exp_pin, f"{model_name}: expected p_in={exp_pin}, got {comp.p_in}"
        assert comp.p_fit == exp_pfit, f"{model_name}: expected p_fit={exp_pfit}, got {comp.p_fit}"


def test_complexity_table():
    """Verify registry.complexity_table produces all registered models."""
    tbl = complexity_table()
    assert len(tbl) >= 50
    methods = {row["Method"] for row in tbl}
    assert "FIN" in methods
    assert "MSWT" in methods
    assert "CAIM" in methods


def test_compute_cell_stats():
    """Verify cell-level bias, sd, rms, w10, noncons calculations."""
    errs_a = np.array([0.05, -0.05, 0.10, -0.15])
    df = pd.DataFrame(
        {
            "Material": ["Mat1"] * 4 + ["Mat1"] * 2,
            "Method": ["M1"] * 4 + ["M1"] * 2,
            "Flag": [1, 1, 1, 1, 0, 0],
            "Error_load": list(errs_a) + [0.01, -0.02],
        }
    )
    cell = compute_cell_stats(df, min_tests=3)
    assert len(cell) == 1
    row = cell.iloc[0]
    assert row["n"] == 4
    assert np.isclose(row["bias"], np.mean(errs_a))
    assert np.isclose(row["sd"], np.std(errs_a, ddof=1))
    assert np.isclose(row["rms"], math.sqrt(np.mean(errs_a**2)))
    # exp(0.05)-1 = 0.051 <= 0.1, exp(-0.05)-1 = -0.048 <= 0.1, exp(0.1)-1 = 0.105 > 0.1, exp(-0.15)-1 = -0.139 > 0.1
    # So 2 of 4 within 10%
    assert np.isclose(row["w10"], 0.5)
    # noncons: e < -ln(1.1) approx -0.0953. Here -0.15 < -0.0953, so 1 of 4.
    assert np.isclose(row["noncons"], 0.25)
    # design rms
    assert np.isclose(row["rms_design"], math.sqrt(np.mean([0.01**2, (-0.02) ** 2])))


def test_equivalence_classes():
    """Verify detection of methods with identical assessment errors."""
    df = pd.DataFrame(
        {
            "Material": ["Mat1"] * 4 + ["Mat1"] * 4,
            "Method": ["M1"] * 4 + ["M2"] * 4,
            "Flag": [1] * 8,
            "Error_load": [0.05, -0.02, 0.08, -0.04] * 2,
        }
    )
    cell = compute_cell_stats(df, min_tests=3)
    equiv = find_equivalence_classes(cell)
    assert len(equiv) == 1
    assert sorted(equiv[0]) == ["M1", "M2"]


def test_family_selection_1se_and_complexity():
    """Verify 1-SE candidate filter and simplest-model selection."""
    ranking = pd.DataFrame(
        [
            {
                "Method": "COMPLEX_BEST",
                "Family": "Critical plane",
                "Z": -0.40,
                "SE": 0.05,
                "K": 5,
                "p_in": 3,
                "p_fit": 4,
            },
            {
                "Method": "SIMPLE_CANDIDATE",
                "Family": "Critical plane",
                "Z": -0.37,  # within -0.40 + 0.05 = -0.35
                "SE": 0.04,
                "K": 5,
                "p_in": 0,
                "p_fit": 2,
            },
            {
                "Method": "POOR_MODEL",
                "Family": "Critical plane",
                "Z": -0.20,  # exceeds -0.35
                "SE": 0.04,
                "K": 5,
                "p_in": 0,
                "p_fit": 2,
            },
        ]
    )
    sel = select_family_candidates(ranking, min_coverage=4)
    assert len(sel) == 3
    cand_map = dict(zip(sel["Method"], sel["Candidate_1SE"]))
    picked_map = dict(zip(sel["Method"], sel["Selected"]))

    assert cand_map["COMPLEX_BEST"] is True
    assert cand_map["SIMPLE_CANDIDATE"] is True
    assert cand_map["POOR_MODEL"] is False

    # SIMPLE_CANDIDATE has p_in=0, p_fit=2 vs 3, 4, so it is selected
    assert picked_map["SIMPLE_CANDIDATE"] is True
    assert picked_map["COMPLEX_BEST"] is False


def test_evaluate_ranking_full(tmp_path):
    """Verify evaluate_ranking and CLI tool execution."""
    materials = ["MatA", "MatB", "MatC", "MatD"]
    methods = ["FIN", "FS", "MSWT", "CAIM", "BP"]
    rows = []
    rng = np.random.default_rng(42)

    from multifat.models.references import IMPLEMENTATION_REVISION

    for mat in materials:
        for meth in methods:
            errs = rng.normal(loc=0.02, scale=0.08, size=6)
            for e in errs:
                rows.append(
                    {
                        "Material": mat,
                        "Method": meth,
                        "Flag": 1,
                        "Error_load": e,
                        "Nf_exp": 1e5,
                        "CP_value": 300.0,
                        "Implementation_revision": IMPLEMENTATION_REVISION,
                    }
                )
            d_errs = rng.normal(loc=0.0, scale=0.04, size=3)
            for e in d_errs:
                rows.append(
                    {
                        "Material": mat,
                        "Method": meth,
                        "Flag": 0,
                        "Error_load": e,
                        "Nf_exp": 1e5,
                        "CP_value": 300.0,
                        "Implementation_revision": IMPLEMENTATION_REVISION,
                    }
                )

    df = pd.DataFrame(rows)
    res = evaluate_ranking(df, min_tests=3, min_coverage=3, n_boot=100)

    assert not res.method_ranking.empty
    assert not res.family_selection.empty
    summary_text = format_ranking_summary(res)
    assert "MULTIFAT CRITERIA RANKING & COMPLEXITY REPORT" in summary_text

    # Test multifat-rank CLI
    res_dir = tmp_path / "results"
    res_dir.mkdir()
    df.to_csv(res_dir / "FIN_MatA.csv", index=False)
    out_dir = tmp_path / "out_ranking"

    rank_main([str(res_dir), "--out-dir", str(out_dir), "--min-tests", "3", "--min-coverage", "1"])
    assert (out_dir / "summary_ranking.csv").exists()
    assert (out_dir / "family_selection.csv").exists()
    assert (out_dir / "ranking_cell_stats.csv").exists()

    # Test multifat-compare with --rank
    from multifat.cli.compare import main as compare_main

    comp_out = tmp_path / "out_compare"
    compare_main(
        [
            "--results_dir",
            str(res_dir),
            "--out_dir",
            str(comp_out),
            "--rank",
            "--rank-min-tests",
            "3",
            "--rank-min-coverage",
            "1",
        ]
    )
    assert (comp_out / "summary_ranking.csv").exists()
    assert (comp_out / "family_selection.csv").exists()
