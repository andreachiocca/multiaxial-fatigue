"""Console entry point for cross-material criteria ranking (issue #57, decision d116).

Usage::

    uv run multifat-rank [RESULTS_DIR] [--out-dir OUT_DIR] [--min-tests 3] [--min-coverage 4]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from multifat.cli.compare import collect_all_results
from multifat.report.ranking import evaluate_ranking


def build_parser() -> argparse.ArgumentParser:
    """Build command line argument parser for multifat-rank."""
    parser = argparse.ArgumentParser(
        prog="multifat-rank",
        description="Cross-material relative RMS load error ranking and complexity KPI tool.",
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        type=Path,
        default=Path("Results"),
        help="Directory containing per-test result CSVs (default: Results)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory to save ranking CSV outputs (default: RESULTS_DIR/Ranking)",
    )
    parser.add_argument(
        "--min-tests",
        type=int,
        default=3,
        help="Minimum assessment failures per (method, material) cell (default: 3)",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=4,
        help="Minimum material coverage K for family selection (default: 4)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for bootstrap standard error estimation (default: 0)",
    )
    parser.add_argument(
        "--n-boot",
        type=int,
        default=1000,
        help="Number of bootstrap iterations (default: 1000)",
    )
    return parser


def format_ranking_summary(res) -> str:
    """Format ranking result into a human-readable text report."""
    lines: list[str] = []
    lines.append("=" * 80)
    lines.append("MULTIFAT CRITERIA RANKING & COMPLEXITY REPORT")
    lines.append("=" * 80)

    # 1. Method ranking table
    lines.append("\n1. Overall Method Ranking (sorted by Family and Accuracy KPI Z):")
    lines.append("-" * 80)
    display_cols = [
        "Method",
        "Family",
        "Z",
        "SE",
        "Coverage",
        "w10_med",
        "noncons_med",
        "p_in",
        "p_fit",
    ]
    formatted = res.method_ranking[display_cols].copy()
    formatted["Z"] = formatted["Z"].map(lambda x: f"{x:+.3f}")
    formatted["SE"] = formatted["SE"].map(lambda x: f"±{x:.3f}")
    formatted["w10_med"] = formatted["w10_med"].map(lambda x: f"{x * 100:.1f}%")
    formatted["noncons_med"] = formatted["noncons_med"].map(lambda x: f"{x * 100:.1f}%")
    lines.append(formatted.to_string(index=False))

    # 2. Equivalence classes
    if res.equivalence_classes:
        lines.append("\n2. Equivalence Classes (identical assessment errors across shared data):")
        lines.append("-" * 80)
        for grp in res.equivalence_classes:
            lines.append(f"  - {' ≡ '.join(grp)}")

    # 3. Family-wise 1-SE selection
    lines.append(
        "\n3. Family-wise 1-SE Selection (simplest among statistically indistinguishable):"
    )
    lines.append("-" * 80)
    if not res.family_selection.empty:
        sel_cols = ["Family", "Method", "Z", "p_in", "p_fit", "Candidate_1SE", "Selected"]
        sel_fmt = res.family_selection[sel_cols].copy()
        sel_fmt["Z"] = sel_fmt["Z"].map(lambda x: f"{x:+.3f}")
        lines.append(sel_fmt.to_string(index=False))

        lines.append("\nPicks per family:")
        for fam, g in res.family_selection[res.family_selection["Selected"]].groupby("Family"):
            p = g.iloc[0]
            lines.append(
                f"  - {fam:15s}: {p['Method']} (Z = {p['Z']:.3f}, p_in = {p['p_in']}, p_fit = {p['p_fit']})"
            )

    # 4. Leave-one-material-out stability
    if res.lomo_stability:
        lines.append("\n4. Leave-One-Material-Out Stability of Family Picks:")
        lines.append("-" * 80)
        for fam, counts in res.lomo_stability.items():
            counts_str = ", ".join(
                f"{m}: {c}" for m, c in sorted(counts.items(), key=lambda x: -x[1])
            )
            lines.append(f"  - {fam:15s}: {counts_str}")

    lines.append("=" * 80)
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> None:
    """CLI entry point for multifat-rank."""
    parser = build_parser()
    args = parser.parse_args(argv)

    results_dir = args.results_dir
    if not results_dir.exists():
        sys.stderr.write(f"Error: Results directory '{results_dir}' does not exist.\n")
        sys.exit(1)

    out_dir = args.out_dir or (results_dir / "Ranking")
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        df = collect_all_results(results_dir)
    except Exception as e:
        sys.stderr.write(f"Error reading result CSVs from {results_dir}: {e}\n")
        sys.exit(1)

    if "Error_load" not in df.columns:
        sys.stderr.write(
            f"Error: Results in {results_dir} do not contain 'Error_load'. "
            "Ensure multifat-run or multifat-pipeline evaluated load-equivalent errors.\n"
        )
        sys.exit(1)

    res = evaluate_ranking(
        df,
        min_tests=args.min_tests,
        min_coverage=args.min_coverage,
        seed=args.seed,
        n_boot=args.n_boot,
    )

    if res.method_ranking.empty:
        sys.stderr.write(f"No qualifying methods met the criteria (min_tests={args.min_tests}).\n")
        sys.exit(1)

    # Export CSVs
    summary_path = out_dir / "summary_ranking.csv"
    res.method_ranking.to_csv(summary_path, index=False)

    selection_path = out_dir / "family_selection.csv"
    res.family_selection.to_csv(selection_path, index=False)

    cells_path = out_dir / "ranking_cell_stats.csv"
    res.cell_stats.drop(columns=["_errors"]).to_csv(cells_path, index=False)

    # Print report
    print(format_ranking_summary(res))
    print(f"\n[INFO] Ranking tables saved to: {out_dir}")


if __name__ == "__main__":
    main()
