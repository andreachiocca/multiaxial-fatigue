"""Build examples/data/load_error_ranking.csv from pipeline output (issue #37).

Per method and material: number, mean, standard deviation and root mean square of the load-equivalent error
Error_load = ln(1/s*) on the assessed failures (Flag = 1), and the fraction within +-20 % of
the load. Regenerate after any model or dataset change:

    uv run multifat-pipeline --material 42CrMo4_Santus --material 42CrMo4_Reis --material AISI316L \
        --material Al7075_T6 --material Al7075_T651_Wang --material Ti6Al4V_Meneghetti --material Ti6Al4V_Berto --material 42CrMo4_FABEST \
        --models all \
        --jobs 4 --results-dir <tmp_dir>
    uv run python examples/build_load_error_ranking.py <tmp_dir>

Check whether the committed CSV and manifest match the current models and workbooks:

    uv run python examples/build_load_error_ranking.py --check

Fatigue-limit datasets (decision d110) enter through their FatigueLimit/<MODEL>_<MAT>_FI.csv
files, run in the same results folder:

    uv run multifat-pipeline Experimental_Data/2_CX18.xlsx ... Experimental_Data/17_Ti6Al4V_Delahay.xlsx \
        --models all --jobs 4 --results-dir <tmp_dir>

A limit state has the same load-equivalent error: FI_load = 1/s - 1 is the stress-scaled fatigue
index (multifat.eval.endurance), so Error_load = ln(1/s) = ln(1 + FI_load), exactly and for any
damage parameter. As the R = -1 uniaxial tension failures of a life dataset are its design set,
the pure R = -1 axial and torsion limit states of a limit dataset are its reference and are left
out: they are sigma_-1 and tau_-1, inputs of most limit criteria, which would reproduce them
exactly (decision d118, revising d110). The column Assessment tells the two
kinds apart (fatigue_life, fatigue_limit), since a limit dataset rests on 3-48 limit states.

Only base methods are kept (no _ext variants). Methods with fewer than 3 finite errors on a
material are left out for that material.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.homepage_data import (  # noqa: E402
    RANKING_CSV,
    RANKING_MATERIALS,
    check_manifest,
    write_manifest,
)

OUT = RANKING_CSV


# Limit datasets left out until their data errors are resolved (issue #24, decision d110); every
# method would otherwise be judged on the data error, not on its criterion.
UNDER_REVIEW = {
    "CX18": "Summary holds the 30NCD16 limits (710 / 450 MPa), not those of its rows (310 / 179)",
    "GGG60": "axial amplitudes about twice the Summary limit (ranges entered as amplitudes?)",
}


def _row(material: str, method: str, errors: np.ndarray, assessment: str) -> dict:
    return dict(
        Material=material,
        Method=method,
        n=int(errors.size),
        Error_load_mean=float(errors.mean()),
        Error_load_sd=float(errors.std(ddof=1)),
        # bias and scatter together: on a torsion-only assessment set every method has
        # the same SD (one tension/torsion ratio each), only its bias differs
        Error_load_rms=float(np.sqrt(np.mean(errors**2))),
        Within_20pct=float(np.mean(np.abs(np.expm1(errors)) <= 0.2)),
        Assessment=assessment,
    )


def is_axial_reference(case: dict, tol: float = 1e-9) -> bool:
    """True for a pure uniaxial R = -1 axial limit state (no mean, axial amplitude only)."""
    return _pure_reversed(case, (1, 1), tol)


def is_torsion_reference(case: dict, tol: float = 1e-9) -> bool:
    """True for a pure R = -1 torsion limit state (no mean, axial-hoop shear amplitude only)."""
    return _pure_reversed(case, (0, 1), tol)


def is_reference(case: dict, tol: float = 1e-9) -> bool:
    """Reference limit state of the ranking: sigma_-1 or tau_-1 itself (d118)."""
    return is_axial_reference(case, tol) or is_torsion_reference(case, tol)


def _pure_reversed(case: dict, slot: tuple[int, int], tol: float) -> bool:
    h = case["harmonics"]
    amplitude = np.hypot(np.asarray(h["S_sin"], float), np.asarray(h["S_cos"], float))
    others = amplitude.copy()
    i, j = slot
    others[i, j] = others[j, i] = 0.0
    return bool(
        float(np.abs(np.asarray(h["S_mean"], float)).max()) <= tol
        and amplitude[i, j] > tol
        and float(others.max()) <= tol
    )


def limit_errors(path: Path, reference: list[bool], specimens: list[str]) -> np.ndarray:
    """Error_load = ln(1 + FI_load) of the assessed limit states of one FI file."""
    frame = pd.read_csv(path)
    # the pipeline writes one row per run-out, in the order of the imported block
    if len(frame) != len(reference) or [str(s) for s in frame["Specimen_id"].fillna("")] != [
        str(s) for s in specimens
    ]:
        raise ValueError(f"{path.name}: rows do not match the run-outs of the workbook")
    fi = pd.to_numeric(frame["FI_load"], errors="coerce").to_numpy(float)
    errors = np.log1p(fi[~np.asarray(reference, bool)])
    return errors[np.isfinite(errors)]


def _limit_states(material: str) -> tuple[list[bool], list[str]]:
    from multifat.cli import run as runmod

    data, resolved = runmod.load_material(material)
    cases = list(data[resolved]["runouts"].values())
    return [is_reference(c) for c in cases], [str(c["meta"].get("specimen_id", "")) for c in cases]


def build(results_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(results_dir.glob("*.csv")):
        frame = pd.read_csv(path)
        if "Error_load" not in frame or "Flag" not in frame:
            continue
        method = str(frame["Method"].iloc[0])
        if method.endswith("_ext"):
            continue
        errors = pd.to_numeric(frame.loc[frame["Flag"] == 1, "Error_load"], errors="coerce")
        errors = errors[np.isfinite(errors)].to_numpy()
        if errors.size < 3:
            continue
        rows.append(_row(str(frame["Material"].iloc[0]), method, errors, "fatigue_life"))
    states: dict[str, tuple[list[bool], list[str]]] = {}
    for path in sorted((results_dir / "FatigueLimit").glob("*_FI.csv")):
        head = pd.read_csv(path, nrows=1)
        if head.empty:
            continue
        material, method = str(head["Material"].iloc[0]), str(head["Method"].iloc[0])
        if method.endswith("_ext") or material in UNDER_REVIEW:
            continue
        if material not in states:
            states[material] = _limit_states(material)
        errors = limit_errors(path, *states[material])
        if errors.size >= 3:
            rows.append(_row(material, method, errors, "fatigue_limit"))
    return pd.DataFrame(rows).sort_values(["Material", "Method"])


def write_ranking(results_dir: Path, out_path: Path | None = None) -> tuple[Path, Path]:
    """Build and save the ranking CSV along with its provenance manifest."""
    target = out_path or OUT
    table = build(results_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(target, index=False, float_format="%.6g", lineterminator="\n")
    # every material that was run, also one whose methods were all skipped or left fewer than 3
    # errors (e.g. Ti6Al4V_Meneghetti, whose master curves are not identifiable, #74): its fingerprint must
    # still be tracked
    materials = sorted(
        {str(m) for m in table["Material"].astype(str)}
        | {
            str(pd.read_csv(path, nrows=1, usecols=["Material"])["Material"].iloc[0])
            for path in [
                *results_dir.glob("*.csv"),
                *(results_dir / "FatigueLimit").glob("*_FI.csv"),
            ]
            if "Material" in pd.read_csv(path, nrows=0).columns
        }
        - set(UNDER_REVIEW)
    )
    methods = sorted(table["Method"].astype(str).unique())
    manifest = write_manifest(target, materials, methods)
    return target, manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results_dir", type=Path, nargs="?", help="pipeline output directory")
    parser.add_argument("--check", action="store_true", help="check committed CSV and manifest")
    args = parser.parse_args(argv)
    if args.check:
        problems = check_manifest(OUT, RANKING_MATERIALS)
        for problem in problems:
            print(problem)
        if problems:
            raise SystemExit(
                "stale homepage ranking data: rerun examples/build_load_error_ranking.py"
            )
        print(f"{OUT.name} is current")
        return
    if args.results_dir is None:
        parser.error("results_dir is required unless --check is passed")
    csv_path, manifest_path = write_ranking(args.results_dir)
    table = pd.read_csv(csv_path)
    print(f"{len(table)} rows -> {csv_path} (manifest: {manifest_path.name})")


if __name__ == "__main__":
    main()
