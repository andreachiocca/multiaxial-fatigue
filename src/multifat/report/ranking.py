"""Relative RMS load error ranking, complexity axis and family selection (issue #57, decision d116).

Implements the cross-material ranking KPI:
  - Accuracy KPI Z_m: average relative log-RMS load error vs the material median.
  - Complexity: explicit (p_in, p_fit) counts reported beside the accuracy KPI.
  - Family-wise 1-SE selection: simplest criterion (min p_in, p_fit) among models
    statistically indistinguishable from the family best (Z <= Z_best + SE_best).
  - Leave-one-material-out stability check.
  - Equivalence class detection.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from multifat.models.complexity import get_model_complexity

LN11 = math.log(1.1)

FAMILY_MAP: dict[str, str] = {
    "cp_methods": "Critical plane",
    "invariant_methods": "Invariant",
    "integral_methods": "Integral",
    "energy_based_methods": "Energy",
}


@dataclass(frozen=True)
class RankingResult:
    """Complete results from the ranking evaluation."""

    cell_stats: pd.DataFrame
    method_ranking: pd.DataFrame
    family_selection: pd.DataFrame
    equivalence_classes: list[list[str]]
    lomo_stability: dict[str, dict[str, int]]


def _get_model_family(method_name: str) -> str:
    """Return the display family name for a method."""
    from multifat.models.registry import available_cp_methods, available_direct_models

    cp_dict = available_cp_methods()
    direct_dict = available_direct_models()
    clean_name = method_name.removesuffix("_ext")
    if clean_name in cp_dict:
        return "Critical plane"
    if clean_name in direct_dict:
        inner = direct_dict[clean_name]
        pkg = type(inner).__module__.split(".")[2]
        return FAMILY_MAP.get(pkg, pkg)
    return "Unknown"


def compute_cell_stats(
    df: pd.DataFrame,
    *,
    min_tests: int = 3,
    exclude_ext: bool = True,
) -> pd.DataFrame:
    """Per-(Method, Material) load error statistics on assessment failures."""
    if "Error_load" not in df or "Flag" not in df:
        raise ValueError("DataFrame must contain 'Error_load' and 'Flag' columns")

    valid = df.copy()
    valid["Error_load_num"] = pd.to_numeric(valid["Error_load"], errors="coerce")
    valid = valid[np.isfinite(valid["Error_load_num"])].copy()
    valid["e"] = valid["Error_load_num"].astype(float)
    valid["Flag_num"] = pd.to_numeric(valid["Flag"], errors="coerce")

    if exclude_ext:
        valid = valid[~valid["Method"].astype(str).str.endswith("_ext")].copy()

    rows: list[dict[str, Any]] = []
    for (meth, material), g in valid.groupby(["Method", "Material"]):
        a = g[g["Flag_num"] == 1]["e"].to_numpy(dtype=float)
        d = g[g["Flag_num"] == 0]["e"].to_numpy(dtype=float)

        if a.size < min_tests:
            continue

        bias = float(np.mean(a))
        sd = float(np.std(a, ddof=1)) if a.size > 1 else 0.0
        rms = float(math.sqrt(np.mean(a**2)))
        w10 = float(np.mean(np.abs(np.expm1(a)) <= 0.10))
        w20 = float(np.mean(np.abs(np.expm1(a)) <= 0.20))
        noncons = float(np.mean(a < -LN11))
        rms_design = float(math.sqrt(np.mean(d**2))) if d.size > 0 else np.nan

        rows.append(
            {
                "Method": str(meth),
                "Material": str(material),
                "n": int(a.size),
                "bias": bias,
                "sd": sd,
                "rms": rms,
                "w10": w10,
                "w20": w20,
                "noncons": noncons,
                "rms_design": rms_design,
                "_errors": a,
            }
        )

    if not rows:
        return pd.DataFrame()

    res = pd.DataFrame(rows)
    # Relative log-RMS against material median
    med_by_mat = res.groupby("Material")["rms"].transform("median")
    res["z"] = np.log(res["rms"] / med_by_mat)
    return res


def compute_method_ranking(
    cell_df: pd.DataFrame,
    *,
    seed: int = 0,
    n_boot: int = 1000,
) -> pd.DataFrame:
    """Aggregate cell statistics across materials into the method ranking table."""
    if cell_df.empty:
        return pd.DataFrame()

    total_materials = int(cell_df["Material"].nunique())
    med_by_mat = cell_df.groupby("Material")["rms"].median().to_dict()

    # Bootstrap standard error SE(Z_m)
    rng = np.random.default_rng(seed)
    boot_se: dict[str, float] = {}

    for meth, g in cell_df.groupby("Method"):
        errors_list = list(g["_errors"])
        mats_list = list(g["Material"])
        k_m = len(errors_list)

        if k_m == 0:
            boot_se[str(meth)] = np.nan
            continue

        boot_z = []
        for _ in range(n_boot):
            draw_z = []
            for errs, material in zip(errors_list, mats_list):
                sampled = rng.choice(errs, size=errs.size, replace=True)
                rms_b = math.sqrt(np.mean(sampled**2))
                draw_z.append(math.log(rms_b / med_by_mat[material]))
            boot_z.append(float(np.mean(draw_z)))
        boot_se[str(meth)] = float(np.std(boot_z))

    records: list[dict[str, Any]] = []
    for meth, g in cell_df.groupby("Method"):
        meth_str = str(meth)
        k = int(g["Material"].nunique())
        z_mean = float(g["z"].mean())
        rms_med = float(g["rms"].median())
        bias_med = float(g["bias"].median())
        w10_med = float(g["w10"].median())
        w20_med = float(g["w20"].median())
        noncons_med = float(g["noncons"].median())

        valid_design = g["rms_design"].dropna()
        if not valid_design.empty:
            gap = float(np.median(np.log(g.loc[valid_design.index, "rms"] / valid_design)))
        else:
            gap = np.nan

        se = boot_se.get(meth_str, np.nan)
        fam = _get_model_family(meth_str)
        comp = get_model_complexity(meth_str)

        records.append(
            {
                "Method": meth_str,
                "Family": fam,
                "Z": z_mean,
                "SE": se,
                "K": k,
                "Coverage": f"{k}/{total_materials}",
                "Coverage_frac": float(k / total_materials),
                "rms_med": rms_med,
                "bias_med": bias_med,
                "w10_med": w10_med,
                "w20_med": w20_med,
                "noncons_med": noncons_med,
                "gap": gap,
                "p_in": comp.p_in,
                "p_fit": comp.p_fit,
                "input_classes": ", ".join(comp.input_classes),
                "fitted_constants": ", ".join(comp.fitted_constants),
            }
        )

    out = pd.DataFrame(records)
    return out.sort_values(["Family", "Z"]).reset_index(drop=True)


def find_equivalence_classes(
    cell_df: pd.DataFrame,
    *,
    tol: float = 1e-6,
) -> list[list[str]]:
    """Identify method equivalence classes where assessment errors match across shared materials."""
    if cell_df.empty:
        return []

    # Map each method to its error signature
    key_dict: dict[str, tuple[float, ...]] = {}
    for meth, g in cell_df.groupby("Method"):
        sorted_g = g.sort_values("Material")
        all_errs: list[float] = []
        for x in sorted_g["_errors"]:
            all_errs.extend(np.round(x / tol) * tol)
        key_dict[str(meth)] = tuple(all_errs)

    # Group by identical tuple
    groups: dict[tuple[float, ...], list[str]] = {}
    for meth, key in key_dict.items():
        groups.setdefault(key, []).append(meth)

    equiv = [sorted(g) for g in groups.values() if len(g) > 1]
    return sorted(equiv, key=lambda g: g[0])


def select_family_candidates(
    ranking_df: pd.DataFrame,
    *,
    min_coverage: int = 4,
    se_multiplier: float = 1.0,
) -> pd.DataFrame:
    """Apply the family-wise 1-SE selection rule to choose simplest top candidates."""
    if ranking_df.empty:
        return pd.DataFrame()

    eligible = ranking_df[ranking_df["K"] >= min_coverage].copy()
    if eligible.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for fam, g in eligible.groupby("Family"):
        best_idx = g["Z"].idxmin()
        best_row = g.loc[best_idx]
        best_z = float(best_row["Z"])
        best_se = float(best_row["SE"])
        threshold = best_z + se_multiplier * best_se

        candidates = g[g["Z"] <= threshold].sort_values(
            ["p_in", "p_fit", "Z"], ascending=[True, True, True]
        )
        picked_method = str(candidates.iloc[0]["Method"])

        for _, r in g.iterrows():
            m_name = str(r["Method"])
            in_candidates = bool(r["Z"] <= threshold)
            is_picked = bool(m_name == picked_method)
            rows.append(
                {
                    **r.to_dict(),
                    "Z_best_family": best_z,
                    "SE_best_family": best_se,
                    "Candidate_1SE": in_candidates,
                    "Selected": is_picked,
                }
            )

    return pd.DataFrame(rows).sort_values(["Family", "Z"]).reset_index(drop=True)


def leave_one_material_out_stability(
    cell_df: pd.DataFrame,
    *,
    min_coverage: int = 4,
    seed: int = 0,
    n_boot: int = 500,
) -> dict[str, dict[str, int]]:
    """Evaluate stability of the family-wise selection when leaving one material out."""
    if cell_df.empty:
        return {}

    all_materials = sorted(cell_df["Material"].unique())
    if len(all_materials) <= 2:
        return {}

    sub_coverage = max(1, min_coverage - 1)
    picks_per_family: dict[str, dict[str, int]] = {}

    for out_mat in all_materials:
        sub_cells = cell_df[cell_df["Material"] != out_mat].copy()
        if sub_cells.empty:
            continue

        # Re-compute z on subset
        sub_med = sub_cells.groupby("Material")["rms"].transform("median")
        sub_cells["z"] = np.log(sub_cells["rms"] / sub_med)

        sub_ranking = compute_method_ranking(sub_cells, seed=seed, n_boot=n_boot)
        sub_selection = select_family_candidates(sub_ranking, min_coverage=sub_coverage)

        if sub_selection.empty:
            continue

        for _, r in sub_selection[sub_selection["Selected"]].iterrows():
            fam = str(r["Family"])
            meth = str(r["Method"])
            family_dict = picks_per_family.setdefault(fam, {})
            family_dict[meth] = family_dict.get(meth, 0) + 1

    return picks_per_family


def evaluate_ranking(
    df: pd.DataFrame,
    *,
    min_tests: int = 3,
    min_coverage: int = 4,
    seed: int = 0,
    n_boot: int = 1000,
) -> RankingResult:
    """Run full ranking pipeline on a results DataFrame."""
    cell_stats = compute_cell_stats(df, min_tests=min_tests)
    method_ranking = compute_method_ranking(cell_stats, seed=seed, n_boot=n_boot)
    family_selection = select_family_candidates(method_ranking, min_coverage=min_coverage)
    equivalence = find_equivalence_classes(cell_stats)
    stability = leave_one_material_out_stability(
        cell_stats, min_coverage=min_coverage, seed=seed, n_boot=min(n_boot, 500)
    )

    return RankingResult(
        cell_stats=cell_stats,
        method_ranking=method_ranking,
        family_selection=family_selection,
        equivalence_classes=equivalence,
        lomo_stability=stability,
    )
