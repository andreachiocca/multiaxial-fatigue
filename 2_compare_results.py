#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors
import numpy as np
import pandas as pd


# =========================
# Utilities
# =========================
def _safe_name(s: str) -> str:
    """Filesystem-safe name (keeps underscores)."""
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in str(s))


_LATEX_SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def _latex_escape(s: str) -> str:
    """Escape a string for LaTeX text mode (keeps it human-readable)."""
    s = str(s)
    return "".join(_LATEX_SPECIALS.get(ch, ch) for ch in s)


def _latex_method_id(method: str) -> str:
    """
    A pgfplots/key-safe identifier for a method.
    - only letters/digits
    - starts with a letter
    - deterministic (collisions handled elsewhere)
    """
    base = re.sub(r"[^A-Za-z0-9]+", "", str(method))
    if not base:
        base = "method"
    if base[0].isdigit():
        base = "m" + base
    return base[:32]


def _normal_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
        return np.full_like(x, np.nan, dtype=float)
    return (1.0 / (sigma * math.sqrt(2.0 * math.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    norm = {c: c.strip().lower().replace(" ", "_") for c in df.columns}
    inv = {v: k for k, v in norm.items()}
    for cand in candidates:
        key = cand.strip().lower().replace(" ", "_")
        if key in inv:
            return inv[key]
    return None


def _infer_from_filename(path: Path) -> Tuple[Optional[str], Optional[str]]:
    stem = path.stem
    parts = stem.split("_")
    if len(parts) < 2:
        return None, None

    def looks_like_material(s: str) -> bool:
        return any(ch.isdigit() for ch in s) and any(ch.isalpha() for ch in s)

    left = "_".join(parts[:-1])
    right = parts[-1]
    if looks_like_material(right):
        return right, left
    if looks_like_material(parts[0]):
        return parts[0], "_".join(parts[1:])
    return None, None


# =========================
# Style system (ColorBrewer-like palette, consistent across plots)
# =========================
# ColorBrewer "Paired" (12) — good separation and print-friendly
_BREWER_PAIRED_12 = [
    "A6CEE3",
    "1F78B4",
    "B2DF8A",
    "33A02C",
    "FB9A99",
    "E31A1C",
    "FDBF6F",
    "FF7F00",
    "CAB2D6",
    "6A3D9A",
    "FFFF99",
    "B15928",
]

# Common pgfplots marks (repeat if more methods than marks)
_PGFPLOTS_MARKS = [
    "*",
    "square*",
    "triangle*",
    "diamond*",
    "pentagon*",
    "x",
    "+",
    "star",
    "o",
    "triangle",
    "square",
    "diamond",
]


def _unique_key_id(base: str, used: set[str]) -> str:
    """Return a unique, deterministic key id (case-insensitive) given a base id."""
    key = base
    k = 1
    while key.lower() in used:
        key = f"{base}{k}"
        k += 1
    used.add(key.lower())
    return key


def pgf_mark_to_mpl(mark: str) -> str:
    """Best-effort mapping from pgfplots mark names to Matplotlib markers."""
    m = str(mark or "").strip()
    return {
        "*": "o",
        "o": "o",
        "square*": "s",
        "square": "s",
        "triangle*": "^",
        "triangle": "^",
        "diamond*": "D",
        "diamond": "D",
        "pentagon*": "p",
        "x": "x",
        "+": "+",
        "star": "*",
    }.get(m, "o")


def _family_cmap_name(family: str) -> str:
    """
    Map a method family name to a ColorBrewer-like sequential palette.
    The default Matplotlib sequential colormaps (Reds, Blues, Greens, …) are based on ColorBrewer.
    """
    f = (family or "UNKNOWN").lower()
    # Your main request:
    # - CP -> red palette
    # - invariants -> blue palette
    if "cp" in f or "critical" in f:
        return "Reds"
    if "invariant" in f:
        return "Blues"
    # Sensible defaults for other families
    if "energy" in f:
        return "Greens"
    if "stress" in f:
        return "Purples"
    if "strain" in f:
        return "Oranges"
    if "unknown" in f:
        return "Greys"
    return "Greys"


def _cmap_sample_hex(cmap_name: str, idx: int, n: int) -> str:
    """
    Sample a sequential colormap to get a readable, print-friendly color.
    Avoid the very light end to keep visibility on white backgrounds.
    """
    cmap = cm.get_cmap(cmap_name)
    if n <= 1:
        t = 0.70
    else:
        # keep colors away from the very light/very dark extremes
        t_values = np.linspace(0.40, 0.90, n)
        t = float(t_values[idx])
    rgb = cmap(t)
    # mcolors.to_hex -> '#RRGGBB'
    return mcolors.to_hex(rgb, keep_alpha=False).lstrip("#").upper()


def build_style_maps(df: pd.DataFrame) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    """
    Build deterministic style mappings using Method_family (method family):

      family_styles: family -> {id, cmap}
      method_styles: method -> {id, family, base, is_ext, hex, mark, color}

    BEHAVIOR:
      - Each FAMILY gets its own sequential palette (e.g. CP=Reds, invariants=Blues).
      - Each BASE METHOD (method name with optional "_ext" stripped) inside the family gets a different shade.
      - If a method has an "_ext" version, the base method and its "_ext" share the SAME color and marker,
        but differ by line style / marker fill:
          * base: solid line, filled marker
          * _ext: densely dashed line, empty marker
      - Mapping is deterministic: sorted families + sorted base methods per family.
    """
    if "Method" not in df.columns:
        raise ValueError("DataFrame must contain a 'Method' column.")
    if "Method_family" not in df.columns:
        df = df.copy()
        df["Method_family"] = "UNKNOWN"
    df["Method_family"] = df["Method_family"].fillna("UNKNOWN").astype(str)

    def _base_method_name(method: str) -> str:
        return method[:-4] if isinstance(method, str) and method.endswith("_ext") else str(method)

    def _is_ext(method: str) -> bool:
        return isinstance(method, str) and method.endswith("_ext")

    # --- Families -> colormap ---
    families_sorted = sorted([f for f in df["Method_family"].dropna().unique().tolist() if isinstance(f, str)])
    fam_used: set[str] = set()
    family_styles: Dict[str, Dict[str, str]] = {}
    for fam in families_sorted:
        base = _latex_method_id(fam)
        fid = _unique_key_id(base, fam_used)
        family_styles[fam] = {"id": fid, "cmap": _family_cmap_name(fam)}

    if "UNKNOWN" not in family_styles:
        fid = _unique_key_id(_latex_method_id("UNKNOWN"), fam_used)
        family_styles["UNKNOWN"] = {"id": fid, "cmap": _family_cmap_name("UNKNOWN")}

    # --- Methods -> per-family shades (by BASE method) + per-base marker ---
    methods_sorted = sorted([m for m in df["Method"].dropna().unique().tolist() if isinstance(m, str)])

    # method -> family (first non-null)
    meth_to_family: Dict[str, str] = {}
    for mth in methods_sorted:
        fam_series = df.loc[df["Method"] == mth, "Method_family"].dropna()
        fam = str(fam_series.iloc[0]) if len(fam_series) else "UNKNOWN"
        if fam not in family_styles:
            fam = "UNKNOWN"
        meth_to_family[mth] = fam

    # family -> base methods
    family_to_bases: Dict[str, List[str]] = {}
    for mth in methods_sorted:
        fam = meth_to_family.get(mth, "UNKNOWN")
        base = _base_method_name(mth)
        family_to_bases.setdefault(fam, []).append(base)

    for fam in family_to_bases:
        family_to_bases[fam] = sorted(set(family_to_bases[fam]))

    # assign colors per base in each family
    base_to_hex: Dict[Tuple[str, str], str] = {}
    for fam, bases in family_to_bases.items():
        cmap_name = family_styles.get(fam, family_styles["UNKNOWN"]).get("cmap", "Greys")
        for j, base in enumerate(bases):
            base_to_hex[(fam, base)] = _cmap_sample_hex(cmap_name, j, len(bases))

    # assign markers per base deterministically across ALL bases
    all_bases_sorted: List[str] = sorted({ _base_method_name(mth) for mth in methods_sorted })
    base_to_mark: Dict[str, str] = {}
    for i, base in enumerate(all_bases_sorted):
        base_to_mark[base] = _PGFPLOTS_MARKS[i % len(_PGFPLOTS_MARKS)]

    meth_used: set[str] = set()
    method_styles: Dict[str, Dict[str, str]] = {}
    for mth in methods_sorted:
        base = _base_method_name(mth)
        fam = meth_to_family.get(mth, "UNKNOWN")
        hex_ = base_to_hex.get((fam, base), _cmap_sample_hex("Greys", 0, 1))
        mark = base_to_mark.get(base, _PGFPLOTS_MARKS[0])

        mid = _unique_key_id(_latex_method_id(mth), meth_used)
        colorname = f"faberMethodColor{mid}"
        method_styles[mth] = {
            "id": mid,
            "family": fam,
            "base": base,
            "is_ext": "1" if _is_ext(mth) else "0",
            "hex": hex_,
            "mark": mark,
            "color": colorname,
        }

    return method_styles, family_styles


def write_faber_style_file(
    method_styles: Dict[str, Dict[str, str]],
    family_styles: Dict[str, Dict[str, str]],
    out_path: Path,
) -> None:
    """
    Write a LaTeX file with color + pgfplots styles.

    Behavior:
      - Colors are defined per METHOD (shade depends on its method family palette).
      - If a method has an "_ext" version, the base method and its "_ext" share the SAME color and marker,
        but differ by:
          * base: solid line, filled marker
          * _ext: densely dashed line, empty marker
      - Boxplots: border-only boxes; outliers are always black 'x'.
      - Global font size is set to \\scriptsize.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines += [
        r"\providecommand{\faberPlotStylesLoaded}{0}",
        r"\ifnum\faberPlotStylesLoaded=0",
        r"  \renewcommand{\faberPlotStylesLoaded}{1}",
        r"  \tikzset{every picture/.style={font=\scriptsize}}",
        "",
    ]

    # Define method colors
    for mth in sorted(method_styles.keys()):
        st = method_styles[mth]
        lines.append(rf"  \definecolor{{{st['color']}}}{{HTML}}{{{st['hex']}}}")

    lines += [
        "",
        r"  \pgfplotsset{",
        r"    every axis/.append style={font=\scriptsize},",
        r"    every axis legend/.append style={font=\scriptsize},",
        r"    % Boxplots: no fill + black x outliers",
        r"    boxplot/every box/.style={solid, fill=none},",
        r"    boxplot/every whisker/.style={solid},",
        r"    boxplot/every median/.style={solid},",
        r"    boxplot/every outlier/.style={only marks, mark=x, draw=black, mark size=2.2pt},",
        r"    boxplot outlier/.style={only marks, mark=x, draw=black, mark size=2.2pt},",
    ]

    # Method styles
    for mth in sorted(method_styles.keys()):
        st = method_styles[mth]
        mid = st["id"]
        col = st["color"]
        mark = st["mark"]
        is_ext = st.get("is_ext", "0") == "1"

        # Base color style
        lines.append(rf"    faber/method/{mid}/.style={{color={col}}},")

        # Scatter style:
        #  - base: filled marker
        #  - _ext: empty marker
        fillopt = col if not is_ext else "white"
        lines.append(
            rf"    faber/method_scatter/{mid}/.style={{faber/method/{mid}, only marks, mark={mark}, mark size=1.6pt, "
            rf"draw={col}, mark options={{solid, fill={fillopt}}}}},"
        )

        # Line style:
        #  - base: solid
        #  - _ext: densely dashed
        linestyle = "densely dashed" if is_ext else "solid"
        lines.append(rf"    faber/method_line/{mid}/.style={{faber/method/{mid}, {linestyle}, line width=1.1pt}},")

        # Boxplot style (border only)
        lines.append(
            rf"    faber/method_box/{mid}/.style={{faber/method/{mid}, boxplot, boxplot/draw direction=y, "
            rf"draw={col}, fill=none}},"
        )

    lines += [
        r"  }",
        r"\fi",
        "",
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")



# =========================
# CSV loading
# =========================
def read_results_csv(csv_path: Path) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[WARN] Could not read {csv_path}: {e}")
        return None

    if df is None or df.empty:
        return None

    df.columns = [c.strip() for c in df.columns]

    col_material = _find_col(df, ["Material", "material", "mat"])
    col_method = _find_col(df, ["Method", "method", "model"])
    col_family = _find_col(df, ["Method_family", "MethodFamily", "Family", "Method_type", "MethodType", "Type"])
    col_metric = _find_col(df, ["Metric_mode", "Metric", "Mode", "Error_metric", "ErrorMetric"])
    col_nfexp = _find_col(df, ["Nf_exp", "Nexp", "Nfexperimental", "Nf_experimental"])
    col_nfpred = _find_col(df, ["Nf_expected", "Npred", "Nfexpected", "Nf_expected_cycles"])
    col_cp = _find_col(df, ["CP_value", "CP", "DamageParameter", "Damage_parameter", "DP", "dp"])
    col_dpfit = _find_col(df, ["DP_fit", "DPfit", "CP_fit", "DP_fit_value"])
    col_err = _find_col(df, ["Error_ln", "error_ln", "Error", "error"])
    col_errlndp = _find_col(df, ["Error_ln_dp", "ErrorLnDP", "Error_lnDP", "Error_lnDP", "Error_dp_ln"])
    # Deprecated/legacy (difference-based) DP error
    col_errdp = _find_col(df, ["Error_dp", "error_dp", "ErrorDP", "Error_DP"])
    col_flag = _find_col(df, ["Flag", "flag"])

    # Infer missing material/method from filename
    if col_material is None or col_method is None:
        mat, meth = _infer_from_filename(csv_path)
        if col_material is None:
            df["Material"] = mat if mat is not None else "UNKNOWN"
            col_material = "Material"
        if col_method is None:
            df["Method"] = meth if meth is not None else "UNKNOWN"
            col_method = "Method"

    # Rename to canonical names
    rename = {}
    if col_material and col_material != "Material":
        rename[col_material] = "Material"
    if col_method and col_method != "Method":
        rename[col_method] = "Method"
    if col_family and col_family != "Method_family":
        rename[col_family] = "Method_family"
    if col_metric and col_metric != "Metric_mode":
        rename[col_metric] = "Metric_mode"
    if col_nfexp and col_nfexp != "Nf_exp":
        rename[col_nfexp] = "Nf_exp"
    if col_nfpred and col_nfpred != "Nf_expected":
        rename[col_nfpred] = "Nf_expected"
    if col_cp and col_cp != "CP_value":
        rename[col_cp] = "CP_value"
    if col_dpfit and col_dpfit != "DP_fit":
        rename[col_dpfit] = "DP_fit"
    if col_err and col_err != "Error_ln":
        rename[col_err] = "Error_ln"
    if col_errlndp and col_errlndp != "Error_ln_dp":
        rename[col_errlndp] = "Error_ln_dp"
    if col_errdp and col_errdp != "Error_dp":
        rename[col_errdp] = "Error_dp"
    if col_flag and col_flag != "Flag":
        rename[col_flag] = "Flag"

    df = df.rename(columns=rename)

    # Numeric conversions

    # Ensure method family exists
    if "Method_family" not in df.columns:
        df["Method_family"] = "UNKNOWN"
    else:
        df["Method_family"] = df["Method_family"].fillna("UNKNOWN").astype(str)

    # Ensure metric mode exists
    if "Metric_mode" not in df.columns:
        df["Metric_mode"] = "Nf"
    else:
        df["Metric_mode"] = df["Metric_mode"].fillna("Nf").astype(str)

    # Numeric conversions
    for c in ["Nf_exp", "Nf_expected", "CP_value", "DP_fit", "Error_ln", "Error_ln_dp", "Error_dp", "Flag"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Compute Error_ln if missing and Nf_expected is available
    if "Error_ln" not in df.columns and "Nf_exp" in df.columns and "Nf_expected" in df.columns:
        df["Error_ln"] = np.nan
        m = (df["Nf_exp"] > 0) & (df["Nf_expected"] > 0) & np.isfinite(df["Nf_exp"]) & np.isfinite(df["Nf_expected"])
        df.loc[m, "Error_ln"] = np.log(df.loc[m, "Nf_exp"] / df.loc[m, "Nf_expected"])

    # Compute DP log-ratio error if missing and DP_fit + CP_value are available
    if "Error_ln_dp" not in df.columns:
        df["Error_ln_dp"] = np.nan
    if "DP_fit" in df.columns and "CP_value" in df.columns:
        m = np.isfinite(df["DP_fit"]) & np.isfinite(df["CP_value"]) & (df["DP_fit"] > 0.0) & (df["CP_value"] > 0.0)
        df.loc[m, "Error_ln_dp"] = np.log(df.loc[m, "CP_value"] / df.loc[m, "DP_fit"])

    # Legacy: compute difference-based Error_dp only if present (or requested by older CSVs)
    if "Error_dp" not in df.columns:
        df["Error_dp"] = np.nan
    if "DP_fit" in df.columns and "CP_value" in df.columns:
        m = np.isfinite(df["DP_fit"]) & np.isfinite(df["CP_value"])
        df.loc[m, "Error_dp"] = df.loc[m, "CP_value"] - df.loc[m, "DP_fit"]

    df["source_file"] = str(csv_path)
    return df


def collect_all_results(results_dir: Path) -> pd.DataFrame:
    csv_files = sorted([p for p in results_dir.rglob("*.csv") if p.is_file() and "Comparison" not in p.parts])
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found under: {results_dir}")

    dfs = []
    for p in csv_files:
        dfi = read_results_csv(p)
        if dfi is not None and not dfi.empty:
            dfs.append(dfi)

    if not dfs:
        raise RuntimeError(f"CSV files found under {results_dir}, but none contained readable data.")
    return pd.concat(dfs, ignore_index=True)


# =========================
# TikZ helpers
# =========================
def _write_dat_1col(path: Path, y: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    with open(path, "w", encoding="utf-8") as f:
        for v in y:
            f.write(f"{v:.12e}\n")


def _write_dat_xy(path: Path, x: np.ndarray, y: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    with open(path, "w", encoding="utf-8") as f:
        for xi, yi in zip(x, y):
            f.write(f"{xi:.12e} {yi:.12e}\n")


def _tikz_header(comment: str) -> str:
    return (
        "% Auto-generated by 2_compare_results.py\n"
        f"% {comment}\n"
        "% NOTE on paths: data file paths are written with '_' (not '\\_') on purpose.\n"
        "% Requires in LaTeX preamble:\n"
        "%   \\usepackage{xcolor}\n"
        "%   \\usepackage{tikz}\n"
        "%   \\usepackage{pgfplots}\n"
        "%   \\pgfplotsset{compat=1.18}\n"
        "%   \\usepgfplotslibrary{statistics}\n\n"
    )


def _factor_pgf_linestyle(i: int) -> str:
    # Ensure Factor 3 vs Factor 5 are distinguishable (densely dashed vs densely dotted)
    styles = ["densely dashed", "densely dotted", "densely dashdotted", "dashdotted", "dashed", "dotted"]
    return styles[i % len(styles)]


def _factor_mpl_linestyle(i: int) -> str:
    styles = ["--", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 3))]
    return styles[i % len(styles)]


# =========================
# Plotters
# =========================
def plot_boxplot_errors(
    df_mat: pd.DataFrame,
    out_png: Path,
    out_tex: Path,
    data_dir: Path,
    material_label: str,
    method_styles: Dict[str, Dict[str, str]],
    *,
    error_col: str = "Error_ln",
    ylabel: Optional[str] = None,
    show_method_legend_table_hint: bool = True,
) -> None:
    methods = sorted(df_mat["Method"].dropna().unique().tolist())

    if ylabel is None:
        if error_col == "Error_ln_dp":
            ylabel = r"$Error_{\ln,dp}=\ln\left(\frac{DP}{DP_{fit}}\right)$"
        else:
            ylabel = r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$"

    # Build (median, method, values) so we can sort boxplots from highest median to lowest.
    items: List[Tuple[float, str, np.ndarray]] = []
    for m in methods:
        vals = df_mat.loc[df_mat["Method"] == m, error_col].dropna().to_numpy(dtype=float)
        if vals.size > 0:
            med = float(np.nanmedian(vals)) if np.any(np.isfinite(vals)) else np.nan
            items.append((med, m, vals))

    if not items:
        print(f"[WARN] No {error_col} for boxplot: {material_label}")
        return

    # Sort: highest median first; NaNs last; stable tie-break by method name.
    items.sort(key=lambda t: (0 if np.isfinite(t[0]) else 1, -t[0] if np.isfinite(t[0]) else 0.0, t[1]))

    labels = [m for _, m, _ in items]
    data = [v for _, _, v in items]

    # Matplotlib (colors consistent with TikZ)
    fig, ax = plt.subplots(figsize=(10, 5))
    flierprops = dict(
        marker="x",
        markerfacecolor="none",
        markeredgecolor="black",
        markeredgewidth=1.0,
        markersize=5,
        linestyle="none",
    )
    positions = np.arange(1, len(data) + 1)
    bp = ax.boxplot(
        data,
        positions=positions,
        labels=labels,
        showfliers=True,
        patch_artist=True,
        flierprops=flierprops,
    )
    for i, (box, mth) in enumerate(zip(bp["boxes"], labels)):
        st = method_styles.get(mth)
        if st:
            col = "#" + st["hex"]
            # Border only (no fill)
            box.set_facecolor("none")
            box.set_edgecolor(col)
            box.set_linewidth(1.2)
            # Match whiskers/caps/median to the same color
            for j in (2 * i, 2 * i + 1):
                if j < len(bp.get("whiskers", [])):
                    bp["whiskers"][j].set_color(col)
                if j < len(bp.get("caps", [])):
                    bp["caps"][j].set_color(col)
            if i < len(bp.get("medians", [])):
                bp["medians"][i].set_color(col)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # TikZ
    dat_files = []
    for m in labels:
        vals = df_mat.loc[df_mat["Method"] == m, error_col].dropna().to_numpy(dtype=float)
        dat = data_dir / "boxplot" / f"{_safe_name(m)}.dat"
        _write_dat_1col(dat, vals)
        dat_files.append((m, dat))

    xticks = ",".join(str(i + 1) for i in range(len(dat_files)))
    xticklabels = ",".join(_latex_escape(m) for m, _ in dat_files)

    tex = (
        _tikz_header(f"Boxplot of {error_col} by method")
        + r"\input{Results/Comparison/faber_plot_styles.tex}"
        + "\n"
        + r"""\begin{tikzpicture}
\begin{axis}[
  width=0.95\linewidth,
  height=0.55\linewidth,
  ylabel={"""
        + ylabel
        + r"""},
  boxplot/draw direction=y,
  xtick={"""
        + xticks
        + r"""},
  xticklabels={"""
        + xticklabels
        + r"""},
  x tick label style={rotate=30, anchor=east},
]
"""
        + r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
        + _latex_escape(material_label)
        + r"};"
        + "\n"
    )

    for i, (m, dat) in enumerate(dat_files, start=1):
        st = method_styles.get(m)
        mid = st["id"] if st else None
        style = f"faber/method_box/{mid}" if mid else "boxplot"
        # IMPORTANT: keep '_' in the path (no LaTeX escaping here by user request)
        tex += (
            rf"\addplot[ {style}, forget plot, boxplot/draw position={i} ] table[y index=0] {{{dat.as_posix()}}};"
            + "\n"
        )

    if show_method_legend_table_hint:
        tex += r"% NOTE: method color/marker mapping is in Results/Comparison/summary_error_stats_table.tex" + "\n"

    tex += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(tex, encoding="utf-8")


def plot_normal_curves(
    df_mat: pd.DataFrame,
    out_png: Path,
    out_tex: Path,
    data_dir: Path,
    material_label: str,
    method_styles: Dict[str, Dict[str, str]],
    pdf_alpha: float,
    pdf_lw: float,
    show_method_legends: bool,
    *,
    error_col: str = "Error_ln",
    xlabel: Optional[str] = None,
) -> None:
    """
    Plot Normal PDF curves of an error metric for each method (mean/std from data).

    Requirements:
      - No grid.
      - y-axis label: "Probability Density Function (PDF)".
      - If show_method_legends is False, add 'forget plot' so methods don't appear in legend.
      - Safe/Not-safe text is present (no background), placed ~20% lower on y.
      - TikZ code starts with:
            \begin{tikzpicture}
            \input{Results/Comparison/faber_plot_styles.tex}
            \begin{axis}[
    """
    if xlabel is None:
        if error_col == "Error_ln_dp":
            xlabel = r"$Error_{\ln,dp}=\ln\left(\frac{DP}{DP_{fit}}\right)$"
        else:
            xlabel = r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$"

    if error_col not in df_mat.columns:
        print(f"[WARN] Missing {error_col} for normal curves: {material_label}")
        return

    # Required paper style: all PDF curves use opacity=0.1
    pdf_alpha = 0.1
    pdf_lw = float(pdf_lw)

    stats: List[Tuple[str, float, float]] = []
    for meth in sorted(df_mat["Method"].dropna().unique().tolist()):
        vals = pd.to_numeric(df_mat.loc[df_mat["Method"] == meth, error_col], errors="coerce").dropna().to_numpy(dtype=float)
        if len(vals) >= 2:
            mu = float(np.nanmean(vals))
            sigma = float(np.nanstd(vals, ddof=1))
            if np.isfinite(mu) and np.isfinite(sigma) and sigma > 0:
                stats.append((meth, mu, sigma))

    if not stats:
        print(f"[WARN] Not enough {error_col} for normal curves: {material_label}")
        return

    xmin = min(mu - 4 * s for _, mu, s in stats)
    xmax = max(mu + 4 * s for _, mu, s in stats)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin == xmax:
        xmin, xmax = -3.0, 3.0
    xs = np.linspace(xmin, xmax, 600)

    # Matplotlib
    fig, ax = plt.subplots(figsize=(10, 5))
    ymax = 0.0
    for meth, mu, sigma in stats:
        ys = _normal_pdf(xs, mu, sigma)
        ymax = max(ymax, float(np.nanmax(ys)))
        st = method_styles.get(meth)
        col = ("#" + st["hex"]) if st else None
        label = f"{meth} (μ={mu:.3g}, σ={sigma:.3g})" if show_method_legends else None
        ls = (0, (3, 1)) if (st and st.get("is_ext") == "1") else "-"
        ax.plot(xs, ys, linewidth=pdf_lw, alpha=0.1, color=col, linestyle=ls, label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability Density Function (PDF)")
    ax.grid(False)
    if show_method_legends:
        ax.legend(fontsize=9)

    # Safe / Not-safe labels (positive = conservative), moved ~20% lower on y
    if np.isfinite(ymax) and ymax > 0:
        ylab = 0.56 * ymax
        x_safe = xmax * 0.72 if xmax > 0 else xmax * 0.92
        x_not = xmin * 0.72 if xmin < 0 else xmin * 0.92
        ax.text(x_not, ylab, "Not-safe", ha="center", va="center")
        ax.text(x_safe, ylab, "Safe", ha="center", va="center")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # TikZ
    curve_files: List[Tuple[str, float, float, Path]] = []
    ymax = 0.0
    for meth, mu, sigma in stats:
        ys = _normal_pdf(xs, mu, sigma)
        ymax = max(ymax, float(np.nanmax(ys)))
        dat = data_dir / "pdf" / f"{_safe_name(meth)}.dat"
        _write_dat_xy(dat, xs, ys)
        curve_files.append((meth, mu, sigma, dat))

    tex = (r"""\begin{tikzpicture}
\input{Results/Comparison/faber_plot_styles.tex}
\begin{axis}[
  width=0.5\textwidth,
  height=0.5\textwidth,
  xlabel={""" + xlabel + r"""},
  ylabel={Probability Density Function (PDF)},
""")

    if show_method_legends:
        tex += r"  legend style={}," + "\n"

    tex += r"]" + "\n"
    tex += (
        r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
        + _latex_escape(material_label)
        + r"};"
        + "\n"
    )

    forget = "" if show_method_legends else ", forget plot"

    for meth, mu, sigma, dat in curve_files:
        st = method_styles.get(meth)
        mid = st["id"] if st else None
        style = f"faber/method_line/{mid}" if mid else ""
        tex += (
            r"\addplot["
            + style
            + forget
            + r", opacity=0.1"
            + r", line width="
            + f"{pdf_lw:.3g}"
            + r"pt] table[x index=0, y index=1] {"
            + dat.as_posix()
            + r"};"
            + "\n"
        )

        if show_method_legends:
            label = f"{_latex_escape(meth)} $(\\mu={mu:.3g},\\ \\sigma={sigma:.3g})$"
            tex += r"\addlegendentry{" + label + r"}" + "\n"

    # Safe / Not-safe labels (positive = conservative), moved ~20% lower on y
    if np.isfinite(ymax) and ymax > 0:
        ylab = 0.56 * ymax
        x_safe = xmax * 0.72 if xmax > 0 else xmax * 0.92
        x_not = xmin * 0.72 if xmin < 0 else xmin * 0.92
        tex += rf"\node[anchor=south] at (axis cs:{x_not:.6g},{ylab:.6g}) {{Not-safe}};" + "\n"
        tex += rf"\node[anchor=south] at (axis cs:{x_safe:.6g},{ylab:.6g}) {{Safe}};" + "\n"

    tex += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(tex.strip(), encoding="utf-8")


def plot_expected_vs_experimental(
    df_mat: pd.DataFrame,
    out_png: Path,
    out_tex: Path,
    data_dir: Path,
    material_label: str,
    method_styles: Dict[str, Dict[str, str]],
    scatter_alpha: float,
    factors: List[float],
    show_method_legends: bool,
) -> None:
    if "Nf_expected" not in df_mat.columns or "Nf_exp" not in df_mat.columns:
        print(f"[WARN] Missing Nf_expected/Nf_exp for scatter: {material_label}")
        return

    # Required paper style: all scatter markers have opacity=0.2 and fill opacity=0.2
    scatter_alpha = 0.2
    factors = sorted(set(float(f) for f in factors if np.isfinite(f) and f > 1.0))

    methods = sorted(df_mat["Method"].dropna().unique().tolist())

    xall = pd.to_numeric(df_mat["Nf_exp"], errors="coerce").to_numpy(dtype=float)
    yall = pd.to_numeric(df_mat["Nf_expected"], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(xall) & np.isfinite(yall) & (xall > 0) & (yall > 0)
    if not np.any(m):
        print(f"[WARN] No positive finite data for scatter: {material_label}")
        return

    xall = xall[m]
    yall = yall[m]

    # Equal range with padding in log-space
    # Also include fixed Safe/Not-safe label positions requested by the paper style.
    x_safe, y_safe = 3_000_000.0, 200.0
    x_not, y_not = 200.0, 3_000_000.0

    lo_raw = float(min(xall.min(), yall.min(), x_safe, y_safe, x_not, y_not))
    hi_raw = float(max(xall.max(), yall.max(), x_safe, y_safe, x_not, y_not))

    logmin = float(np.log10(lo_raw))
    logmax = float(np.log10(hi_raw))
    span = (logmax - logmin) if logmax > logmin else 1.0
    pad = 0.15 * span + 0.05
    lo = 10 ** (logmin - pad)
    hi = 10 ** (logmax + pad)

    # Greys for factors: keep consistent, but distinguish by line style as well
    grey_levels = [0.45, 0.65, 0.80, 0.30, 0.55]
    factor_grey: Dict[float, float] = {f: grey_levels[i % len(grey_levels)] for i, f in enumerate(factors)}

    # --- Matplotlib ---
    fig, ax = plt.subplots(figsize=(7, 7))

    for meth in methods:
        sub = df_mat[df_mat["Method"] == meth]
        x = pd.to_numeric(sub["Nf_exp"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["Nf_expected"], errors="coerce").to_numpy(dtype=float)
        mm = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if np.any(mm):
            st = method_styles.get(meth)
            col = ("#" + st["hex"]) if st else None
            label = meth if show_method_legends else "_nolegend_"
            mk = pgf_mark_to_mpl(st["mark"]) if st else "o"
            is_ext = bool(st and st.get("is_ext") == "1")
            fc = "none" if is_ext or not col else col
            ax.scatter(
                x[mm],
                y[mm],
                label=label,
                alpha=0.2,
                s=30,
                marker=mk,
                facecolors=fc,
                edgecolors=col,
                linewidths=0.9,
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    xs = np.array([lo, hi], dtype=float)

    # Bisector (solid black) in legend
    ax.plot(xs, xs, color="black", linewidth=1.8, label="Bisector")

    # Factors: distinguish styles (dense dashed vs dense dotted ...)
    for i, f in enumerate(factors):
        g = str(factor_grey[f])  # grayscale string
        f_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        ls = _factor_mpl_linestyle(i)
        ax.plot(xs, f * xs, linestyle=ls, color=g, linewidth=1.3, label=f_label)
        ax.plot(xs, xs / f, linestyle=ls, color=g, linewidth=1.3, label="_nolegend_")

    # Safe / Not-safe
    ax.text(
        x_safe,
        y_safe,
        "Safe",
        rotation=45,
        color="black",
        fontsize=11,
        ha="right",
        va="bottom",
    )
    ax.text(
        x_not,
        y_not,
        "Not-safe",
        rotation=45,
        color="black",
        fontsize=11,
        ha="left",
        va="top",
    )

    ax.set_xlabel(r"Experimental number of cycles to failure - $N_f$")
    ax.set_ylabel(r"Expected number of cycles to failure - $N_{f,e}$")
    ax.legend(fontsize=9, loc="lower left")  # only Bisector + Factors (and methods if enabled)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # --- TikZ export ---
    scat_files = []
    for meth in methods:
        sub = df_mat[df_mat["Method"] == meth]
        x = pd.to_numeric(sub["Nf_exp"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["Nf_expected"], errors="coerce").to_numpy(dtype=float)
        mm = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if np.any(mm):
            dat = data_dir / "scatter" / f"{_safe_name(meth)}.dat"
            _write_dat_xy(dat, x[mm], y[mm])
            scat_files.append((meth, dat))

    def pgf_gray(level01: float) -> str:
        p = int(round(level01 * 100))
        return f"gray!{p}"


    tex = (r"""\begin{tikzpicture}
\input{Results/Comparison/faber_plot_styles.tex}
\begin{axis}[
  width=0.5\textwidth,
  height=0.5\textwidth,
  xmode=log,
  ymode=log,
  xmin="""
        + f"{lo:.6e}"
        + r""", xmax="""
        + f"{hi:.6e}"
        + r""",
  ymin="""
        + f"{lo:.6e}"
        + r""", ymax="""
        + f"{hi:.6e}"
        + r""",
  axis equal image,
  grid=both,
  xlabel={Experimental number of cycles to failure - $N_f$},
  ylabel={Expected number of cycles to failure - $N_{f,e}$},
  legend style={at={(0.02,0.02)}, anchor=south west, font=\tiny},
]
"""
        + r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
        + _latex_escape(material_label)
        + r"};"
        + "\n"
    )

    forget = "" if show_method_legends else ", forget plot"

    for meth, dat in scat_files:
        st = method_styles.get(meth)
        mid = st["id"] if st else None
        style = f"faber/method_scatter/{mid}" if mid else "only marks"
        tex += (
            r"\addplot["
            + style
            + forget
            + r", opacity=0.2, fill opacity=0.2"
            + r"] table[x index=0, y index=1] {"
            + dat.as_posix()
            + r"};"
            + "\n"
        )
        if show_method_legends:
            tex += r"\addlegendentry{" + _latex_escape(meth) + r"}" + "\n"

    # Bisector with legend entry
    tex += (
        r"\addplot[black, line width=1.2pt] coordinates {("
        + f"{lo:.6e}"
        + ","
        + f"{lo:.6e}"
        + ")("
        + f"{hi:.6e}"
        + ","
        + f"{hi:.6e}"
        + ")};"
        "\n"
        r"\addlegendentry{Bisector}"
        "\n"
    )

    # Factors with distinguishable styles + legend entry only once per factor
    for i, f in enumerate(factors):
        g01 = factor_grey[f]
        gpgf = pgf_gray(g01)
        f_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        ls = _factor_pgf_linestyle(i)

        # y = f x (legend)
        tex += (
            r"\addplot[" + gpgf + ", " + ls + r", line width=0.9pt] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(f * lo):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(f * hi):.6e}"
            + ")};"
            "\n"
            r"\addlegendentry{" + f_label + r"}"
            "\n"
        )
        # y = x/f (no legend)
        tex += (
            r"\addplot[" + gpgf + ", " + ls + r", line width=0.9pt, forget plot] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(lo / f):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(hi / f):.6e}"
            + ")};"
            "\n"
        )

    # Safe / Not-safe labels
    tex += (
        r"\node[rotate=45, anchor=south east] at (axis cs:"
        + f"{x_safe:.6e}"
        + ","
        + f"{y_safe:.6e}"
        + r") {Safe};"
        "\n"
        r"\node[rotate=45, anchor=north west] at (axis cs:"
        + f"{x_not:.6e}"
        + ","
        + f"{y_not:.6e}"
        + r") {Not-safe};"
        "\n"
    )

    tex += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(tex.strip(), encoding="utf-8")


def plot_dp_point_vs_fit(
    df_mat: pd.DataFrame,
    out_png: Path,
    out_tex: Path,
    data_dir: Path,
    material_label: str,
    method_styles: Dict[str, Dict[str, str]],
    scatter_alpha: float,
    factors: List[float],
    show_method_legends: bool,
) -> None:
    """DP parity plot for DP_DIFF mode.

    Plots:
      y = DP_point (damage parameter from model evaluation)
      x = DP_fit   (calibration best-fit evaluated at Nf_exp)

    Bisector + factor bands are meaningful again.
    """
    if "DP_fit" not in df_mat.columns or "CP_value" not in df_mat.columns:
        print(f"[WARN] Missing DP_fit/CP_value for DP parity: {material_label}")
        return

    # Match paper style (same as Nf parity scatter)
    scatter_alpha = 0.2
    factors = sorted(set(float(f) for f in factors if np.isfinite(f) and f > 1.0))

    methods = sorted(df_mat["Method"].dropna().unique().tolist())

    xall = pd.to_numeric(df_mat["DP_fit"], errors="coerce").to_numpy(dtype=float)
    yall = pd.to_numeric(df_mat["CP_value"], errors="coerce").to_numpy(dtype=float)
    m = np.isfinite(xall) & np.isfinite(yall) & (xall > 0) & (yall > 0)
    if not np.any(m):
        print(f"[WARN] No positive finite DP data for parity: {material_label}")
        return

    xall = xall[m]
    yall = yall[m]

    lo_raw = float(min(xall.min(), yall.min()))
    hi_raw = float(max(xall.max(), yall.max()))

    logmin = float(np.log10(lo_raw))
    logmax = float(np.log10(hi_raw))
    span = (logmax - logmin) if logmax > logmin else 1.0
    pad = 0.15 * span + 0.05
    lo = 10 ** (logmin - pad)
    hi = 10 ** (logmax + pad)

    # Safe / Not-safe labels (DP-point higher than fit -> conservative)
    x_safe = 10 ** (logmin + 0.20 * span)
    y_safe = 10 ** (logmin + 0.80 * span)
    x_not = 10 ** (logmin + 0.80 * span)
    y_not = 10 ** (logmin + 0.20 * span)

    lo = float(min(lo, x_safe, y_safe, x_not, y_not))
    hi = float(max(hi, x_safe, y_safe, x_not, y_not))

    # Greys for factors
    grey_levels = [0.45, 0.65, 0.80, 0.30, 0.55]
    factor_grey: Dict[float, float] = {f: grey_levels[i % len(grey_levels)] for i, f in enumerate(factors)}

    # --- Matplotlib ---
    fig, ax = plt.subplots(figsize=(7, 7))

    for meth in methods:
        sub = df_mat[df_mat["Method"] == meth]
        x = pd.to_numeric(sub["DP_fit"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["CP_value"], errors="coerce").to_numpy(dtype=float)
        mm = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if np.any(mm):
            st = method_styles.get(meth)
            col = ("#" + st["hex"]) if st else None
            label = meth if show_method_legends else "_nolegend_"
            mk = pgf_mark_to_mpl(st["mark"]) if st else "o"
            is_ext = bool(st and st.get("is_ext") == "1")
            fc = "none" if is_ext or not col else col
            ax.scatter(
                x[mm],
                y[mm],
                label=label,
                alpha=0.2,
                s=30,
                marker=mk,
                facecolors=fc,
                edgecolors=col,
                linewidths=0.9,
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    xs = np.array([lo, hi], dtype=float)
    ax.plot(xs, xs, color="black", linewidth=1.8, label="Bisector")

    for i, f in enumerate(factors):
        g = str(factor_grey[f])
        f_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        ls = _factor_mpl_linestyle(i)
        ax.plot(xs, f * xs, linestyle=ls, color=g, linewidth=1.3, label=f_label)
        ax.plot(xs, xs / f, linestyle=ls, color=g, linewidth=1.3, label="_nolegend_")

    ax.text(x_safe, y_safe, "Safe", rotation=45, color="black", fontsize=11, ha="left", va="top")
    ax.text(x_not, y_not, "Not-safe", rotation=45, color="black", fontsize=11, ha="right", va="bottom")

    ax.set_xlabel(r"Calibration best-fit damage parameter - $DP_{fit}$")
    ax.set_ylabel(r"Damage parameter from loading conditions - $DP$")
    ax.legend(fontsize=9, loc="lower left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # --- TikZ export ---
    scat_files = []
    for meth in methods:
        sub = df_mat[df_mat["Method"] == meth]
        x = pd.to_numeric(sub["DP_fit"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["CP_value"], errors="coerce").to_numpy(dtype=float)
        mm = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if np.any(mm):
            dat = data_dir / "scatter" / f"{_safe_name(meth)}.dat"
            _write_dat_xy(dat, x[mm], y[mm])
            scat_files.append((meth, dat))

    def pgf_gray(level01: float) -> str:
        p = int(round(level01 * 100))
        return f"gray!{p}"

    tex = (r"""\begin{tikzpicture}
\input{Results/Comparison/faber_plot_styles.tex}
\begin{axis}[
  width=0.5\textwidth,
  height=0.5\textwidth,
  xmode=log,
  ymode=log,
  xmin="""
        + f"{lo:.6e}"
        + r""", xmax="""
        + f"{hi:.6e}"
        + r""",
  ymin="""
        + f"{lo:.6e}"
        + r""", ymax="""
        + f"{hi:.6e}"
        + r""",
  axis equal image,
  grid=both,
  xlabel={Calibration best-fit damage parameter - $DP_{fit}$},
  ylabel={Damage parameter from loading conditions - $DP$},
  legend style={at={(0.02,0.02)}, anchor=south west, font=\tiny},
	]
	"""
	    + r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
	    + _latex_escape(material_label)
	    + r"};"
	    + "\n"
	)

    forget = "" if show_method_legends else ", forget plot"

    for meth, dat in scat_files:
        st = method_styles.get(meth)
        mid = st["id"] if st else None
        style = f"faber/method_scatter/{mid}" if mid else "only marks"
        tex += (
            r"\addplot["
            + style
            + forget
            + r", opacity=0.2, fill opacity=0.2"
            + r"] table[x index=0, y index=1] {"
            + dat.as_posix()
            + r"};"
            + "\n"
        )
        if show_method_legends:
            tex += r"\addlegendentry{" + _latex_escape(meth) + r"}" + "\n"

    # Bisector
    tex += (
        r"\addplot[black, line width=1.2pt] coordinates {("
        + f"{lo:.6e}"
        + ","
        + f"{lo:.6e}"
        + ")("
        + f"{hi:.6e}"
        + ","
        + f"{hi:.6e}"
        + ")};"
        "\n"
        r"\addlegendentry{Bisector}"
        "\n"
    )

    # Factors
    for i, f in enumerate(factors):
        g01 = factor_grey[f]
        gpgf = pgf_gray(g01)
        f_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        ls = _factor_pgf_linestyle(i)

        # NOTE: pgfplots expects a single brace group for coordinates.
        # Avoid writing double braces (e.g. `coordinates {{...}}`) which can trigger
        # `Argument of \\pgfplots@addplotimpl@coordinates has an extra }`.
        tex += (
            r"\addplot[" + gpgf + ", " + ls + r", line width=0.9pt] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(f * lo):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(f * hi):.6e}"
            + ")};"
            "\n"
            r"\addlegendentry{" + f_label + r"}"
            "\n"
        )
        tex += (
            r"\addplot[" + gpgf + ", " + ls + r", line width=0.9pt, forget plot] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(lo / f):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(hi / f):.6e}"
            + ")};"
            "\n"
        )

    # Safe / Not-safe labels
    tex += (
        r"\node[rotate=45, anchor=north west] at (axis cs:"
        + f"{x_safe:.6e}"
        + ","
        + f"{y_safe:.6e}"
        + r") {Safe};"
        "\n"
        r"\node[rotate=45, anchor=south east] at (axis cs:"
        + f"{x_not:.6e}"
        + ","
        + f"{y_not:.6e}"
        + r") {Not-safe};"
        "\n"
    )

    tex += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(tex.strip(), encoding="utf-8")


def compute_stats_table(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-(material, method) error statistics.

    The error metric depends on the CSV flag Metric_mode:
      - Metric_mode != 'DP_DIFF'  -> Error_ln
      - Metric_mode == 'DP_DIFF'  -> Error_ln_dp
    """
    rows = []
    if "Metric_mode" not in df.columns:
        df = df.copy()
        df["Metric_mode"] = "Nf"

    for (mat, meth, mode), g in df.groupby(["Material", "Method", "Metric_mode"], dropna=False):
        mode_s = str(mode or "Nf")
        if mode_s.strip().upper() == "DP_DIFF":
            metric = "Error_ln_dp"
            err = pd.to_numeric(g.get("Error_ln_dp"), errors="coerce").to_numpy(dtype=float)
        else:
            metric = "Error_ln"
            err = pd.to_numeric(g.get("Error_ln"), errors="coerce").to_numpy(dtype=float)

        err = err[np.isfinite(err)]
        rows.append(
            {
                "Material": mat,
                "Method": meth,
                "Metric": metric,
                "n": int(err.size),
                "Error_mean": float(np.mean(err)) if err.size else np.nan,
                "Error_std": float(np.std(err, ddof=1)) if err.size > 1 else np.nan,
            }
        )

    return pd.DataFrame(rows).sort_values(["Material", "Method", "Metric"])


def write_stats_table_with_styles(
    stats: pd.DataFrame,
    method_styles: Dict[str, Dict[str, str]],
    out_csv: Path,
    out_tex: Path,
) -> None:
    """Write a single global LaTeX table (independent from material folders).

    Contents:
      - associated color + marker per method (consistent with plots)
      - mean/std of the active error metric per (material, method)
        * Metric='Error_ln' -> Error_ln = ln(Nf / Nf,e)
        * Metric='Error_ln_dp' -> Error_ln_dp = ln(DP / DP_{fit})
    """
    stats = stats.copy()
    stats["Style_id"] = stats["Method"].map(lambda m: method_styles.get(m, {}).get("id", ""))
    stats["Color_hex"] = stats["Method"].map(lambda m: method_styles.get(m, {}).get("hex", ""))
    stats["PGF_color"] = stats["Method"].map(lambda m: method_styles.get(m, {}).get("color", ""))
    stats["PGF_mark"] = stats["Method"].map(lambda m: method_styles.get(m, {}).get("mark", ""))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    stats.to_csv(out_csv, index=False)

    # LaTeX table
    out_tex.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines += [
        "% Auto-generated by 2_compare_results.py",
        "% Table with: method styles + per-material error stats.",
        "% Requires in LaTeX preamble:",
        "%   \\usepackage{booktabs}",
        "%   \\usepackage{tabularx}",
        "%   % $\\mu$ and $\\sigma$ are computed for the active Error parameter (Error$_{\\ln}$ or Error$_{\\ln,dp}$).",
        "%   \\usepackage{tikz}",
        "%   \\usepackage{pgfplots}",
        "%   \\usepackage{xcolor}",
        "",
        r"\input{Results/Comparison/faber_plot_styles.tex}",
        "",
        r"\begin{tabularx}{\textwidth}{lXlccrr}",
        r"\toprule",
        r"Material & Method & Metric & Color & Mark & $\mu_{Error}$ & $\sigma_{Error}$\\",
        r"\midrule",
    ]

    def fmt(x: float) -> str:
        if not np.isfinite(x):
            return "--"
        # compact but stable (avoid scientific for typical values)
        return f"{x:.4g}"

    for _, r in stats.sort_values(["Method", "Material", "Metric"]).iterrows():
        mat = _latex_escape(r["Material"])
        meth = _latex_escape(r["Method"])
        meth_raw = str(r["Method"])
        is_ext = meth_raw.endswith("_ext")
        metric_raw = str(r.get("Metric", "Error_ln"))
        if metric_raw == "Error_ln_dp":
            metric_cell = r"$Error_{\ln,dp}$"
        else:
            metric_cell = r"$Error_{\ln}$"
        col = r["PGF_color"]
        mark = r["PGF_mark"]
        mu = fmt(float(r["Error_mean"])) if pd.notna(r["Error_mean"]) else "--"
        sig = fmt(float(r["Error_std"])) if pd.notna(r["Error_std"]) else "--"

        # Table swatches:
        #   - Color column: show ONLY line style (solid vs densely dashed) with the method color.
        #   - Mark column: show ONLY the marker (filled for base, white-filled for _ext) with no line segment.
        line_style = "densely dashed" if is_ext else "solid"
        color_swatch = (
            r"\tikz[baseline=-0.6ex]{\draw[" + col + r"," + line_style + r", line width=0.9pt] (0,0) -- (0.5,0);}"  # noqa: E501
            if isinstance(col, str) and col
            else "--"
        )
        marker_swatch = (
            r"\tikz[baseline=-0.6ex]{"
            + r"\draw[" + col
            + r", mark=" + mark
            + r", mark options={solid, fill=" + ("white" if is_ext else col) + r"}] plot coordinates {(0,0)};"
            + r"}"
            if isinstance(col, str) and col and isinstance(mark, str) and mark
            else "--"
        )

        mu_cell = (rf"$ {mu} $") if mu != "--" else "--"
        sig_cell = (rf"$ {sig} $") if sig != "--" else "--"
        lines.append(f"{mat} & {meth} & {metric_cell} & {color_swatch} & {marker_swatch} & {mu_cell} & {sig_cell}\\\\")

    lines += [r"\bottomrule", r"\end{tabularx}", ""]
    out_tex.write_text("\n".join(lines), encoding="utf-8")


# =========================
# Main
# =========================
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(allow_abbrev=False)  # IMPORTANT for Jupyter (-f kernel.json)
    ap.add_argument("--results_dir", type=str, default="Results")
    ap.add_argument("--out_dir", type=str, default="Results/Comparison")
    ap.add_argument("--flag_filter", type=str, default="all", choices=["all", "design", "other"])

    # scatter opacity (clustering visibility)
    ap.add_argument("--scatter_alpha", type=float, default=0.18)

    # normal curves aesthetics
    ap.add_argument("--pdf_lw", type=float, default=2.0)
    ap.add_argument("--pdf_alpha", type=float, default=0.85)

    # factor bands (default: 3 and 5)
    ap.add_argument("--factors", type=float, nargs="*", default=[3.0, 5.0])

    # Legends (methods) can be disabled because you now get a global table
    ap.add_argument("--show_method_legends", action="store_true", help="Include method legend entries in plots.")

    args, _unknown = ap.parse_known_args(argv)

    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scatter_alpha = float(np.clip(args.scatter_alpha, 0.0, 1.0))
    pdf_alpha = float(np.clip(args.pdf_alpha, 0.0, 1.0))
    pdf_lw = float(max(args.pdf_lw, 0.1))
    factors = [float(f) for f in (args.factors or []) if np.isfinite(f) and f > 1.0]

    df = collect_all_results(results_dir)

    # Minimal required columns (support both Nf-based and DP-based modes)
    for needed in ["Material", "Method", "Nf_exp"]:
        if needed not in df.columns:
            raise RuntimeError(f"Missing required column '{needed}' after reading CSVs.")
    if "CP_value" not in df.columns:
        raise RuntimeError("Missing required column 'CP_value' (damage parameter) after reading CSVs.")

    if "Metric_mode" not in df.columns:
        df["Metric_mode"] = "Nf"
    else:
        df["Metric_mode"] = df["Metric_mode"].fillna("Nf").astype(str)

    if "Error_ln" not in df.columns:
        df["Error_ln"] = np.nan
        if "Nf_expected" in df.columns:
            m = (df["Nf_exp"] > 0) & (df["Nf_expected"] > 0) & np.isfinite(df["Nf_exp"]) & np.isfinite(df["Nf_expected"])
            df.loc[m, "Error_ln"] = np.log(df.loc[m, "Nf_exp"] / df.loc[m, "Nf_expected"])

    # DP log-ratio error (used when Metric_mode == 'DP_DIFF')
    if "Error_ln_dp" not in df.columns:
        df["Error_ln_dp"] = np.nan
    if "DP_fit" in df.columns:
        m = np.isfinite(df["DP_fit"]) & np.isfinite(df["CP_value"]) & (df["DP_fit"] > 0.0) & (df["CP_value"] > 0.0)
        df.loc[m, "Error_ln_dp"] = np.log(df.loc[m, "CP_value"] / df.loc[m, "DP_fit"])

    if "Error_dp" not in df.columns:
        df["Error_dp"] = np.nan
    if "DP_fit" in df.columns:
        m = np.isfinite(df["DP_fit"]) & np.isfinite(df["CP_value"])
        df.loc[m, "Error_dp"] = df.loc[m, "CP_value"] - df.loc[m, "DP_fit"]

    if args.flag_filter != "all" and "Flag" in df.columns:
        if args.flag_filter == "design":
            df = df[df["Flag"] == 0]
        elif args.flag_filter == "other":
            df = df[df["Flag"] == 1]

    # Build global palette:
    # - colors per Method_family
    # - markers per Method
    method_styles, family_styles = build_style_maps(df)

    # Write global style file (used by every plot)
    write_faber_style_file(method_styles, family_styles, out_dir / "faber_plot_styles.tex")

    # Global stats (for your single "legend table")
    stats = compute_stats_table(df)
    stats.to_csv(out_dir / "summary_error_stats.csv", index=False)
    write_stats_table_with_styles(
        stats=stats,
        method_styles=method_styles,
        out_csv=out_dir / "summary_error_stats_with_styles.csv",
        out_tex=out_dir / "summary_error_stats_table.tex",
    )

    materials = sorted(df["Material"].dropna().unique().tolist())
    if not materials:
        raise RuntimeError("No materials found in loaded CSVs.")

    for mat in materials:
        df_mat = df[df["Material"] == mat].copy()
        if df_mat.empty:
            continue

        # Split by metric mode
        mode_u = df_mat["Metric_mode"].astype(str).str.strip().str.upper()
        df_nf = df_mat[mode_u != "DP_DIFF"].copy()
        df_dp = df_mat[mode_u == "DP_DIFF"].copy()

        safe_mat = _safe_name(mat)
        mat_dir = out_dir / safe_mat
        mat_dir.mkdir(parents=True, exist_ok=True)

        tikz_dir = mat_dir / "tikz_data"
        tikz_dir.mkdir(parents=True, exist_ok=True)
        tikz_nf = tikz_dir / "nf"
        tikz_dp = tikz_dir / "dp"
        tikz_nf.mkdir(parents=True, exist_ok=True)
        tikz_dp.mkdir(parents=True, exist_ok=True)

        # -------------------------
        # Nf-based plots (standard)
        # -------------------------
        if not df_nf.empty:
            plot_boxplot_errors(
                df_mat=df_nf,
                out_png=mat_dir / f"{safe_mat}_boxplot_error.png",
                out_tex=mat_dir / f"{safe_mat}_boxplot_error.tex",
                data_dir=tikz_nf,
                material_label=mat,
                method_styles=method_styles,
                error_col="Error_ln",
                ylabel=r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$",
            )

            plot_normal_curves(
                df_mat=df_nf,
                out_png=mat_dir / f"{safe_mat}_normal_error_pdf.png",
                out_tex=mat_dir / f"{safe_mat}_normal_error_pdf.tex",
                data_dir=tikz_nf,
                material_label=mat,
                method_styles=method_styles,
                pdf_lw=pdf_lw,
                pdf_alpha=pdf_alpha,
                show_method_legends=args.show_method_legends,
                error_col="Error_ln",
                xlabel=r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$",
            )

            plot_expected_vs_experimental(
                df_mat=df_nf,
                out_png=mat_dir / f"{safe_mat}_expected_vs_experimental.png",
                out_tex=mat_dir / f"{safe_mat}_expected_vs_experimental.tex",
                data_dir=tikz_nf,
                material_label=mat,
                method_styles=method_styles,
                scatter_alpha=scatter_alpha,
                factors=factors,
                show_method_legends=args.show_method_legends,
            )

        # -------------------------
        # DP-based plots (DP_DIFF)
        # -------------------------
        if not df_dp.empty:
            plot_boxplot_errors(
                df_mat=df_dp,
                out_png=mat_dir / f"{safe_mat}_dp_boxplot_error.png",
                out_tex=mat_dir / f"{safe_mat}_dp_boxplot_error.tex",
                data_dir=tikz_dp,
                material_label=mat,
                method_styles=method_styles,
                error_col="Error_ln_dp",
                ylabel=r"$Error_{\ln,dp}=\ln\left(\frac{DP}{DP_{fit}}\right)$",
            )

            plot_normal_curves(
                df_mat=df_dp,
                out_png=mat_dir / f"{safe_mat}_dp_error_pdf.png",
                out_tex=mat_dir / f"{safe_mat}_dp_error_pdf.tex",
                data_dir=tikz_dp,
                material_label=mat,
                method_styles=method_styles,
                pdf_lw=pdf_lw,
                pdf_alpha=pdf_alpha,
                show_method_legends=args.show_method_legends,
                error_col="Error_ln_dp",
                xlabel=r"$Error_{\ln,dp}=\ln\left(\frac{DP}{DP_{fit}}\right)$",
            )

            plot_dp_point_vs_fit(
                df_mat=df_dp,
                out_png=mat_dir / f"{safe_mat}_dp_point_vs_fit.png",
                out_tex=mat_dir / f"{safe_mat}_dp_point_vs_fit.tex",
                data_dir=tikz_dp,
                material_label=mat,
                method_styles=method_styles,
                scatter_alpha=scatter_alpha,
                factors=factors,
                show_method_legends=args.show_method_legends,
            )

    print(f"Done. Outputs written to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()