#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from multifat.eval.metrics import log10_dp_ratio
from multifat.models.references import DISABLED_METHODS, IMPLEMENTATION_REVISION
from multifat.models.comparison import DEFAULT_CAPABILITIES, capabilities_for_name
import numpy as np
import pandas as pd


# =========================
# Utilities
# =========================
def _fs_safe_stem(text: str) -> str:
    """Return a name that is safe on disk (underscores are kept)."""
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in str(text))


_TEX_ESCAPES = {
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


def _tex_escape(text: str) -> str:
    """Escape a string so that LaTeX prints it literally in text mode."""
    text = str(text)
    return "".join(_TEX_ESCAPES.get(ch, ch) for ch in text)


def _tex_identifier(method: str) -> str:
    """
    Derive an identifier for a method that is safe inside pgfplots keys:
    alphanumeric only, starting with a letter, and reproducible from run to run
    (duplicates are resolved by the caller).
    """
    stem = re.sub(r"[^A-Za-z0-9]+", "", str(method))
    if not stem:
        stem = "method"
    if stem[0].isdigit():
        stem = "m" + stem
    return stem[:32]


def _gaussian_density(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    if not np.isfinite(mu) or not np.isfinite(sigma) or sigma <= 0:
        return np.full_like(x, np.nan, dtype=float)
    return (1.0 / (sigma * math.sqrt(2.0 * math.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _locate_column(frame: pd.DataFrame, names: List[str]) -> Optional[str]:
    normalized = {col: col.strip().lower().replace(" ", "_") for col in frame.columns}
    lookup = {norm: col for col, norm in normalized.items()}
    for name in names:
        probe = name.strip().lower().replace(" ", "_")
        if probe in lookup:
            return lookup[probe]
    return None


def _derive_labels_from_stem(source: Path) -> Tuple[Optional[str], Optional[str]]:
    stem = source.stem
    chunks = stem.split("_")
    if len(chunks) < 2:
        return None, None

    def _is_material_like(token: str) -> bool:
        return any(ch.isdigit() for ch in token) and any(ch.isalpha() for ch in token)

    head = "_".join(chunks[:-1])
    tail = chunks[-1]
    if _is_material_like(tail):
        return tail, head
    if _is_material_like(chunks[0]):
        return chunks[0], "_".join(chunks[1:])
    return None, None


# =========================
# Styling (ColorBrewer-like palette, shared by every figure)
# =========================
# ColorBrewer "Paired" (12 entries) -- strong separation and print friendly
_PAIRED_PALETTE_12 = [
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

# Usual pgfplots markers (cycled when there are more methods than markers)
_PGF_MARKERS = [
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


def _dedupe_key(base: str, taken: set[str]) -> str:
    """Produce a deterministic key that is unique case-insensitively, starting from a base id."""
    candidate = base
    counter = 1
    while candidate.lower() in taken:
        candidate = f"{base}{counter}"
        counter += 1
    taken.add(candidate.lower())
    return candidate


def pgf_mark_to_mpl(mark: str) -> str:
    """Translate a pgfplots marker name into the closest Matplotlib marker."""
    token = str(mark or "").strip()
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
    }.get(token, "o")


def _colormap_for_family(family: str) -> str:
    """
    Pick a ColorBrewer-like sequential colormap for a method family.
    Matplotlib's default sequential maps (Reds, Blues, Greens, ...) already follow ColorBrewer.
    """
    key = (family or "UNKNOWN").lower()
    # Main requirement:
    # - CP -> red shades
    # - invariants -> blue shades
    if "cp" in key or "critical" in key:
        return "Reds"
    if "invariant" in key:
        return "Blues"
    # Reasonable fallbacks for the remaining families
    if "energy" in key:
        return "Greens"
    if "integral" in key:
        return "Purples"
    if "stress" in key:
        return "Purples"
    if "strain" in key:
        return "Oranges"
    if "unknown" in key:
        return "Greys"
    return "Greys"


def _colormap_hex_at(palette: str, idx: int, count: int) -> str:
    """
    Sample a sequential colormap to obtain a legible, print-friendly colour.
    The lightest end is avoided so that lines stay visible on white backgrounds.
    """
    cmap = plt.get_cmap(palette)
    if count <= 1:
        position = 0.70
    else:
        # keep colours away from both extreme ends of the ramp
        samples = np.linspace(0.40, 0.90, count)
        position = float(samples[idx])
    rgb = cmap(position)
    # mcolors.to_hex -> '#RRGGBB'
    return mcolors.to_hex(rgb, keep_alpha=False).lstrip("#").upper()


def build_style_maps(
    df: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    """
    Build reproducible style mappings driven by Method_family (the method family):

      family_styles: family -> {id, cmap}
      method_styles: method -> {id, family, base, is_ext, hex, mark, color}

    RULES:
      - Every FAMILY owns a sequential palette (e.g. CP=Reds, invariants=Blues).
      - Every BASE METHOD (the method name with an optional "_ext" suffix removed)
        inside a family receives a distinct shade.
      - When a method has an "_ext" variant, both share the SAME colour and marker,
        and are told apart by line style / marker fill:
          * base: solid line, filled marker
          * _ext: densely dashed line, hollow marker
      - The mapping is reproducible: sorted families plus sorted base methods per family.
    """
    if "Method" not in df.columns:
        raise ValueError("DataFrame must contain a 'Method' column.")
    if "Method_family" not in df.columns:
        df = df.copy()
        df["Method_family"] = "UNKNOWN"
    df["Method_family"] = df["Method_family"].fillna("UNKNOWN").astype(str)

    def _strip_ext(method: str) -> str:
        return method[:-4] if isinstance(method, str) and method.endswith("_ext") else str(method)

    def _has_ext(method: str) -> bool:
        return isinstance(method, str) and method.endswith("_ext")

    # --- Families -> colormap ---
    ordered_families = sorted(
        [f for f in df["Method_family"].dropna().unique().tolist() if isinstance(f, str)]
    )
    family_ids_taken: set[str] = set()
    family_styles: Dict[str, Dict[str, str]] = {}
    for family in ordered_families:
        stem = _tex_identifier(family)
        fam_id = _dedupe_key(stem, family_ids_taken)
        family_styles[family] = {"id": fam_id, "cmap": _colormap_for_family(family)}

    if "UNKNOWN" not in family_styles:
        fam_id = _dedupe_key(_tex_identifier("UNKNOWN"), family_ids_taken)
        family_styles["UNKNOWN"] = {
            "id": fam_id,
            "cmap": _colormap_for_family("UNKNOWN"),
        }

    # --- Methods -> per-family shades (keyed by BASE method) + per-base marker ---
    ordered_methods = sorted(
        [m for m in df["Method"].dropna().unique().tolist() if isinstance(m, str)]
    )

    # method -> family (first non-null value)
    method_family_of: Dict[str, str] = {}
    for method in ordered_methods:
        family_values = df.loc[df["Method"] == method, "Method_family"].dropna()
        family = str(family_values.iloc[0]) if len(family_values) else "UNKNOWN"
        if family not in family_styles:
            family = "UNKNOWN"
        method_family_of[method] = family

    # family -> base methods
    bases_by_family: Dict[str, List[str]] = {}
    for method in ordered_methods:
        family = method_family_of.get(method, "UNKNOWN")
        base_name = _strip_ext(method)
        bases_by_family.setdefault(family, []).append(base_name)

    for family in bases_by_family:
        bases_by_family[family] = sorted(set(bases_by_family[family]))

    # assign a shade to each base inside its family
    hex_by_base: Dict[Tuple[str, str], str] = {}
    for family, bases in bases_by_family.items():
        palette = family_styles.get(family, family_styles["UNKNOWN"]).get("cmap", "Greys")
        for shade_index, base_name in enumerate(bases):
            hex_by_base[(family, base_name)] = _colormap_hex_at(palette, shade_index, len(bases))

    # assign markers per base, reproducibly, across ALL bases
    ordered_all_bases: List[str] = sorted({_strip_ext(method) for method in ordered_methods})
    mark_by_base: Dict[str, str] = {}
    for position, base_name in enumerate(ordered_all_bases):
        mark_by_base[base_name] = _PGF_MARKERS[position % len(_PGF_MARKERS)]

    method_ids_taken: set[str] = set()
    method_styles: Dict[str, Dict[str, str]] = {}
    for method in ordered_methods:
        base_name = _strip_ext(method)
        family = method_family_of.get(method, "UNKNOWN")
        hex_value = hex_by_base.get((family, base_name), _colormap_hex_at("Greys", 0, 1))
        marker = mark_by_base.get(base_name, _PGF_MARKERS[0])

        key = _dedupe_key(_tex_identifier(method), method_ids_taken)
        color_name = f"faberMethodColor{key}"
        method_styles[method] = {
            "id": key,
            "family": family,
            "base": base_name,
            "is_ext": "1" if _has_ext(method) else "0",
            "hex": hex_value,
            "mark": marker,
            "color": color_name,
        }

    return method_styles, family_styles


def write_faber_style_file(
    method_styles: Dict[str, Dict[str, str]],
    family_styles: Dict[str, Dict[str, str]],
    out_path: Path,
) -> None:
    """
    Emit a LaTeX file holding the colour and pgfplots styles.

    Behaviour:
      - Colours are declared per METHOD (the shade follows its method family palette).
      - When a method has an "_ext" variant, base and "_ext" share the SAME colour and marker,
        and are told apart by:
          * base: solid line, filled marker
          * _ext: densely dashed line, hollow marker
      - Boxplots: outlined boxes only; outliers are always a black 'x'.
      - The global font size is \\scriptsize.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    out_lines: List[str] = []
    out_lines += [
        r"\providecommand{\faberPlotStylesLoaded}{0}",
        r"\ifnum\faberPlotStylesLoaded=0",
        r"  \renewcommand{\faberPlotStylesLoaded}{1}",
        r"  \tikzset{every picture/.style={font=\scriptsize}}",
        "",
    ]

    # Declare one colour per method
    for method in sorted(method_styles.keys()):
        style = method_styles[method]
        out_lines.append(rf"  \definecolor{{{style['color']}}}{{HTML}}{{{style['hex']}}}")

    out_lines += [
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

    # Per-method styles
    for method in sorted(method_styles.keys()):
        style = method_styles[method]
        key = style["id"]
        colour = style["color"]
        marker = style["mark"]
        extended = style.get("is_ext", "0") == "1"

        # Base colour style
        out_lines.append(rf"    faber/method/{key}/.style={{color={colour}}},")

        # Scatter style:
        #  - base: filled marker
        #  - _ext: hollow marker
        fill_option = colour if not extended else "white"
        out_lines.append(
            rf"    faber/method_scatter/{key}/.style={{faber/method/{key}, only marks, mark={marker}, mark size=1.6pt, "
            rf"draw={colour}, mark options={{solid, fill={fill_option}}}}},"  # noqa: E501
        )

        # Line style:
        #  - base: solid
        #  - _ext: densely dashed
        dash = "densely dashed" if extended else "solid"
        out_lines.append(
            rf"    faber/method_line/{key}/.style={{faber/method/{key}, {dash}, line width=1.1pt}},"
        )

        # Boxplot style (outline only)
        out_lines.append(
            rf"    faber/method_box/{key}/.style={{faber/method/{key}, boxplot, boxplot/draw direction=y, "
            rf"draw={colour}, fill=none}},"
        )

    out_lines += [
        r"  }",
        r"\fi",
        "",
    ]

    out_path.write_text("\n".join(out_lines), encoding="utf-8")


# =========================
# CSV loading
# =========================
def read_results_csv(csv_path: Path) -> Optional[pd.DataFrame]:
    try:
        frame = pd.read_csv(csv_path)
    except Exception as exc:
        print(f"[WARN] Unable to read {csv_path}: {exc}")
        return None

    if frame is None or frame.empty:
        return None

    frame.columns = [c.strip() for c in frame.columns]

    idx_material = _locate_column(frame, ["Material", "material", "mat"])
    idx_method = _locate_column(frame, ["Method", "method", "model"])
    idx_family = _locate_column(
        frame,
        [
            "Method_family",
            "MethodFamily",
            "Family",
            "Method_type",
            "MethodType",
            "Type",
        ],
    )
    idx_metric = _locate_column(
        frame, ["Metric_mode", "Metric", "Mode", "Error_metric", "ErrorMetric"]
    )
    idx_nfexp = _locate_column(frame, ["Nf_exp", "Nexp", "Nfexperimental", "Nf_experimental"])
    idx_nfpred = _locate_column(frame, ["Nf_expected", "Npred", "Nfexpected", "Nf_expected_cycles"])
    idx_cp = _locate_column(
        frame, ["CP_value", "CP", "DamageParameter", "Damage_parameter", "DP", "dp"]
    )
    idx_dpfit = _locate_column(frame, ["DP_fit", "DPfit", "CP_fit", "DP_fit_value"])
    idx_err = _locate_column(frame, ["Error_ln", "error_ln", "Error", "error"])
    idx_errlndp = _locate_column(
        frame, ["Error_ln_dp", "ErrorLnDP", "Error_lnDP", "Error_lnDP", "Error_dp_ln"]
    )
    # Deprecated/legacy (difference-based) DP error
    idx_errdp = _locate_column(frame, ["Error_dp", "error_dp", "ErrorDP", "Error_DP"])
    idx_flag = _locate_column(frame, ["Flag", "flag"])

    # Recover missing material/method from the file name
    if idx_material is None or idx_method is None:
        material_guess, method_guess = _derive_labels_from_stem(csv_path)
        if idx_material is None:
            frame["Material"] = material_guess if material_guess is not None else "UNKNOWN"
            idx_material = "Material"
        if idx_method is None:
            frame["Method"] = method_guess if method_guess is not None else "UNKNOWN"
            idx_method = "Method"

    # Rename to the canonical column names
    renames = {}
    if idx_material and idx_material != "Material":
        renames[idx_material] = "Material"
    if idx_method and idx_method != "Method":
        renames[idx_method] = "Method"
    if idx_family and idx_family != "Method_family":
        renames[idx_family] = "Method_family"
    if idx_metric and idx_metric != "Metric_mode":
        renames[idx_metric] = "Metric_mode"
    if idx_nfexp and idx_nfexp != "Nf_exp":
        renames[idx_nfexp] = "Nf_exp"
    if idx_nfpred and idx_nfpred != "Nf_expected":
        renames[idx_nfpred] = "Nf_expected"
    if idx_cp and idx_cp != "CP_value":
        renames[idx_cp] = "CP_value"
    if idx_dpfit and idx_dpfit != "DP_fit":
        renames[idx_dpfit] = "DP_fit"
    if idx_err and idx_err != "Error_ln":
        renames[idx_err] = "Error_ln"
    if idx_errlndp and idx_errlndp != "Error_ln_dp":
        renames[idx_errlndp] = "Error_ln_dp"
    if idx_errdp and idx_errdp != "Error_dp":
        renames[idx_errdp] = "Error_dp"
    if idx_flag and idx_flag != "Flag":
        renames[idx_flag] = "Flag"

    frame = frame.rename(columns=renames)

    # Make sure a method family column is available
    if "Method_family" not in frame.columns:
        frame["Method_family"] = "UNKNOWN"
    else:
        frame["Method_family"] = frame["Method_family"].fillna("UNKNOWN").astype(str)

    # Make sure a metric mode column is available
    if "Metric_mode" not in frame.columns:
        frame["Metric_mode"] = "Nf"
    else:
        frame["Metric_mode"] = frame["Metric_mode"].fillna("Nf").astype(str)

    if "Error_log10_dp" not in frame:
        frame["Error_log10_dp"] = np.nan

    # Numeric conversions
    for column in [
        "Nf_exp",
        "Nf_expected",
        "CP_value",
        "DP_fit",
        "Error_ln",
        "Error_ln_dp",
        "Error_log10_dp",
        "Error_dp",
        "Flag",
    ]:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

    # Derive Error_ln when it is missing and Nf_expected is present
    if (
        "Error_ln" not in frame.columns
        and "Nf_exp" in frame.columns
        and "Nf_expected" in frame.columns
    ):
        frame["Error_ln"] = np.nan
        usable = (
            (frame["Nf_exp"] > 0)
            & (frame["Nf_expected"] > 0)
            & np.isfinite(frame["Nf_exp"])
            & np.isfinite(frame["Nf_expected"])
        )
        frame.loc[usable, "Error_ln"] = np.log(
            frame.loc[usable, "Nf_exp"] / frame.loc[usable, "Nf_expected"]
        )

    # Derive the DP log-ratio error when it is missing and DP_fit + CP_value exist
    if "Error_ln_dp" not in frame.columns:
        frame["Error_ln_dp"] = np.nan
    if "DP_fit" in frame.columns and "CP_value" in frame.columns:
        usable = (
            np.isfinite(frame["DP_fit"])
            & np.isfinite(frame["CP_value"])
            & (frame["DP_fit"] > 0.0)
            & (frame["CP_value"] > 0.0)
        )
        frame.loc[usable, "Error_ln_dp"] = np.log(
            frame.loc[usable, "CP_value"] / frame.loc[usable, "DP_fit"]
        )
        frame["Error_log10_dp"] = log10_dp_ratio(frame["DP_fit"], frame["CP_value"])

    # Legacy: only fill the difference-based Error_dp (older CSVs may still carry it)
    if "Error_dp" not in frame.columns:
        frame["Error_dp"] = np.nan
    if "DP_fit" in frame.columns and "CP_value" in frame.columns:
        usable = np.isfinite(frame["DP_fit"]) & np.isfinite(frame["CP_value"])
        frame.loc[usable, "Error_dp"] = frame.loc[usable, "CP_value"] - frame.loc[usable, "DP_fit"]

    frame["source_file"] = str(csv_path)
    return frame


def collect_all_results(results_dir: Path) -> pd.DataFrame:
    csv_paths = sorted(
        [
            path
            for path in results_dir.rglob("*.csv")
            if path.is_file() and not {"Comparison", "FatigueLimit"} & set(path.parts)
        ]
    )
    if not csv_paths:
        raise FileNotFoundError(f"No CSV files located under: {results_dir}")

    frames = []
    for path in csv_paths:
        frame = read_results_csv(path)
        if frame is not None and not frame.empty:
            frames.append(frame)

    if not frames:
        raise RuntimeError(f"CSV files exist under {results_dir}, but none held readable data.")
    return pd.concat(frames, ignore_index=True)


def annotate_comparison_capabilities(frame: pd.DataFrame) -> pd.DataFrame:
    """Identify valid comparison quantities without fitting auxiliary native indices."""
    frame = frame.copy()
    caps_by_name = {}
    for name in frame["Method"].dropna().unique():
        try:
            caps_by_name[name] = capabilities_for_name(name)
        except KeyError:
            caps_by_name[name] = DEFAULT_CAPABILITIES

    prediction = frame["Method"].map(
        lambda name: caps_by_name.get(name, DEFAULT_CAPABILITIES).life_prediction
    )
    dp = frame["Method"].map(
        lambda name: caps_by_name.get(name, DEFAULT_CAPABILITIES).dp_comparison
    )
    mode = frame["Metric_mode"].astype(str).str.strip().str.upper()
    native = prediction == "native"
    frame["Life_prediction_source"] = prediction
    frame["Life_comparison_supported"] = (native & (mode == "NF_DIRECT")) | (
        ~native & (mode != "DP_DIFF") & (mode != "NF_DIRECT")
    )
    dp_fit = pd.to_numeric(
        frame.get("DP_fit", pd.Series(np.nan, index=frame.index)), errors="coerce"
    )
    native_dp_usable = np.isfinite(dp_fit) & (dp_fit > 0)
    frame["DP_comparison_supported"] = dp & ((mode != "NF_DIRECT") | native_dp_usable)
    frame["DP_reference_source"] = frame["Method"].map(
        lambda name: caps_by_name.get(name, DEFAULT_CAPABILITIES).dp_reference
    )
    return frame


def write_native_life_summary(frame: pd.DataFrame, output_root: Path) -> None:
    """Count native infinite-life predictions for failed and censored observations."""
    selected = frame[frame["Life_prediction_source"] == "native"]
    dest = output_root / "native_life_prediction_summary.csv"
    if selected.empty:
        dest.unlink(missing_ok=True)
        return
    rows = []
    for (material, method), group in selected.groupby(["Material", "Method"]):
        native = pd.to_numeric(
            group.get("Nf_native", pd.Series(np.nan, index=group.index)),
            errors="coerce",
        )
        flag = pd.to_numeric(group["Flag"], errors="coerce")
        failures = flag != RUNOUT_FLAG
        failed_inf = int((failures & np.isposinf(native)).sum())
        runout_inf = int(((~failures) & np.isposinf(native)).sum())
        rows.append(
            {
                "Material": material,
                "Method": method,
                "Failure_observations": int(failures.sum()),
                "Failed_specimens_predicted_infinite_life": failed_inf,
                "Runout_observations": int((~failures).sum()),
                "Runouts_predicted_infinite_life": runout_inf,
                "Undefined_native_life": int(native.isna().sum()),
                "Life_source": "Nf_native",
                "Life_comparison_supported": bool(group["Life_comparison_supported"].all()),
                "DP_comparison_supported": False,
            }
        )
    pd.DataFrame(rows).to_csv(dest, index=False)


# =========================
# TikZ helpers
# =========================
def _dump_column(dest: Path, values: np.ndarray) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    with open(dest, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value:.12e}\n")


def _dump_pairs(dest: Path, xs: np.ndarray, ys: np.ndarray) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    usable = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[usable]
    ys = ys[usable]
    with open(dest, "w", encoding="utf-8") as handle:
        for x_val, y_val in zip(xs, ys):
            handle.write(f"{x_val:.12e} {y_val:.12e}\n")


def _tex_header(note: str) -> str:
    return (
        "% Auto-generated by 2_compare_results.py\n"
        f"% {note}\n"
        "% NOTE on paths: data file paths are written with '_' (not '\\_') on purpose.\n"
        "% Requires in LaTeX preamble:\n"
        "%   \\usepackage{xcolor}\n"
        "%   \\usepackage{tikz}\n"
        "%   \\usepackage{pgfplots}\n"
        "%   \\pgfplotsset{compat=1.18}\n"
        "%   \\usepgfplotslibrary{statistics}\n\n"
    )


def _factor_pgf_dash(idx: int) -> str:
    # Keep Factor 3 and Factor 5 visually distinct (densely dashed vs densely dotted)
    dashes = [
        "densely dashed",
        "densely dotted",
        "densely dashdotted",
        "dashdotted",
        "dashed",
        "dotted",
    ]
    return dashes[idx % len(dashes)]


def _factor_mpl_dash(idx: int) -> str:
    dashes = ["--", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 3))]
    return dashes[idx % len(dashes)]


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
    method_names = sorted(df_mat["Method"].dropna().unique().tolist())

    if ylabel is None:
        if error_col == "Error_log10_dp":
            ylabel = r"$\log_{10}\left(\frac{DP_e}{DP}\right)$"
        elif error_col == "Error_ln_dp":
            ylabel = r"$\ln(DP/DP_{fit})$"
        else:
            ylabel = r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$"

    # Store (median, method, values) so the boxes can be ordered from the highest median down.
    records: List[Tuple[float, str, np.ndarray]] = []
    for method in method_names:
        values = df_mat.loc[df_mat["Method"] == method, error_col].dropna().to_numpy(dtype=float)
        if values.size > 0:
            median_value = float(np.nanmedian(values)) if np.any(np.isfinite(values)) else np.nan
            records.append((median_value, method, values))

    if not records:
        print(f"[WARN] No {error_col} available for the boxplot: {material_label}")
        return

    # Ordering: highest median first, NaNs last, ties broken by method name.
    records.sort(
        key=lambda rec: (
            0 if np.isfinite(rec[0]) else 1,
            -rec[0] if np.isfinite(rec[0]) else 0.0,
            rec[1],
        )
    )

    tick_labels = [method for _, method, _ in records]
    samples = [values for _, _, values in records]

    # Matplotlib (colours kept in sync with the TikZ output)
    fig, ax = plt.subplots(figsize=(10, 5))
    outlier_style = dict(
        marker="x",
        markerfacecolor="none",
        markeredgecolor="black",
        markeredgewidth=1.0,
        markersize=5,
        linestyle="none",
    )
    slots = np.arange(1, len(samples) + 1)
    box_art = ax.boxplot(
        samples,
        positions=slots,
        showfliers=True,
        patch_artist=True,
        flierprops=outlier_style,
    )
    ax.set_xticks(slots, tick_labels)
    for position, (box, method) in enumerate(zip(box_art["boxes"], tick_labels)):
        style = method_styles.get(method)
        if style:
            color_hex = "#" + style["hex"]
            # Outline only (no fill)
            box.set_facecolor("none")
            box.set_edgecolor(color_hex)
            box.set_linewidth(1.2)
            # Whiskers, caps and median use the same colour
            for whisker_idx in (2 * position, 2 * position + 1):
                if whisker_idx < len(box_art.get("whiskers", [])):
                    box_art["whiskers"][whisker_idx].set_color(color_hex)
                if whisker_idx < len(box_art.get("caps", [])):
                    box_art["caps"][whisker_idx].set_color(color_hex)
            if position < len(box_art.get("medians", [])):
                box_art["medians"][position].set_color(color_hex)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # TikZ
    data_files = []
    for method in tick_labels:
        values = df_mat.loc[df_mat["Method"] == method, error_col].dropna().to_numpy(dtype=float)
        dat_path = data_dir / "boxplot" / f"{_fs_safe_stem(method)}.dat"
        _dump_column(dat_path, values)
        data_files.append((method, dat_path))

    xticks = ",".join(str(idx + 1) for idx in range(len(data_files)))
    x_tick_labels = ",".join(_tex_escape(method) for method, _ in data_files)

    body = (
        _tex_header(f"Boxplot of {error_col} by method")
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
        + x_tick_labels
        + r"""},
  x tick label style={rotate=30, anchor=east},
]
"""
        + r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
        + _tex_escape(material_label)
        + r"};"
        + "\n"
    )

    for position, (method, dat_path) in enumerate(data_files, start=1):
        style = method_styles.get(method)
        key = style["id"] if style else None
        plot_style = f"faber/method_box/{key}" if key else "boxplot"
        # IMPORTANT: keep '_' in the path (by request, no LaTeX escaping happens here)
        body += (
            rf"\addplot[ {plot_style}, forget plot, boxplot/draw position={position} ] table[y index=0] {{{dat_path.as_posix()}}};"
            + "\n"
        )

    if show_method_legend_table_hint:
        body += (
            r"% NOTE: method color/marker mapping is in Results/Comparison/summary_error_stats_table.tex"
            + "\n"
        )

    body += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(body, encoding="utf-8")


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
    r"""
    Draw a Normal PDF curve of an error metric for every method (mean/std taken from the data).

    Constraints:
      - No grid.
      - y-axis label: "Probability Density Function (PDF)".
      - When show_method_legends is False, emit 'forget plot' so the methods stay out of the legend.
      - The Safe/Not-safe annotations appear (without background), roughly 20% lower in y.
      - The TikZ output starts with:
            \begin{tikzpicture}
            \input{Results/Comparison/faber_plot_styles.tex}
            \begin{axis}[
    """
    if xlabel is None:
        if error_col == "Error_log10_dp":
            xlabel = r"$\log_{10}\left(\frac{DP_e}{DP}\right)$"
        elif error_col == "Error_ln_dp":
            xlabel = r"$\ln(DP/DP_{fit})$"
        else:
            xlabel = r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$"

    if error_col not in df_mat.columns:
        print(f"[WARN] {error_col} is missing, cannot draw normal curves: {material_label}")
        return

    # Damage labels stay neutral and follow the sign of the requested ratio.
    negative_label, positive_label = "Not-safe", "Safe"
    if error_col == "Error_log10_dp":
        negative_label, positive_label = r"$DP_e < DP$", r"$DP_e > DP$"
    elif error_col == "Error_ln_dp":
        negative_label, positive_label = r"$DP < DP_e$", r"$DP > DP_e$"

    # Paper style: every PDF curve is drawn with opacity=0.1
    pdf_alpha = 0.1
    pdf_lw = float(pdf_lw)

    stats: List[Tuple[str, float, float]] = []
    for method in sorted(df_mat["Method"].dropna().unique().tolist()):
        values = (
            pd.to_numeric(df_mat.loc[df_mat["Method"] == method, error_col], errors="coerce")
            .dropna()
            .to_numpy(dtype=float)
        )
        if len(values) >= 2:
            mu = float(np.nanmean(values))
            sigma = float(np.nanstd(values, ddof=1))
            if np.isfinite(mu) and np.isfinite(sigma) and sigma > 0:
                stats.append((method, mu, sigma))

    if not stats:
        print(f"[WARN] Not enough {error_col} values for normal curves: {material_label}")
        return

    x_min = min(mu - 4 * s for _, mu, s in stats)
    x_max = max(mu + 4 * s for _, mu, s in stats)
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        x_min, x_max = -3.0, 3.0
    xs = np.linspace(x_min, x_max, 600)

    # Matplotlib
    fig, ax = plt.subplots(figsize=(10, 5))
    y_max = 0.0
    for method, mu, sigma in stats:
        ys = _gaussian_density(xs, mu, sigma)
        y_max = max(y_max, float(np.nanmax(ys)))
        style = method_styles.get(method)
        color_hex = ("#" + style["hex"]) if style else None
        label = f"{method} (μ={mu:.3g}, σ={sigma:.3g})" if show_method_legends else None
        dash = (0, (3, 1)) if (style and style.get("is_ext") == "1") else "-"
        ax.plot(
            xs,
            ys,
            linewidth=pdf_lw,
            alpha=0.1,
            color=color_hex,
            linestyle=dash,
            label=label,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability Density Function (PDF)")
    ax.grid(False)
    if show_method_legends:
        ax.legend(fontsize=9)

    # Safe / Not-safe annotations (positive = conservative), shifted ~20% down in y
    if np.isfinite(y_max) and y_max > 0:
        label_y = 0.56 * y_max
        safe_x = x_max * 0.72 if x_max > 0 else x_max * 0.92
        unsafe_x = x_min * 0.72 if x_min < 0 else x_min * 0.92
        ax.text(unsafe_x, label_y, negative_label, ha="center", va="center")
        ax.text(safe_x, label_y, positive_label, ha="center", va="center")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # TikZ
    curve_files: List[Tuple[str, float, float, Path]] = []
    y_max = 0.0
    for method, mu, sigma in stats:
        ys = _gaussian_density(xs, mu, sigma)
        y_max = max(y_max, float(np.nanmax(ys)))
        dat_path = data_dir / "pdf" / f"{_fs_safe_stem(method)}.dat"
        _dump_pairs(dat_path, xs, ys)
        curve_files.append((method, mu, sigma, dat_path))

    body = (
        r"""\begin{tikzpicture}
\input{Results/Comparison/faber_plot_styles.tex}
\begin{axis}[
  width=0.5\textwidth,
  height=0.5\textwidth,
  xlabel={"""
        + xlabel
        + r"""},
  ylabel={Probability Density Function (PDF)},
"""
    )

    if show_method_legends:
        body += r"  legend style={}," + "\n"

    body += r"]" + "\n"
    body += (
        r"\node[anchor=north west, fill=white, inner sep=2pt, rounded corners=1pt] at (axis description cs:0.02,0.98) {"
        + _tex_escape(material_label)
        + r"};"
        + "\n"
    )

    forget = "" if show_method_legends else ", forget plot"

    for method, mu, sigma, dat_path in curve_files:
        style = method_styles.get(method)
        key = style["id"] if style else None
        plot_style = f"faber/method_line/{key}" if key else ""
        body += (
            r"\addplot["
            + plot_style
            + forget
            + r", opacity=0.1"
            + r", line width="
            + f"{pdf_lw:.3g}"
            + r"pt] table[x index=0, y index=1] {"
            + dat_path.as_posix()
            + r"};"
            + "\n"
        )

        if show_method_legends:
            label = f"{_tex_escape(method)} $(\\mu={mu:.3g},\\ \\sigma={sigma:.3g})$"
            body += r"\addlegendentry{" + label + r"}" + "\n"

    # Safe / Not-safe annotations (positive = conservative), shifted ~20% down in y
    if np.isfinite(y_max) and y_max > 0:
        label_y = 0.56 * y_max
        safe_x = x_max * 0.72 if x_max > 0 else x_max * 0.92
        unsafe_x = x_min * 0.72 if x_min < 0 else x_min * 0.92
        body += (
            rf"\node[anchor=south] at (axis cs:{unsafe_x:.6g},{label_y:.6g}) {{{negative_label}}};"
            + "\n"
        )
        body += (
            rf"\node[anchor=south] at (axis cs:{safe_x:.6g},{label_y:.6g}) {{{positive_label}}};"
            + "\n"
        )

    body += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(body.strip(), encoding="utf-8")


RUNOUT_FLAG = 2  # censored rows: Nf_exp is a lower bound on life (see report/export_csv.py)
RUNOUT_ARROW = 1.8  # arrow length expressed as a factor on the log axis


def _runout_points(
    configured_runouts: Optional[pd.DataFrame], method: str, x_field: str, y_field: str
):
    if (
        configured_runouts is None
        or configured_runouts.empty
        or x_field not in configured_runouts.columns
    ):
        return np.array([]), np.array([])
    subset = configured_runouts[configured_runouts["Method"] == method]
    xs = pd.to_numeric(subset[x_field], errors="coerce").to_numpy(dtype=float)
    ys = pd.to_numeric(subset[y_field], errors="coerce").to_numpy(dtype=float)
    usable = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
    return xs[usable], ys[usable]


def _runout_bounds(configured_runouts, method_names, x_field, y_field, side) -> np.ndarray:
    """Coordinates (arrow tips included) that the equal-axis range has to contain."""
    collected = []
    for method in method_names:
        xs, ys = _runout_points(configured_runouts, method, x_field, y_field)
        if xs.size:
            collected += [xs, ys, xs * RUNOUT_ARROW**side]
    return np.concatenate(collected) if collected else np.array([])


def _draw_runouts_mpl(
    ax, configured_runouts, method_names, style_map, x_field, y_field, side
) -> None:
    """Hollow markers with an arrow pointing towards the side where the true value lies."""
    first = True
    for method in method_names:
        xs, ys = _runout_points(configured_runouts, method, x_field, y_field)
        if not xs.size:
            continue
        style = style_map.get(method)
        color_hex = ("#" + style["hex"]) if style else "black"
        ax.scatter(
            xs,
            ys,
            marker="o",
            facecolors="none",
            edgecolors=color_hex,
            alpha=0.5,
            s=30,
            linewidths=0.9,
            label="Run-out" if first else "_nolegend_",
        )
        first = False
        for x_val, y_val in zip(xs, ys):
            ax.annotate(
                "",
                xy=(x_val * RUNOUT_ARROW**side, y_val),
                xytext=(x_val, y_val),
                arrowprops=dict(arrowstyle="->", color=color_hex, lw=0.7, alpha=0.5),
            )


def _runouts_as_tikz(
    configured_runouts, method_names, style_map, x_field, y_field, side, data_dir
) -> str:
    body = ""
    for method in method_names:
        xs, ys = _runout_points(configured_runouts, method, x_field, y_field)
        if not xs.size:
            continue
        dat_path = data_dir / "runouts" / f"{_fs_safe_stem(method)}.dat"
        _dump_pairs(dat_path, xs, ys)
        style = style_map.get(method)
        color = (
            f"{{rgb,255:red,{int(style['hex'][0:2], 16)};green,{int(style['hex'][2:4], 16)};blue,{int(style['hex'][4:6], 16)}}}"
            if style
            else "black"
        )
        body += (
            r"\addplot[only marks, mark=o, draw=" + color + r", opacity=0.5, forget plot] "
            r"table[x index=0, y index=1] {" + dat_path.as_posix() + r"};" + "\n"
        )
        for x_val, y_val in zip(xs, ys):
            body += (
                r"\draw[-latex, draw="
                + color
                + r", opacity=0.5] (axis cs:"
                + f"{x_val:.6e},{y_val:.6e}"
                + r") -- (axis cs:"
                + f"{x_val * RUNOUT_ARROW**side:.6e},{y_val:.6e}"
                + r");"
                + "\n"
            )
    return body


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
    df_run: Optional[pd.DataFrame] = None,
) -> None:
    if "Nf_expected" not in df_mat.columns or "Nf_exp" not in df_mat.columns:
        print(f"[WARN] Nf_expected/Nf_exp missing, no scatter plot: {material_label}")
        return

    # Paper style: every scatter marker uses opacity=0.2 and fill opacity=0.2
    scatter_alpha = 0.2
    factors = sorted(set(float(f) for f in factors if np.isfinite(f) and f > 1.0))

    method_names = sorted(df_mat["Method"].dropna().unique().tolist())

    x_all = pd.to_numeric(df_mat["Nf_exp"], errors="coerce").to_numpy(dtype=float)
    y_all = pd.to_numeric(df_mat["Nf_expected"], errors="coerce").to_numpy(dtype=float)
    usable = np.isfinite(x_all) & np.isfinite(y_all) & (x_all > 0) & (y_all > 0)
    if not np.any(usable):
        print(f"[WARN] No positive finite values for the scatter: {material_label}")
        return

    x_all = x_all[usable]
    y_all = y_all[usable]
    extra = _runout_bounds(df_run, method_names, "Nf_exp", "Nf_expected", +1)
    if extra.size:
        x_all = np.concatenate([x_all, extra])

    # Equal range with padding applied in log-space.
    # The fixed Safe/Not-safe anchor positions required by the paper style are included too.
    safe_x, safe_y = 3_000_000.0, 200.0
    unsafe_x, unsafe_y = 200.0, 3_000_000.0

    raw_low = float(min(x_all.min(), y_all.min(), safe_x, safe_y, unsafe_x, unsafe_y))
    raw_high = float(max(x_all.max(), y_all.max(), safe_x, safe_y, unsafe_x, unsafe_y))

    log_min = float(np.log10(raw_low))
    log_max = float(np.log10(raw_high))
    span = (log_max - log_min) if log_max > log_min else 1.0
    padding = 0.15 * span + 0.05
    lo = 10 ** (log_min - padding)
    hi = 10 ** (log_max + padding)

    # Greys for the factor bands: kept consistent, differentiated by line style as well
    grey_levels = [0.45, 0.65, 0.80, 0.30, 0.55]
    factor_grey: Dict[float, float] = {
        f: grey_levels[idx % len(grey_levels)] for idx, f in enumerate(factors)
    }

    # --- Matplotlib ---
    fig, ax = plt.subplots(figsize=(7, 7))

    for method in method_names:
        subset = df_mat[df_mat["Method"] == method]
        xs = pd.to_numeric(subset["Nf_exp"], errors="coerce").to_numpy(dtype=float)
        ys = pd.to_numeric(subset["Nf_expected"], errors="coerce").to_numpy(dtype=float)
        keep = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
        if np.any(keep):
            style = method_styles.get(method)
            color_hex = ("#" + style["hex"]) if style else None
            label = method if show_method_legends else "_nolegend_"
            marker = pgf_mark_to_mpl(style["mark"]) if style else "o"
            extended = bool(style and style.get("is_ext") == "1")
            face = "none" if extended or not color_hex else color_hex
            ax.scatter(
                xs[keep],
                ys[keep],
                label=label,
                alpha=0.2,
                s=30,
                marker=marker,
                facecolors=face,
                edgecolors=color_hex,
                linewidths=0.9,
            )

    _draw_runouts_mpl(ax, df_run, method_names, method_styles, "Nf_exp", "Nf_expected", +1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    diagonal = np.array([lo, hi], dtype=float)

    # Bisector (solid black) with a legend entry
    ax.plot(diagonal, diagonal, color="black", linewidth=1.8, label="Bisector")

    # Factor bands: distinguishable styles (dense dashed vs dense dotted ...)
    for idx, f in enumerate(factors):
        grey_text = str(factor_grey[f])  # grayscale string
        factor_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        dash = _factor_mpl_dash(idx)
        ax.plot(
            diagonal,
            f * diagonal,
            linestyle=dash,
            color=grey_text,
            linewidth=1.3,
            label=factor_label,
        )
        ax.plot(
            diagonal,
            diagonal / f,
            linestyle=dash,
            color=grey_text,
            linewidth=1.3,
            label="_nolegend_",
        )

    # Safe / Not-safe
    ax.text(
        safe_x,
        safe_y,
        "Safe",
        rotation=45,
        color="black",
        fontsize=11,
        ha="right",
        va="bottom",
    )
    ax.text(
        unsafe_x,
        unsafe_y,
        "Not-safe",
        rotation=45,
        color="black",
        fontsize=11,
        ha="left",
        va="top",
    )

    ax.set_xlabel(r"Experimental number of cycles to failure - $N_f$")
    ax.set_ylabel(r"Expected number of cycles to failure - $N_{f,e}$")
    ax.legend(fontsize=9, loc="lower left")  # only Bisector + Factors (and methods when enabled)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # --- TikZ export ---
    scatter_files = []
    for method in method_names:
        subset = df_mat[df_mat["Method"] == method]
        xs = pd.to_numeric(subset["Nf_exp"], errors="coerce").to_numpy(dtype=float)
        ys = pd.to_numeric(subset["Nf_expected"], errors="coerce").to_numpy(dtype=float)
        keep = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
        if np.any(keep):
            dat_path = data_dir / "scatter" / f"{_fs_safe_stem(method)}.dat"
            _dump_pairs(dat_path, xs[keep], ys[keep])
            scatter_files.append((method, dat_path))

    def _grey_to_pgf(level01: float) -> str:
        percent = int(round(level01 * 100))
        return f"gray!{percent}"

    body = (
        r"""\begin{tikzpicture}
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
        + _tex_escape(material_label)
        + r"};"
        + "\n"
    )

    forget = "" if show_method_legends else ", forget plot"

    for method, dat_path in scatter_files:
        style = method_styles.get(method)
        key = style["id"] if style else None
        plot_style = f"faber/method_scatter/{key}" if key else "only marks"
        body += (
            r"\addplot["
            + plot_style
            + forget
            + r", opacity=0.2, fill opacity=0.2"
            + r"] table[x index=0, y index=1] {"
            + dat_path.as_posix()
            + r"};"
            + "\n"
        )
        if show_method_legends:
            body += r"\addlegendentry{" + _tex_escape(method) + r"}" + "\n"

    # Bisector with its legend entry
    body += (
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

    # Factor bands with distinguishable styles, one legend entry per factor
    for idx, f in enumerate(factors):
        grey_level = factor_grey[f]
        pgf_grey = _grey_to_pgf(grey_level)
        factor_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        dash = _factor_pgf_dash(idx)

        # y = f x (legend entry)
        body += (
            r"\addplot["
            + pgf_grey
            + ", "
            + dash
            + r", line width=0.9pt] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(f * lo):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(f * hi):.6e}"
            + ")};"
            "\n"
            r"\addlegendentry{" + factor_label + r"}"
            "\n"
        )
        # y = x/f (no legend entry)
        body += (
            r"\addplot["
            + pgf_grey
            + ", "
            + dash
            + r", line width=0.9pt, forget plot] coordinates {("
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

    body += _runouts_as_tikz(
        df_run, method_names, method_styles, "Nf_exp", "Nf_expected", +1, data_dir
    )

    # Safe / Not-safe labels
    body += (
        r"\node[rotate=45, anchor=south east] at (axis cs:"
        + f"{safe_x:.6e}"
        + ","
        + f"{safe_y:.6e}"
        + r") {Safe};"
        "\n"
        r"\node[rotate=45, anchor=north west] at (axis cs:"
        + f"{unsafe_x:.6e}"
        + ","
        + f"{unsafe_y:.6e}"
        + r") {Not-safe};"
        "\n"
    )

    body += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(body.strip(), encoding="utf-8")


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
    df_run: Optional[pd.DataFrame] = None,
) -> None:
    """Damage-parameter parity plot valid for every slope/metric mode.

    Draws:
      y = DP_point (damage parameter produced by the model evaluation)
      x = DP_fit   (calibration best-fit evaluated at Nf_exp)

    With these axes the bisector and the factor bands are meaningful again.
    """
    if "DP_fit" not in df_mat.columns or "CP_value" not in df_mat.columns:
        print(f"[WARN] DP_fit/CP_value missing, no DP parity plot: {material_label}")
        return

    # Follow the paper style (same as the Nf parity scatter)
    scatter_alpha = 0.2
    factors = sorted(set(float(f) for f in factors if np.isfinite(f) and f > 1.0))

    method_names = sorted(df_mat["Method"].dropna().unique().tolist())

    x_all = pd.to_numeric(df_mat["DP_fit"], errors="coerce").to_numpy(dtype=float)
    y_all = pd.to_numeric(df_mat["CP_value"], errors="coerce").to_numpy(dtype=float)
    usable = np.isfinite(x_all) & np.isfinite(y_all) & (x_all > 0) & (y_all > 0)
    if not np.any(usable):
        print(f"[WARN] No positive finite DP values for the parity plot: {material_label}")
        return

    x_all = x_all[usable]
    y_all = y_all[usable]
    extra = _runout_bounds(df_run, method_names, "DP_fit", "CP_value", -1)
    if extra.size:
        x_all = np.concatenate([x_all, extra])

    raw_low = float(min(x_all.min(), y_all.min()))
    raw_high = float(max(x_all.max(), y_all.max()))

    log_min = float(np.log10(raw_low))
    log_max = float(np.log10(raw_high))
    span = (log_max - log_min) if log_max > log_min else 1.0
    padding = 0.15 * span + 0.05
    lo = 10 ** (log_min - padding)
    hi = 10 ** (log_max + padding)

    # Neutral DP annotations also work when the master curve is flat.
    safe_x = 10 ** (log_min + 0.20 * span)
    safe_y = 10 ** (log_min + 0.80 * span)
    unsafe_x = 10 ** (log_min + 0.80 * span)
    unsafe_y = 10 ** (log_min + 0.20 * span)

    lo = float(min(lo, safe_x, safe_y, unsafe_x, unsafe_y))
    hi = float(max(hi, safe_x, safe_y, unsafe_x, unsafe_y))

    # Greys for the factor bands
    grey_levels = [0.45, 0.65, 0.80, 0.30, 0.55]
    factor_grey: Dict[float, float] = {
        f: grey_levels[idx % len(grey_levels)] for idx, f in enumerate(factors)
    }

    # --- Matplotlib ---
    fig, ax = plt.subplots(figsize=(7, 7))

    for method in method_names:
        subset = df_mat[df_mat["Method"] == method]
        xs = pd.to_numeric(subset["DP_fit"], errors="coerce").to_numpy(dtype=float)
        ys = pd.to_numeric(subset["CP_value"], errors="coerce").to_numpy(dtype=float)
        keep = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
        if np.any(keep):
            style = method_styles.get(method)
            color_hex = ("#" + style["hex"]) if style else None
            label = method if show_method_legends else "_nolegend_"
            marker = pgf_mark_to_mpl(style["mark"]) if style else "o"
            extended = bool(style and style.get("is_ext") == "1")
            face = "none" if extended or not color_hex else color_hex
            ax.scatter(
                xs[keep],
                ys[keep],
                label=label,
                alpha=0.2,
                s=30,
                marker=marker,
                facecolors=face,
                edgecolors=color_hex,
                linewidths=0.9,
            )

    _draw_runouts_mpl(ax, df_run, method_names, method_styles, "DP_fit", "CP_value", -1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    diagonal = np.array([lo, hi], dtype=float)
    ax.plot(diagonal, diagonal, color="black", linewidth=1.8, label="Bisector")

    for idx, f in enumerate(factors):
        grey_text = str(factor_grey[f])
        factor_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        dash = _factor_mpl_dash(idx)
        ax.plot(
            diagonal,
            f * diagonal,
            linestyle=dash,
            color=grey_text,
            linewidth=1.3,
            label=factor_label,
        )
        ax.plot(
            diagonal,
            diagonal / f,
            linestyle=dash,
            color=grey_text,
            linewidth=1.3,
            label="_nolegend_",
        )

    ax.text(
        safe_x,
        safe_y,
        r"$DP > DP_e$",
        rotation=45,
        color="black",
        fontsize=11,
        ha="left",
        va="top",
    )
    ax.text(
        unsafe_x,
        unsafe_y,
        r"$DP < DP_e$",
        rotation=45,
        color="black",
        fontsize=11,
        ha="right",
        va="bottom",
    )

    ax.set_xlabel(r"Calibration best-fit damage parameter - $DP_{fit}$")
    ax.set_ylabel(r"Damage parameter from loading conditions - $DP$")
    ax.legend(fontsize=9, loc="lower left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # --- TikZ export ---
    scatter_files = []
    for method in method_names:
        subset = df_mat[df_mat["Method"] == method]
        xs = pd.to_numeric(subset["DP_fit"], errors="coerce").to_numpy(dtype=float)
        ys = pd.to_numeric(subset["CP_value"], errors="coerce").to_numpy(dtype=float)
        keep = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
        if np.any(keep):
            dat_path = data_dir / "scatter" / f"{_fs_safe_stem(method)}.dat"
            _dump_pairs(dat_path, xs[keep], ys[keep])
            scatter_files.append((method, dat_path))

    def _grey_to_pgf(level01: float) -> str:
        percent = int(round(level01 * 100))
        return f"gray!{percent}"

    body = (
        r"""\begin{tikzpicture}
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
        + _tex_escape(material_label)
        + r"};"
        + "\n"
    )

    forget = "" if show_method_legends else ", forget plot"

    for method, dat_path in scatter_files:
        style = method_styles.get(method)
        key = style["id"] if style else None
        plot_style = f"faber/method_scatter/{key}" if key else "only marks"
        body += (
            r"\addplot["
            + plot_style
            + forget
            + r", opacity=0.2, fill opacity=0.2"
            + r"] table[x index=0, y index=1] {"
            + dat_path.as_posix()
            + r"};"
            + "\n"
        )
        if show_method_legends:
            body += r"\addlegendentry{" + _tex_escape(method) + r"}" + "\n"

    # Bisector
    body += (
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

    # Factor bands
    for idx, f in enumerate(factors):
        grey_level = factor_grey[f]
        pgf_grey = _grey_to_pgf(grey_level)
        factor_label = f"Factor {int(f) if abs(f-round(f))<1e-9 else f:g}"
        dash = _factor_pgf_dash(idx)

        # NOTE: pgfplots wants a single brace group for the coordinates.
        # Never emit doubled braces (e.g. `coordinates {{...}}`), which triggers
        # `Argument of \\pgfplots@addplotimpl@coordinates has an extra }`.
        body += (
            r"\addplot["
            + pgf_grey
            + ", "
            + dash
            + r", line width=0.9pt] coordinates {("
            + f"{lo:.6e}"
            + ","
            + f"{(f * lo):.6e}"
            + ")("
            + f"{hi:.6e}"
            + ","
            + f"{(f * hi):.6e}"
            + ")};"
            "\n"
            r"\addlegendentry{" + factor_label + r"}"
            "\n"
        )
        body += (
            r"\addplot["
            + pgf_grey
            + ", "
            + dash
            + r", line width=0.9pt, forget plot] coordinates {("
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

    body += _runouts_as_tikz(
        df_run, method_names, method_styles, "DP_fit", "CP_value", -1, data_dir
    )

    # Safe / Not-safe labels
    body += (
        r"\node[rotate=45, anchor=north west] at (axis cs:"
        + f"{safe_x:.6e}"
        + ","
        + f"{safe_y:.6e}"
        + r") {$DP > DP_e$};"
        "\n"
        r"\node[rotate=45, anchor=south east] at (axis cs:"
        + f"{unsafe_x:.6e}"
        + ","
        + f"{unsafe_y:.6e}"
        + r") {$DP < DP_e$};"
        "\n"
    )

    body += r"""\end{axis}
\end{tikzpicture}
"""
    out_tex.parent.mkdir(parents=True, exist_ok=True)
    out_tex.write_text(body.strip(), encoding="utf-8")


def compute_stats_table(df: pd.DataFrame, plot_kind: str = "both") -> pd.DataFrame:
    """DP log10 statistics for every slope, plus life errors where they make sense."""
    df = df.copy()
    if "Metric_mode" not in df:
        df["Metric_mode"] = "Nf"
    if "DP_fit" in df and "CP_value" in df:
        df["Error_log10_dp"] = log10_dp_ratio(df["DP_fit"], df["CP_value"])
    records = []
    for (material, method_name, mode_name), group in df.groupby(
        ["Material", "Method", "Metric_mode"], dropna=False
    ):
        dp_ok = (
            bool(group["DP_comparison_supported"].any())
            if "DP_comparison_supported" in group
            else True
        )
        life_ok = (
            bool(group["Life_comparison_supported"].any())
            if "Life_comparison_supported" in group
            else str(mode_name).strip().upper() != "DP_DIFF"
        )
        metric_names = []
        if dp_ok and plot_kind in ("both", "damage"):
            metric_names.append("Error_log10_dp")
        if life_ok and plot_kind in ("both", "life"):
            metric_names.append("Error_ln")
        for metric_name in metric_names:
            errors = pd.to_numeric(
                group.get(metric_name, pd.Series(np.nan, index=group.index)),
                errors="coerce",
            ).to_numpy(dtype=float)
            errors = errors[np.isfinite(errors)]
            records.append(
                dict(
                    Material=material,
                    Method=method_name,
                    Metric=metric_name,
                    n=int(errors.size),
                    Error_mean=float(np.mean(errors)) if errors.size else np.nan,
                    Error_std=(float(np.std(errors, ddof=1)) if errors.size > 1 else np.nan),
                )
            )
    return pd.DataFrame(
        records,
        columns=["Material", "Method", "Metric", "n", "Error_mean", "Error_std"],
    ).sort_values(["Material", "Method", "Metric"])


def write_stats_table_with_styles(
    stats: pd.DataFrame,
    method_styles: Dict[str, Dict[str, str]],
    out_csv: Path,
    out_tex: Path,
) -> None:
    """Write one global LaTeX table (independent from the per-material folders).

    Contents:
      - the colour and marker associated with each method (matching the figures)
      - mean/std of the active error metric per (material, method)
        * Metric='Error_ln' -> Error_ln = ln(Nf / Nf,e)
        * Metric='Error_log10_dp' -> log10(DP_e / DP), all slopes
    """
    table = stats.copy()
    table["Style_id"] = table["Method"].map(lambda name: method_styles.get(name, {}).get("id", ""))
    table["Color_hex"] = table["Method"].map(
        lambda name: method_styles.get(name, {}).get("hex", "")
    )
    table["PGF_color"] = table["Method"].map(
        lambda name: method_styles.get(name, {}).get("color", "")
    )
    table["PGF_mark"] = table["Method"].map(
        lambda name: method_styles.get(name, {}).get("mark", "")
    )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)

    # LaTeX table
    out_tex.parent.mkdir(parents=True, exist_ok=True)

    output_lines: List[str] = []
    output_lines += [
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

    def _format_value(value: float) -> str:
        if not np.isfinite(value):
            return "--"
        # compact yet stable (scientific notation avoided for typical values)
        return f"{value:.4g}"

    for _, row in table.sort_values(["Method", "Material", "Metric"]).iterrows():
        material_cell = _tex_escape(row["Material"])
        method_cell = _tex_escape(row["Method"])
        method_raw = str(row["Method"])
        extended = method_raw.endswith("_ext")
        metric_raw = str(row.get("Metric", "Error_ln"))
        if metric_raw == "Error_log10_dp":
            metric_cell = r"$\log_{10}(DP_e/DP)$"
        elif metric_raw == "Error_ln_dp":
            metric_cell = r"$Error_{\ln,dp}$"
        else:
            metric_cell = r"$Error_{\ln}$"
        color_name = row["PGF_color"]
        marker_name = row["PGF_mark"]
        mu = _format_value(float(row["Error_mean"])) if pd.notna(row["Error_mean"]) else "--"
        sigma = _format_value(float(row["Error_std"])) if pd.notna(row["Error_std"]) else "--"

        # Table swatches:
        #   - Colour column: shows ONLY the line style (solid vs densely dashed) in the method colour.
        #   - Marker column: shows ONLY the marker (filled for base, white-filled for _ext), no line.
        line_style = "densely dashed" if extended else "solid"
        color_swatch = (
            r"\tikz[baseline=-0.6ex]{\draw["
            + color_name
            + r","
            + line_style
            + r", line width=0.9pt] (0,0) -- (0.5,0);}"  # noqa: E501
            if isinstance(color_name, str) and color_name
            else "--"
        )
        marker_swatch = (
            r"\tikz[baseline=-0.6ex]{"
            + r"\draw["
            + color_name
            + r", mark="
            + marker_name
            + r", mark options={solid, fill="
            + ("white" if extended else color_name)
            + r"}] plot coordinates {(0,0)};"
            + r"}"
            if isinstance(color_name, str)
            and color_name
            and isinstance(marker_name, str)
            and marker_name
            else "--"
        )

        mu_cell = (rf"$ {mu} $") if mu != "--" else "--"
        sigma_cell = (rf"$ {sigma} $") if sigma != "--" else "--"
        output_lines.append(
            f"{material_cell} & {method_cell} & {metric_cell} & {color_swatch} & {marker_swatch} & {mu_cell} & {sigma_cell}\\\\"
        )

    output_lines += [r"\bottomrule", r"\end{tabularx}", ""]
    out_tex.write_text("\n".join(output_lines), encoding="utf-8")


# =========================
# Main
# =========================
def main(argv=None) -> None:
    parser = argparse.ArgumentParser(allow_abbrev=False)  # needed for Jupyter (-f kernel.json)
    parser.add_argument("--results_dir", type=str, default="Results")
    parser.add_argument("--out_dir", type=str, default="Results/Comparison")
    parser.add_argument(
        "--include-legacy",
        action="store_true",
        help="Also keep pre-audit results; methods that were disabled remain excluded.",
    )
    parser.add_argument(
        "--flag_filter", type=str, default="all", choices=["all", "design", "other"]
    )
    parser.add_argument(
        "--hide-runouts",
        action="store_true",
        help="Skip run-outs (Flag 2) when drawing figures.",
    )

    parser.add_argument(
        "--common-methods",
        nargs="+",
        default=None,
        help="Restrict the comparison to the failure observations assessed by every listed method "
        "(upstream 904c619); writes selected_common_observations.csv and the subset reports.",
    )
    parser.add_argument("--comparison-metric", choices=["life", "damage"], default="life")
    parser.add_argument(
        "--plot-kind",
        choices=["life", "damage", "both"],
        default="both",
        help="Generate life-based plots, DP-based plots, or both (default).",
    )

    # opacity of the scatter markers (helps seeing clusters)
    parser.add_argument("--scatter_alpha", type=float, default=0.18)

    # appearance of the normal curves
    parser.add_argument("--pdf_lw", type=float, default=2.0)
    parser.add_argument("--pdf_alpha", type=float, default=0.85)

    # factor bands (defaults: 3 and 5)
    parser.add_argument("--factors", type=float, nargs="*", default=[3.0, 5.0])

    # Method legends are optional because a global table is produced anyway
    parser.add_argument(
        "--show_method_legends",
        action="store_true",
        help="Add method entries to the plot legends.",
    )

    cli_args, _unrecognized = parser.parse_known_args(argv)

    results_root = Path(cli_args.results_dir)
    output_root = Path(cli_args.out_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    scatter_alpha = float(np.clip(cli_args.scatter_alpha, 0.0, 1.0))
    pdf_alpha = float(np.clip(cli_args.pdf_alpha, 0.0, 1.0))
    pdf_lw = float(max(cli_args.pdf_lw, 0.1))
    factors = [float(f) for f in (cli_args.factors or []) if np.isfinite(f) and f > 1.0]

    frame = collect_all_results(results_root)
    disabled = frame["Method"].astype(str).str.removesuffix("_ext").isin(DISABLED_METHODS)
    if disabled.any():
        print(
            f"[WARN] Dropping {int(disabled.sum())} rows from disabled methods; see docs/method-audit.md"
        )
    frame = frame.loc[~disabled].copy()
    revision = frame.get("Implementation_revision", pd.Series("legacy", index=frame.index))
    legacy = revision != IMPLEMENTATION_REVISION
    if legacy.any():
        print(
            f"[WARN] {int(legacy.sum())} rows predate this audit; rerun the analysis to pick up the corrected equations."
        )
        if not cli_args.include_legacy:
            frame = frame.loc[~legacy].copy()
    if frame.empty:
        raise RuntimeError(
            "No current results available. Run multifat-run (or scripts/1_run_material.py), or pass --include-legacy to also use historical outputs."
        )

    frame = annotate_comparison_capabilities(frame)
    write_native_life_summary(frame, output_root)
    capability_report = frame[
        [
            "Material",
            "Method",
            "Metric_mode",
            "Life_prediction_source",
            "DP_reference_source",
            "Life_comparison_supported",
            "DP_comparison_supported",
        ]
    ].drop_duplicates()
    capability_report.to_csv(output_root / "comparison_capabilities.csv", index=False)
    if cli_args.common_methods:
        supported = (
            "Life_comparison_supported"
            if cli_args.comparison_metric == "life"
            else "DP_comparison_supported"
        )
        unavailable = [
            name
            for name in cli_args.common_methods
            if not frame.loc[frame["Method"] == name, supported].any()
        ]
        if unavailable:
            raise ValueError(
                f"Methods without a valid {cli_args.comparison_metric} comparison: "
                f"{', '.join(unavailable)}"
            )
    if cli_args.plot_kind == "life":
        supported = frame["Life_comparison_supported"]
    elif cli_args.plot_kind == "damage":
        supported = frame["DP_comparison_supported"]
    else:
        supported = frame["Life_comparison_supported"] | frame["DP_comparison_supported"]
    excluded = int((~supported).sum())
    if excluded:
        print(f"[INFO] Excluded {excluded} rows without a valid {cli_args.plot_kind} comparison")
    frame = frame.loc[supported].copy()
    if frame.empty:
        print("[INFO] No compatible methods for the selected plot kind.")
        return

    # Minimal set of required columns (works for both the Nf-based and DP-based modes)
    for needed in ["Material", "Method", "Nf_exp"]:
        if needed not in frame.columns:
            raise RuntimeError(f"Required column '{needed}' is missing after reading the CSVs.")
    if "CP_value" not in frame.columns:
        raise RuntimeError(
            "Required column 'CP_value' (damage parameter) is missing after reading the CSVs."
        )

    if "Metric_mode" not in frame.columns:
        frame["Metric_mode"] = "Nf"
    else:
        frame["Metric_mode"] = frame["Metric_mode"].fillna("Nf").astype(str)

    if "Error_ln" not in frame.columns:
        frame["Error_ln"] = np.nan
        if "Nf_expected" in frame.columns:
            usable = (
                (frame["Nf_exp"] > 0)
                & (frame["Nf_expected"] > 0)
                & np.isfinite(frame["Nf_exp"])
                & np.isfinite(frame["Nf_expected"])
            )
            frame.loc[usable, "Error_ln"] = np.log(
                frame.loc[usable, "Nf_exp"] / frame.loc[usable, "Nf_expected"]
            )

    # DP log-ratio error (used when Metric_mode == 'DP_DIFF')
    if "Error_ln_dp" not in frame.columns:
        frame["Error_ln_dp"] = np.nan
    if "DP_fit" in frame.columns:
        usable = (
            np.isfinite(frame["DP_fit"])
            & np.isfinite(frame["CP_value"])
            & (frame["DP_fit"] > 0.0)
            & (frame["CP_value"] > 0.0)
        )
        frame.loc[usable, "Error_ln_dp"] = np.log(
            frame.loc[usable, "CP_value"] / frame.loc[usable, "DP_fit"]
        )
        frame["Error_log10_dp"] = log10_dp_ratio(frame["DP_fit"], frame["CP_value"])

    if "Error_dp" not in frame.columns:
        frame["Error_dp"] = np.nan
    if "DP_fit" in frame.columns:
        usable = np.isfinite(frame["DP_fit"]) & np.isfinite(frame["CP_value"])
        frame.loc[usable, "Error_dp"] = frame.loc[usable, "CP_value"] - frame.loc[usable, "DP_fit"]

    # Run-outs (Flag 2) are lower bounds on life: shown with arrows, kept out of the statistics.
    if "Flag" in frame.columns:
        runout_mask = pd.to_numeric(frame["Flag"], errors="coerce") == RUNOUT_FLAG
        runouts_frame = frame.loc[runout_mask].copy()
        frame = frame.loc[~runout_mask].copy()
    else:
        runouts_frame = frame.iloc[0:0].copy()
    if cli_args.hide_runouts:
        runouts_frame = runouts_frame.iloc[0:0]

    if cli_args.flag_filter != "all" and "Flag" in frame.columns:
        if cli_args.flag_filter == "design":
            frame = frame[frame["Flag"] == 0]
        elif cli_args.flag_filter == "other":
            frame = frame[frame["Flag"] == 1]

    if cli_args.common_methods:
        from multifat.report.subsets import common_subset, write_subset_reports

        write_subset_reports(frame, output_root, cli_args.common_methods)
        frame = common_subset(frame, cli_args.common_methods, cli_args.comparison_metric)
        if frame.empty:
            raise RuntimeError(f"No common observations for {cli_args.common_methods}")
        frame.to_csv(output_root / "selected_common_observations.csv", index=False)
        obs_per_method = frame.groupby("Method").size().iloc[0]
        (output_root / "subset_definition.txt").write_text(
            f"methods: {' '.join(cli_args.common_methods)}\n"
            f"metric: {cli_args.comparison_metric}\n"
            f"flag filter: {cli_args.flag_filter}\n"
            f"observations per method: {obs_per_method} (failures with a finite metric in every method)\n"
            f"materials: {', '.join(sorted(frame['Material'].unique()))}\n",
            encoding="utf-8",
        )
        runouts_frame = runouts_frame[runouts_frame["Method"].isin(cli_args.common_methods)]

    # Build the global palette:
    # - colours per Method_family
    # - markers per Method
    method_styles, family_styles = build_style_maps(frame)

    # Emit the global style file (used by every figure)
    write_faber_style_file(method_styles, family_styles, output_root / "faber_plot_styles.tex")

    # Global statistics (feeding the single "legend table")
    stats = compute_stats_table(frame, cli_args.plot_kind)
    stats.to_csv(output_root / "summary_error_stats.csv", index=False)
    write_stats_table_with_styles(
        stats=stats,
        method_styles=method_styles,
        out_csv=output_root / "summary_error_stats_with_styles.csv",
        out_tex=output_root / "summary_error_stats_table.tex",
    )

    materials = sorted(frame["Material"].dropna().unique().tolist())
    if not materials:
        raise RuntimeError("No materials found in the loaded CSVs.")

    for material in materials:
        material_frame = frame[frame["Material"] == material].copy()
        if material_frame.empty:
            continue
        material_runouts = runouts_frame[runouts_frame["Material"] == material]

        # Filter each output by the declared method capabilities.
        life_frame = material_frame[material_frame["Life_comparison_supported"]].copy()
        damage_frame = material_frame[material_frame["DP_comparison_supported"]].copy()
        runouts_nf = material_runouts[material_runouts["Life_comparison_supported"]].copy()
        runouts_dp = material_runouts[material_runouts["DP_comparison_supported"]].copy()
        if cli_args.plot_kind == "life":
            damage_frame = damage_frame.iloc[0:0]
        elif cli_args.plot_kind == "damage":
            life_frame = life_frame.iloc[0:0]

        stem = _fs_safe_stem(material)
        material_dir = output_root / stem
        material_dir.mkdir(parents=True, exist_ok=True)

        tikz_dir = material_dir / "tikz_data"
        tikz_dir.mkdir(parents=True, exist_ok=True)
        tikz_nf_dir = tikz_dir / "nf"
        tikz_dp_dir = tikz_dir / "dp"
        tikz_nf_dir.mkdir(parents=True, exist_ok=True)
        tikz_dp_dir.mkdir(parents=True, exist_ok=True)

        # -------------------------
        # Life-based figures (standard)
        # -------------------------
        if not life_frame.empty:
            plot_boxplot_errors(
                df_mat=life_frame,
                out_png=material_dir / f"{stem}_boxplot_error.png",
                out_tex=material_dir / f"{stem}_boxplot_error.tex",
                data_dir=tikz_nf_dir,
                material_label=material,
                method_styles=method_styles,
                error_col="Error_ln",
                ylabel=r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$",
            )

            plot_normal_curves(
                df_mat=life_frame,
                out_png=material_dir / f"{stem}_normal_error_pdf.png",
                out_tex=material_dir / f"{stem}_normal_error_pdf.tex",
                data_dir=tikz_nf_dir,
                material_label=material,
                method_styles=method_styles,
                pdf_lw=pdf_lw,
                pdf_alpha=pdf_alpha,
                show_method_legends=cli_args.show_method_legends,
                error_col="Error_ln",
                xlabel=r"$Error_{\ln}=\ln\left(\frac{N_f}{N_{f,e}}\right)$",
            )

            plot_expected_vs_experimental(
                df_mat=life_frame,
                out_png=material_dir / f"{stem}_expected_vs_experimental.png",
                out_tex=material_dir / f"{stem}_expected_vs_experimental.tex",
                data_dir=tikz_nf_dir,
                material_label=material,
                method_styles=method_styles,
                scatter_alpha=scatter_alpha,
                factors=factors,
                show_method_legends=cli_args.show_method_legends,
                df_run=runouts_nf,
            )

        # -------------------------
        # Damage-based figures (all slopes)
        # -------------------------
        if not damage_frame.empty:
            plot_boxplot_errors(
                df_mat=damage_frame,
                out_png=material_dir / f"{stem}_dp_boxplot_error.png",
                out_tex=material_dir / f"{stem}_dp_boxplot_error.tex",
                data_dir=tikz_dp_dir,
                material_label=material,
                method_styles=method_styles,
                error_col="Error_log10_dp",
                ylabel=r"$\log_{10}\left(\frac{DP_e}{DP}\right)$",
            )

            plot_normal_curves(
                df_mat=damage_frame,
                out_png=material_dir / f"{stem}_dp_error_pdf.png",
                out_tex=material_dir / f"{stem}_dp_error_pdf.tex",
                data_dir=tikz_dp_dir,
                material_label=material,
                method_styles=method_styles,
                pdf_lw=pdf_lw,
                pdf_alpha=pdf_alpha,
                show_method_legends=cli_args.show_method_legends,
                error_col="Error_log10_dp",
                xlabel=r"$\log_{10}\left(\frac{DP_e}{DP}\right)$",
            )

            plot_dp_point_vs_fit(
                df_mat=damage_frame,
                out_png=material_dir / f"{stem}_dp_point_vs_fit.png",
                out_tex=material_dir / f"{stem}_dp_point_vs_fit.tex",
                data_dir=tikz_dp_dir,
                material_label=material,
                method_styles=method_styles,
                scatter_alpha=scatter_alpha,
                factors=factors,
                show_method_legends=cli_args.show_method_legends,
                df_run=runouts_dp,
            )

    print(f"Done. Outputs written to: {output_root.resolve()}")


if __name__ == "__main__":
    main()
