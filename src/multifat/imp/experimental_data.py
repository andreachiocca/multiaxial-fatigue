from __future__ import annotations

from pathlib import Path
import io
from typing import Any, Dict, List, Optional, Tuple, Set
import numpy as np
import math
import re
import pandas as pd
from openpyxl import load_workbook

from multifat.imp.overlays import apply_overlay_to_workbook

#: Default for the ``overlay`` switch of the importers; the command-line tools set it from
#: ``--no-overlay`` so that every workbook load in one run sees the same data.
OVERLAYS_ENABLED = True

from multifat.imp.schema import (
    DATASET_OPTIONAL,
    DATASET_REQUIRED,
    SUMMARY_CONSTANTS,
    SUMMARY_FIELDS,
    ValidationReport,
    WorkbookError,
    find_col,
    is_no_value,
    is_nullish,
    is_yes_value,
    match_constant_columns,
    material_id,
    workbook_files,
    norm_text,
    to_float_or_none,
    to_float_or_zero,
)


def build_stress_tensor(sig_rad: float, sig_ax: float, tau: float) -> np.ndarray:
    return np.array([[sig_rad, tau, 0.0], [tau, sig_ax, 0.0], [0.0, 0.0, 0.0]], dtype=float)


def stress_to_strain_iso(stress: np.ndarray, E: float | None, nu: float | None) -> np.ndarray:
    if E is None or nu is None:  # strains undefined; stress-only models can still run
        return np.full((3, 3), np.nan)
    sig_xx, sig_yy, sig_zz = stress[0, 0], stress[1, 1], stress[2, 2]
    tau_xy, tau_yz, tau_zx = stress[0, 1], stress[1, 2], stress[2, 0]

    G = E / (2.0 * (1.0 + nu))

    eps_xx = (sig_xx - nu * (sig_yy + sig_zz)) / E
    eps_yy = (sig_yy - nu * (sig_xx + sig_zz)) / E
    eps_zz = (sig_zz - nu * (sig_xx + sig_yy)) / E

    eps_xy = tau_xy / (2.0 * G)
    eps_yz = tau_yz / (2.0 * G)
    eps_zx = tau_zx / (2.0 * G)

    return np.array(
        [[eps_xx, eps_xy, eps_zx], [eps_xy, eps_yy, eps_yz], [eps_zx, eps_yz, eps_zz]],
        dtype=float,
    )


# ---- Summary helpers for the "Average" named cell rule ----
def _get_named_cell_location(wb, defined_name: str) -> Optional[Tuple[str, str]]:
    if defined_name not in wb.defined_names:
        return None
    dn = wb.defined_names[defined_name]
    dests = list(dn.destinations)
    if not dests:
        return None
    return dests[0][0], dests[0][1]


def _find_row_with_text(ws, text: str, search_rows: int = 300) -> Optional[int]:
    tnorm = norm_text(text)
    for r in range(1, search_rows + 1):
        for cell in ws[r]:
            if norm_text(cell.value) == tnorm:
                return r
    return None


def _find_header_row_and_col_indices(
    ws, header_contains: str, search_rows: int = 80
) -> Optional[Tuple[int, int]]:
    """Header row and column of the first cell containing `header_contains`, else None."""
    target = header_contains.lower()
    for r in range(1, search_rows + 1):
        for c_idx, cell in enumerate(ws[r], start=1):
            if target in norm_text(cell.value):
                return r, c_idx
    return None


def _declared_param_names() -> set[str]:
    """Constants declared by registered models (required or optional)."""
    from multifat.models.registry import declared_param_names

    return declared_param_names()


def _check_material(fname: str, m: Dict[str, Any], rep: ValidationReport) -> None:
    """Plausibility warnings on Summary values (never raises)."""
    if m["E_MPa"] is None or m["nu"] is None:
        rep.warnings.append(
            "Summary: E or nu missing; strains are undefined and models that need them are skipped"
        )
    elif not 1e4 <= m["E_MPa"] <= 5e5:
        rep.warnings.append(f"Summary: E = {m['E_MPa']} MPa outside 1e4-5e5 (units?)")
    sy, su = m["Yield_strength_MPa"], m["Tensile_strength_MPa"]
    if sy is not None and su is not None and sy > su:
        rep.warnings.append(f"Summary: yield strength {sy} > tensile strength {su}")
    s1, t1 = m["Tensile_fatigue_limit_MPa"], m["Torsional_fatigue_limit_MPa"]
    if s1 is not None and t1 is not None and not t1 < s1:
        rep.warnings.append(f"Summary: torsional fatigue limit {t1} >= tensile limit {s1}")
    names = {
        "Yield_strength_MPa": "Sy",
        "Tensile_strength_MPa": "Su",
        "Tensile_fatigue_limit_MPa": "Sigm1",
        "Torsional_fatigue_limit_MPa": "Taum1",
        "k_FS": "k_FS",
        "k_FI": "k_FI",
    }
    missing = [v for k, v in names.items() if m.get(k) is None]
    if missing:
        rep.info.append(
            "Summary values missing (models needing them are skipped): " + ", ".join(missing)
        )


#### Extract experimental data from Excel files
def import_experimental_data_grouped_by_material(
    folder: str | Path = "Experimental_Data", verbose: bool = False, overlay: bool | None = None
) -> Dict[str, Any]:
    """
    Import experimental datasets grouped by material (= Excel filename stem).

    Output format (UNCHANGED structure):
      data[material_name] = {
        "material": {...},   # now includes additional material properties
        "cases": {
          case_id: {
            "meta": {"Nf_cycles": float, "row_index": int, "source_file": str},
            node: {
              loadstep: {"stress": 3x3 np.ndarray, "strain": 3x3 np.ndarray}
            }
          },
          ...
        }
      }

    Conventions:
      - node = 1 for plain specimens
      - loadstep 1 = max, loadstep 2 = min
      - Filter Dataset: Notched == No; Run out == No gives the failures (cases),
        Run out == Yes the run-outs (kept apart, never used in fits)
      - Dataset numeric Nulls -> 0.0
      - Most material properties come from the named "Average" row; revised R = -1 fatigue
        limits come from compatible Summary source rows.
    """
    folder = Path(folder)
    files = workbook_files(folder)
    if not files:
        raise FileNotFoundError(f"No .xlsx files found in folder: {folder.resolve()}")
    data: Dict[str, Any] = {}
    for fp in files:
        block = load_material_workbook(fp, verbose=verbose, overlay=overlay)
        for m in block:
            if m in data:
                raise WorkbookError(
                    fp.name, "", "", f"material ID '{m}' already used by another file"
                )
        data.update(block)
    return data


def _alloy_token(name: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).split("(")[0].lower())


def _summary_sources(ws, header_row: int, avg_row: int) -> list[dict]:
    """Source rows of the Summary (between the header and the Average row), by header."""
    heads = [str(c.value).strip() if c.value is not None else "" for c in ws[header_row]]
    rows = []
    for row_number, values in enumerate(
        ws.iter_rows(min_row=header_row + 1, max_row=avg_row - 1, values_only=True),
        start=header_row + 1,
    ):
        if len(values) > 3 and values[3] is not None and str(values[3]).strip().lower() != "null":
            rows.append({"_row": row_number, **{h: v for h, v in zip(heads, values) if h}})
    return rows


def _summary_fatigue_limits(
    material: str, ws, header_row: int, avg_row: int
) -> tuple[list[str], dict[str, str], dict[str, float], dict[str, list[dict]]]:
    """Select R = -1 amplitudes from revised Summary rows, retaining every contributing source."""
    headers = {norm_text(c.value): str(c.value).strip() for c in ws[header_row] if c.value}
    if "stress ratio r" not in headers:
        return [], {}, {}, {}
    ratio_column = headers["stress ratio r"]
    material_column = headers.get("material name", "Material name")
    cycles_column = next(
        (name for key, name in headers.items() if key.startswith("number of cycles")), None
    )
    sources = _summary_sources(ws, header_row, avg_row)
    declared = [str(r.get(material_column, "")).strip() for r in sources]
    identity: list[str] = []
    if any(_alloy_token(n) != _alloy_token(material) for n in declared):
        identity.append(
            f"file material {material} conflicts with Summary names {declared}; "
            "confirm identity before pooling"
        )
    issues: dict[str, str] = {}
    selected: dict[str, float] = {}
    provenance: dict[str, list[dict]] = {}
    for key, label in (
        ("Sigm1", "tensile fatigue limit"),
        ("Taum1", "torsional fatigue limit"),
    ):
        column = next((name for head, name in headers.items() if label in head), None)
        if column is None:
            continue
        present = []
        compatible = []
        for row in sources:
            raw_value = row.get(column)
            value = to_float_or_none(raw_value)
            if value is None:
                continue
            present.append(row)
            ratio = to_float_or_none(row.get(ratio_column))
            if (
                math.isfinite(value)
                and value > 0
                and not (
                    isinstance(raw_value, str)
                    and raw_value.lstrip().startswith((">", "<", "\u2265", "\u2264"))
                )
                and ratio is not None
                and math.isfinite(ratio)
                and math.isclose(ratio, -1.0, rel_tol=0.0, abs_tol=1e-6)
                and _alloy_token(row.get(material_column)) == _alloy_token(material)
            ):
                compatible.append(
                    {
                        "excel_row": row["_row"],
                        "value_MPa": value,
                        "stress_ratio": ratio,
                        "reference_cycles": to_float_or_none(row.get(cycles_column)),
                        "material_name": row.get(material_column),
                        "reference": row.get("Article DOI / Reference"),
                        "authors": row.get("Authors"),
                    }
                )
        if compatible:
            selected[key] = math.fsum(item["value_MPa"] for item in compatible) / len(compatible)
            provenance[key] = compatible
        elif present:
            issues[key] = (
                "no compatible R = -1 source fatigue limit for this material; "
                "the Summary average is not used"
            )
    return identity, issues, selected, provenance


def load_material_workbook(
    fp: str | Path,
    verbose: bool = False,
    report: ValidationReport | None = None,
    overlay: bool | None = None,
) -> Dict[str, Any]:
    """Import one workbook; returns {material_name: block} as the folder importer does.

    Contract violations raise WorkbookError (see docs/workbook-format.md); soft issues are
    recorded in `report` when given. With `overlay` the declared corrections of the file
    (imp/overlays.py) are applied in memory; `overlay=False` gives the workbook as supplied.
    """
    fp = Path(fp)
    material_name = material_id(fp.stem)  # file stem without an ordering prefix
    rep = report if report is not None else ValidationReport(file=fp.name)
    rep.material = material_name

    wb = load_workbook(fp, data_only=True)
    for sheet in ("Summary", "Dataset"):
        if sheet not in wb.sheetnames:
            raise WorkbookError(fp.name, sheet, "", f"missing sheet '{sheet}'")
    if overlay is None:
        overlay = OVERLAYS_ENABLED
    overlay_id = apply_overlay_to_workbook(wb, fp) if overlay else None
    dataset_source: Any = fp
    if overlay_id is not None:
        rep.info.append(
            f"assumptions applied: {overlay_id['name']} ({overlay_id['operations']} operations)"
        )
        buffer = io.BytesIO()
        wb.save(buffer)
        buffer.seek(0)
        dataset_source = buffer
    ws = wb["Summary"]

    # Summary columns by header substring; only E and nu are required at load time.
    col_map: Dict[str, Optional[int]] = {}
    header_row: Optional[int] = None
    for key, (header_substring, required) in SUMMARY_FIELDS.items():
        found = _find_header_row_and_col_indices(ws, header_substring, search_rows=80)
        if found is None:
            if required:
                raise WorkbookError(
                    fp.name,
                    "Summary",
                    "header",
                    f"no column header containing '{header_substring}'",
                )
            col_map[key] = None
            continue
        r, c = found
        col_map[key] = c
        header_row = r if header_row is None else max(header_row, r)

    avg_loc = _get_named_cell_location(wb, "Average")
    avg_row = None
    if avg_loc is not None:
        avg_sheet, avg_cell = avg_loc
        ws_avg = wb[avg_sheet]
        cell_obj = ws_avg[avg_cell]
        val_num = to_float_or_none(cell_obj.value)
        avg_row = (
            int(round(val_num)) if (val_num is not None and val_num > 0) else int(cell_obj.row)
        )
    else:
        avg_row = _find_row_with_text(ws, "Average", search_rows=300)

    if avg_row is None:
        raise WorkbookError(fp.name, "Summary", "", "no row labelled 'Average' (material values)")

    def _value(key: str) -> Optional[float]:
        c = col_map.get(key)
        return None if c is None else to_float_or_none(ws.cell(row=avg_row, column=c).value)

    # --- Extract material properties from Summary at avg_row ---
    E_file = _value("E_MPa")
    nu_file = _value("nu")

    # Missing E or nu leaves strains undefined (stress-only models still run); a supplied
    # implausible value is a contract violation.
    if E_file is not None and E_file <= 0:
        raise WorkbookError(
            fp.name, "Summary", f"row {avg_row}", f"Young's modulus is not positive (got {E_file})"
        )
    if nu_file is not None and not (0.0 < nu_file < 0.6):
        raise WorkbookError(
            fp.name,
            "Summary",
            f"row {avg_row}",
            f"Poisson's ratio is missing/invalid or not plausible (got {nu_file})",
        )

    # Optional constants: whole-header match (units suffix ignored), see schema.SUMMARY_CONSTANTS.
    headers = [cell.value for cell in ws[header_row]] if header_row is not None else []
    const_cols = match_constant_columns(headers, SUMMARY_CONSTANTS)

    def _read_opt(key: str) -> Optional[float]:
        c = const_cols.get(key)
        return None if c is None else to_float_or_none(ws.cell(row=avg_row, column=c).value)

    # Any other column named after a constant declared by a registered model.
    known = set(SUMMARY_CONSTANTS) | {"E", "nu", "Sy", "Su", "Sigm1", "Taum1", "k_FS", "k_FI"}
    extra_names = {n: (n.lower(),) for n in _declared_param_names() - known}
    extra = {}
    for key, c in match_constant_columns(headers, extra_names).items():
        v = to_float_or_none(ws.cell(row=avg_row, column=c).value)
        if v is not None:
            extra[key] = v

    if verbose:
        print(f"[{fp.name}] material='{material_name}' -> E={E_file} MPa, nu={nu_file}")

    # Dataset
    df = pd.read_excel(dataset_source, sheet_name="Dataset", engine="openpyxl")

    cols: Dict[str, Optional[str]] = {}
    for key, tokens in DATASET_REQUIRED.items():
        try:
            cols[key] = find_col(list(df.columns), tokens)
        except KeyError:
            raise WorkbookError(
                fp.name, "Dataset", "header", f"no column whose header contains all of {tokens}"
            ) from None
    # Optional out-of-phase columns: phases in degrees with respect to axial.
    # If not present in a workbook, we fall back to 0° (in-phase).
    for key, tokens in DATASET_OPTIONAL.items():
        try:
            cols[key] = find_col(list(df.columns), tokens)
        except KeyError:
            cols[key] = None
    col_Nf, col_notch, col_run = cols["Nf"], cols["notched"], cols["runout"]
    col_mr, col_ar, col_max, col_aax = (
        cols[k] for k in ("mean_rad", "amp_rad", "mean_ax", "amp_ax")
    )
    col_ms, col_as = cols["mean_sh"], cols["amp_sh"]
    col_pr, col_pt, col_pa = cols["phase_rad"], cols["phase_tor"], cols["phase_ax"]
    revised = col_pa is not None
    freq_cols = {
        "radial": cols["freq_rad"],
        "axial": cols["freq_ax"],
        "torsional": cols["freq_tor"],
    }

    material: Dict[str, Any] = {
        "E_MPa": E_file,
        "nu": nu_file,
        "Yield_strength_MPa": _value("Yield_strength_MPa"),
        "Tensile_strength_MPa": _value("Tensile_strength_MPa"),
        "Tensile_fatigue_limit_MPa": _value("Tensile_fatigue_limit_MPa"),
        "Torsional_fatigue_limit_MPa": _value("Torsional_fatigue_limit_MPa"),
        "Tensile_fatigue_R0_MPa": _read_opt("Sig0"),
        "Torsional_fatigue_R0_MPa": _read_opt("Tau0"),
        "k_FS": _value("k_FS"),
        "k_FI": _value("k_FI"),
    }
    # Optional fatigue properties and model constants
    for key in SUMMARY_CONSTANTS:
        if key not in ("Sig0", "Tau0"):
            material[key] = _read_opt(key)
    material["extra"] = extra
    identity_issues, property_issues, selected_limits, limit_sources = _summary_fatigue_limits(
        material_name, ws, header_row, avg_row
    )
    limit_keys = {
        "Sigm1": "Tensile_fatigue_limit_MPa",
        "Taum1": "Torsional_fatigue_limit_MPa",
    }
    revised_summary = any(norm_text(c.value) == "stress ratio r" for c in ws[header_row])
    for key, field in limit_keys.items():
        if key in selected_limits:
            material[field] = selected_limits[key]
            rows = ", ".join(str(item["excel_row"]) for item in limit_sources[key])
            rep.info.append(f"Summary: {key} is the mean of compatible R = -1 source rows {rows}")
        else:
            if revised_summary and material[field] is not None:
                property_issues.setdefault(
                    key, "no usable R = -1 source fatigue limit; the Summary average is not used"
                )
            if key in property_issues:
                material["flagged_" + field] = material[field]
                material[field] = None
    _check_material(fp.name, material, rep)
    rep.warnings += [f"Summary: {t}" for t in identity_issues]
    rep.warnings += [f"Summary: {k}: {t}" for k, t in property_issues.items()]
    data: Dict[str, Any] = {
        material_name: {
            "material": material,
            "cases": {},
            "runouts": {},
            "source": {
                "workbook": fp.as_posix(),
                "format": "revised" if revised else "legacy",
                "summary_average_row": avg_row,
            },
            "identity_issues": identity_issues,
            "property_issues": property_issues,
            "fatigue_limit_sources": limit_sources,
            "excluded_observations": [],
            "overlay": overlay_id,
        }
    }
    excluded = data[material_name]["excluded_observations"]

    dropped = dict.fromkeys(
        ["notched", "run-out", "unrecognised Notched", "unrecognised Run out", "Nf <= 0"], 0
    )
    coerced = 0
    stress_cols = [col_mr, col_ar, col_max, col_aax, col_ms, col_as]
    case_id = 0
    runout_id = 0
    for idx, row in df.iterrows():
        if not is_no_value(row[col_notch]):
            dropped["notched" if is_yes_value(row[col_notch]) else "unrecognised Notched"] += 1
            continue
        is_runout = is_yes_value(row[col_run])
        if not (is_runout or is_no_value(row[col_run])):
            dropped["unrecognised Run out"] += 1
            continue

        Nf = to_float_or_zero(row[col_Nf])
        if Nf <= 0:
            dropped["Nf <= 0"] += 1
            continue

        raw_specimen = row[cols["specimen"]] if cols["specimen"] is not None else None
        specimen = "" if is_nullish(raw_specimen) else str(raw_specimen).strip()
        status = "runout" if is_runout else "failure"

        def exclude(reason: str) -> None:
            excluded.append(
                {
                    "excel_row": int(idx) + 2,
                    "specimen_id": specimen,
                    "reason": reason,
                    "test_status": status,
                }
            )

        n_empty = sum(to_float_or_none(row[c]) is None for c in stress_cols)
        if revised and n_empty:
            # Revised format: an empty stress cell is unknown loading, not an inactive component.
            exclude("missing stress component")
            continue
        coerced += n_empty
        mr = to_float_or_zero(row[col_mr])
        ar = to_float_or_zero(row[col_ar])
        ma = to_float_or_zero(row[col_max])
        aa = to_float_or_zero(row[col_aax])
        ms = to_float_or_zero(row[col_ms])
        aS = to_float_or_zero(row[col_as])
        if min(ar, aa, aS) < 0:
            rep.warnings.append(f"Dataset row {int(idx) + 2}: negative stress amplitude")

        # -----------------------------
        # Out-of-phase support (single-harmonic sinusoidal loading)
        #
        # Excel provides phase shifts (deg) for radial and torsional loads
        # w.r.t. axial. We store a compact *harmonic* representation:
        #   S(t) = S_mean + S_sin*sin(ωt) + S_cos*cos(ωt)
        # and the same for strain (linear elastic here).
        #
        # This preserves the existing 2-loadstep tensors (sigma_max/min)
        # for backward compatibility, while enabling an evaluator-side
        # correction for shear amplitude on any candidate plane.
        # -----------------------------
        phi_r_deg = to_float_or_zero(row[col_pr]) if col_pr is not None else 0.0
        phi_t_deg = to_float_or_zero(row[col_pt]) if col_pt is not None else 0.0
        if col_pa is not None:  # phases are relative to the axial load
            phi_a_deg = to_float_or_zero(row[col_pa])
            phi_r_deg -= phi_a_deg
            phi_t_deg -= phi_a_deg

        freq_text = {
            k: ("" if is_nullish(row[c]) else str(row[c]))
            for k, c in freq_cols.items()
            if c is not None
        }
        active = {"radial": ar, "axial": aa, "torsional": aS}
        known = [
            to_float_or_none(row[c])
            for k, c in freq_cols.items()
            if c is not None and active[k] != 0 and to_float_or_none(row[c]) is not None
        ]
        if known and not np.allclose(known, known[0]):
            exclude("different loading frequencies on active components (single-harmonic model)")
            continue
        p_mean = to_float_or_zero(row[cols["mean_p"]]) if cols["mean_p"] is not None else 0.0
        p_amp = to_float_or_zero(row[cols["amp_p"]]) if cols["amp_p"] is not None else 0.0
        if (p_mean or p_amp) and not (mr or ar):
            exclude("internal pressure not converted to hoop stress")
            continue

        phi_r = math.radians(float(phi_r_deg))
        phi_t = math.radians(float(phi_t_deg))

        sigma_mean = build_stress_tensor(mr, ma, ms)

        # Decompose each sinusoid with phase shift into sin/cos parts:
        #   A*sin(ωt+φ) = (A*cosφ)*sin(ωt) + (A*sinφ)*cos(ωt)
        sigma_sin = build_stress_tensor(ar * math.cos(phi_r), aa, aS * math.cos(phi_t))
        sigma_cos = build_stress_tensor(ar * math.sin(phi_r), 0.0, aS * math.sin(phi_t))

        eps_mean = stress_to_strain_iso(sigma_mean, E_file, nu_file)
        eps_sin = stress_to_strain_iso(sigma_sin, E_file, nu_file)
        eps_cos = stress_to_strain_iso(sigma_cos, E_file, nu_file)

        # max/min
        sigma_max = build_stress_tensor(mr + ar, ma + aa, ms + aS)
        sigma_min = build_stress_tensor(mr - ar, ma - aa, ms - aS)

        eps_max = stress_to_strain_iso(sigma_max, E_file, nu_file)
        eps_min = stress_to_strain_iso(sigma_min, E_file, nu_file)

        block = {
            "meta": {
                "Nf_cycles": float(Nf),
                "row_index": int(idx),
                "source_file": fp.name,
                "phase_rad_deg": float(phi_r_deg),
                "phase_tor_deg": float(phi_t_deg),
                "test_status": status,
                "specimen_id": specimen,
                "excel_row": int(idx) + 2,
                "frequency_text": freq_text,
            },
            1: {
                1: {"stress": sigma_max, "strain": eps_max},
                2: {"stress": sigma_min, "strain": eps_min},
            },
            "harmonics": {
                "S_mean": sigma_mean,
                "S_sin": sigma_sin,
                "S_cos": sigma_cos,
                "E_mean": eps_mean,
                "E_sin": eps_sin,
                "E_cos": eps_cos,
            },
        }
        if is_runout:
            data[material_name]["runouts"][runout_id] = block
            runout_id += 1
            continue
        data[material_name]["cases"][case_id] = block
        case_id += 1

    rep.info.append(f"Dataset rows: {len(df)}, kept: {case_id} failures, {runout_id} run-outs")
    for reason, n in dropped.items():
        if n:
            (rep.warnings if reason.startswith("unrecognised") else rep.info).append(
                f"rows dropped ({reason}): {n}"
            )
    if coerced:
        rep.info.append(f"empty/Null stress cells read as 0: {coerced}")
    for reason in dict.fromkeys(e["reason"] for e in excluded):
        n = sum(e["reason"] == reason for e in excluded)
        rep.warnings.append(f"rows excluded ({reason}): {n}")
    if case_id == 0 and runout_id == 0:
        rep.errors.append("no usable test rows (plain, Nf > 0)")
    elif case_id == 0:
        rep.info.append(f"no plain failures: {runout_id} run-outs only")

    if verbose:
        print(f"[{material_name}] cases kept: {case_id}")

    return data


#### Extract tensors
def extract_tensor(
    data: Dict[str, Any],
    material: str,
    case_id: int,
    node: int = 1,
    include_harmonics: bool = False,
) -> (
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    | Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]
):
    """
    Extract tensors from the nested experimental-data dict.

    Convention (legacy / kept for backward compatibility):
      loadstep 1 = max (step 1)
      loadstep 2 = min (step 2)

    Out-of-phase support
    -------------------
    If ``include_harmonics=True``, this also returns a 5th element: a dict
    with a compact single-harmonic representation:

      S(t) = S_mean + S_sin*sin(ωt) + S_cos*cos(ωt)
      E(t) = E_mean + E_sin*sin(ωt) + E_cos*cos(ωt)

    When phase columns are not present in the Excel sheet, the importer stores
    a default in-phase cycle (S_cos=0, E_cos=0).
    """
    case = data[material]["cases"][case_id]

    S0 = np.asarray(case[node][1]["stress"], dtype=float)
    E0 = np.asarray(case[node][1]["strain"], dtype=float)
    S1 = np.asarray(case[node][2]["stress"], dtype=float)
    E1 = np.asarray(case[node][2]["strain"], dtype=float)

    # Optional sanity check
    for name, A in [("S0", S0), ("E0", E0), ("S1", S1), ("E1", E1)]:
        if A.shape != (3, 3):
            raise ValueError(f"{name} must be (3,3), got {A.shape}")

    if not include_harmonics:
        return S0, E0, S1, E1

    harmonics = case.get("harmonics", None)
    if harmonics is None:
        z = np.zeros((3, 3), dtype=float)
        harmonics = {
            "S_mean": 0.5 * (S0 + S1),
            "S_sin": 0.5 * (S0 - S1),
            "S_cos": z.copy(),
            "E_mean": 0.5 * (E0 + E1),
            "E_sin": 0.5 * (E0 - E1),
            "E_cos": z.copy(),
        }

    return S0, E0, S1, E1, harmonics


#### Extract material data
def get_material_properties(data: Dict[str, Any], material: str) -> Dict[str, Any]:
    """
    Return the material properties dictionary for a given material name.

    This keeps extract_tensor working unchanged, because the structure of data is preserved:
      data[material]["material"] is still where material properties live.
    """
    if material not in data:
        raise KeyError(f"Material '{material}' not found. Available: {list(data.keys())}")
    return data[material]["material"]


#### Extract sub-dataset for uniaxial tension with R=-1
def extract_uniaxial_tension(
    data: Dict[str, Any],
    *,
    node: int = 1,
    rtol_R: float = 1e-3,
    atol_R: float = 1e-6,
    rtol_uniax: float = 1e-3,
    atol_uniax: float = 1e-6,
) -> Dict[str, Any]:
    """
    Extract experimental points corresponding to *uniaxial tensile loading* with *R = -1*.

    This works solely from the already-stored max/min stress tensors in `data` and DOES NOT
    change the data structure.

    Assumptions:
      - AXIAL stress is stored in sigma_yy (stress[1,1])
      - loadstep 1 = max, loadstep 2 = min

    Selection criteria (per case):
      1) R_axial = sigma_yy(min) / sigma_yy(max) ≈ -1
      2) Uniaxial check: all stress components except sigma_yy are near zero
      3) sigma_yy(max) > 0 (tension in the maximum load state)

    Parameters
    ----------
    data : dict
        Output of import_experimental_data_grouped_by_material(...)
    node : int, default 1
        Node number (plain specimen -> 1).
    rtol_R, atol_R : float
        Relative/absolute tolerances for the R=-1 condition.
    rtol_uniax, atol_uniax : float
        Relative/absolute tolerances for the uniaxial condition.

    Returns
    -------
    out : dict
        {
          material_name: {
            "material": <material properties dict>,
            "points": [
              {
                "case_id": int,
                "Nf_cycles": float,
                "sigma_max_MPa": float,
                "sigma_min_MPa": float,
                "sigma_a_MPa": float,
                "eps_a": float,
              },
              ...
            ]
          },
          ...
        }
    """
    out: Dict[str, Any] = {}

    for material_name, mat_block in data.items():
        mat_props = mat_block.get("material", {})
        cases = mat_block.get("cases", {})

        points: List[Dict[str, float]] = []

        for case_id, case in cases.items():
            # Extract tensors directly from the existing structure
            Smax = np.asarray(case[node][1]["stress"], dtype=float)  # loadstep 1: max
            Emin = np.asarray(
                case[node][2]["strain"], dtype=float
            )  # not used directly, but kept for symmetry
            Smin = np.asarray(case[node][2]["stress"], dtype=float)  # loadstep 2: min
            Emax = np.asarray(case[node][1]["strain"], dtype=float)

            # Axial component (AXIAL -> sigma_yy)
            smax = float(Smax[1, 1])
            smin = float(Smin[1, 1])

            # Need a nonzero max to define R robustly
            if abs(smax) < atol_R:
                continue

            # Must be tensile in maximum state
            if smax <= 0.0:
                continue

            # Check R ≈ -1
            R_ax = smin / smax
            if not np.isclose(R_ax, -1.0, rtol=rtol_R, atol=atol_R):
                continue

            # Uniaxial check: all other stress components near zero relative to |smax|
            S_other = Smax.copy()
            S_other[1, 1] = 0.0  # remove axial component
            max_other = float(np.max(np.abs(S_other)))
            if max_other > (rtol_uniax * abs(smax) + atol_uniax):
                continue

            # Build amplitudes (general definition; for R=-1 this will be ~smax)
            sigma_a = 0.5 * (smax - smin)

            # Axial strain amplitude from the already computed strain tensors
            emax = float(Emax[1, 1])
            emin = float(Emin[1, 1])
            eps_a = 0.5 * (emax - emin)

            Nf = float(case["meta"]["Nf_cycles"])

            points.append(
                {
                    "case_id": int(case_id),
                    "Nf_cycles": Nf,
                    "sigma_max_MPa": smax,
                    "sigma_min_MPa": smin,
                    "sigma_a_MPa": float(sigma_a),
                    "eps_a": float(eps_a),
                }
            )

        if points:
            out[material_name] = {
                "material": mat_props,  # keep all material properties for convenience
                "points": points,
            }

    return out


def exclude_subdataset(
    data: Dict[str, Any],
    excluded: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Return a new dataset with the cases listed in `excluded` removed from `data`.

    The sub-datasets come from extraction functions such as
    `extract_uniaxial_tension_Rminus1`, which returns
    ``excluded[material_name]["points"] = [{"case_id": int, ...}, ...]``.

    The returned dataset keeps the structure of `data`, so `extract_tensor` and
    the other utilities are unchanged, and a material with no case left is kept
    as an empty block.

    Parameters
    ----------
    data : dict
        Output of import_experimental_data_grouped_by_material(...)
    excluded : dict
        Output of an extraction function, e.g. extract_uniaxial_tension_Rminus1(...)

    Returns
    -------
    filtered_data : dict
        Same structure as `data`, but with excluded case_ids removed per material.
    """
    # Build a lookup: material -> set(case_ids)
    to_remove: Dict[str, Set[int]] = {}
    for material_name, blk in excluded.items():
        pts = blk.get("points", [])
        ids = set()
        for p in pts:
            if "case_id" in p:
                ids.add(int(p["case_id"]))
        if ids:
            to_remove[material_name] = ids

    # Copy while filtering cases
    filtered: Dict[str, Any] = {}
    for material_name, mat_block in data.items():
        mat_props = mat_block.get("material", {})
        cases = mat_block.get("cases", {})

        remove_ids = to_remove.get(material_name, set())

        new_cases = {}
        for cid, case in cases.items():
            if int(cid) in remove_ids:
                continue
            new_cases[int(cid)] = case

        filtered[material_name] = {
            "material": mat_props,
            "cases": new_cases,
        }

    return filtered

