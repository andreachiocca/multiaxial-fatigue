from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set
import numpy as np
import pandas as pd
import math
from openpyxl import load_workbook


#### Extract experimental data from Excel files
def import_experimental_data_grouped_by_material(
    folder: str | Path = "Experimental_Data",
    verbose: bool = False
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
      - Filter Dataset: Notched == No and Run out == No
      - Dataset numeric Nulls -> 0.0
      - Material properties are read from "Summary" at the row number given by the named cell "Average"
    """
    folder = Path(folder)
    files = sorted(folder.glob("*.xlsx"))
    if not files:
        raise FileNotFoundError(f"No .xlsx files found in folder: {folder.resolve()}")

    # -----------------------------
    # Helpers
    # -----------------------------
    def norm_text(s: Any) -> str:
        return " ".join(str(s).strip().split()).lower()

    def is_nullish(x: Any) -> bool:
        if x is None:
            return True
        if isinstance(x, float) and np.isnan(x):
            return True
        if isinstance(x, (np.floating,)):
            return bool(np.isnan(float(x)))
        if isinstance(x, str):
            t = x.strip().lower()
            return t in {"", "null", "nan", "na", "n/a", "-", "—"}
        return False

    def _parse_numeric_string(s: str) -> Tuple[str, bool]:
        s = s.strip()
        if s == "":
            return "", False

        if s.startswith((">", "<", "≥", "≤")):
            s = s[1:].strip()

        is_percent = False
        if "%" in s:
            s = s.replace("%", "")
            is_percent = True

        s = s.replace("\u00A0", "").replace(" ", "")

        if "," in s and "." in s:
            last_comma = s.rfind(",")
            last_dot = s.rfind(".")
            if last_comma > last_dot:
                s = s.replace(".", "")
                s = s.replace(",", ".")
            else:
                s = s.replace(",", "")
        elif "," in s and "." not in s:
            parts = s.split(",")
            if len(parts) == 2:
                left, right = parts[0], parts[1]
                if "e" in right.lower():
                    s = left + "." + right
                else:
                    if left == "0" or len(right) != 3:
                        s = left + "." + right
                    else:
                        s = left + right
            else:
                s = s.replace(",", "")

        return s, is_percent

    def to_float_or_none(x: Any) -> Optional[float]:
        if is_nullish(x):
            return None
        if isinstance(x, (int, float, np.integer, np.floating)) and not (isinstance(x, float) and np.isnan(x)):
            return float(x)

        s_norm, is_percent = _parse_numeric_string(str(x))
        if s_norm == "":
            return None

        try:
            val = float(s_norm)
        except Exception:
            v = pd.to_numeric(pd.Series([s_norm]), errors="coerce").iloc[0]
            if isinstance(v, float) and np.isnan(v):
                return None
            val = float(v)

        if is_percent:
            val /= 100.0
        return float(val)

    def to_float_or_zero(x: Any) -> float:
        v = to_float_or_none(x)
        return 0.0 if v is None else float(v)

    def find_col(columns: List[str], must_contain: List[str]) -> str:
        ncols = [(c, norm_text(c)) for c in columns]
        for original, normalized in ncols:
            if all(token.lower() in normalized for token in must_contain):
                return original
        raise KeyError(f"Could not find column tokens {must_contain}. Available: {list(columns)}")

    def is_no_value(x: Any) -> bool:
        if is_nullish(x):
            return False
        if isinstance(x, (bool, np.bool_)):
            return (x is False)
        if isinstance(x, (int, np.integer)):
            return int(x) == 0
        return str(x).strip().lower().startswith("no")

    def build_stress_tensor(sig_rad: float, sig_ax: float, tau: float) -> np.ndarray:
        return np.array([[sig_rad, tau,    0.0],
                         [tau,    sig_ax, 0.0],
                         [0.0,    0.0,    0.0]], dtype=float)

    def stress_to_strain_iso(stress: np.ndarray, E: float, nu: float) -> np.ndarray:
        sig_xx, sig_yy, sig_zz = stress[0, 0], stress[1, 1], stress[2, 2]
        tau_xy, tau_yz, tau_zx = stress[0, 1], stress[1, 2], stress[2, 0]

        G = E / (2.0 * (1.0 + nu))

        eps_xx = (sig_xx - nu * (sig_yy + sig_zz)) / E
        eps_yy = (sig_yy - nu * (sig_xx + sig_zz)) / E
        eps_zz = (sig_zz - nu * (sig_xx + sig_yy)) / E

        eps_xy = tau_xy / (2.0 * G)
        eps_yz = tau_yz / (2.0 * G)
        eps_zx = tau_zx / (2.0 * G)

        return np.array([[eps_xx, eps_xy, eps_zx],
                         [eps_xy, eps_yy, eps_yz],
                         [eps_zx, eps_yz, eps_zz]], dtype=float)

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

    def _find_header_row_and_col_indices(ws, header_contains: str, search_rows: int = 80) -> Tuple[int, int]:
        """
        Find the header row and column index where a cell contains `header_contains` (case-insensitive).
        Returns (header_row, col_idx). Raises if not found.
        """
        target = header_contains.lower()
        for r in range(1, search_rows + 1):
            for c_idx, cell in enumerate(ws[r], start=1):
                if target in norm_text(cell.value):
                    return r, c_idx
        raise KeyError(f"Could not find header containing '{header_contains}' in Summary sheet.")

    def _find_header_row_and_cols(ws) -> Tuple[int, Dict[str, int]]:
        """
        Build a mapping from desired property names -> column index, for Summary sheet.
        We search for headers by "contains" matching, robust to minor formatting differences.
        Returns (header_row_used, col_index_map).
        """
        # Desired properties and the header substrings we search for
        desired = {
            "E_MPa": "Young's modulus",
            "nu": "Poisson's ratio",
            "Yield_strength_MPa": "Yield strength",
            "Tensile_strength_MPa": "Tensile strength",
            "Tensile_fatigue_limit_MPa": "Tensile fatigue limit",
            "Torsional_fatigue_limit_MPa": "Torsional fatigue limit",
            "k_FS": "k FS",
            "k_FI": "k FI",
        }

        col_map: Dict[str, int] = {}
        header_row_used: Optional[int] = None

        # Find each column independently; keep the max header row found as "header_row_used"
        for key, header_substring in desired.items():
            r, c = _find_header_row_and_col_indices(ws, header_substring, search_rows=80)
            col_map[key] = c
            header_row_used = r if header_row_used is None else max(header_row_used, r)

        if header_row_used is None:
            raise KeyError("Could not locate Summary headers.")
        return header_row_used, col_map

    # Dataset column token map
    colmap_tokens = {
        "Nf":       ["number", "cycles", "failure"],
        "notched":  ["notched"],
        "runout":   ["run", "out"],
        "mean_rad": ["mean", "radial", "stress"],
        "amp_rad":  ["amplit", "radial", "stress"],
        "mean_ax":  ["mean", "axial", "stress"],
        "amp_ax":   ["amplit", "axial", "stress"],
        "mean_sh":  ["mean", "shear", "stress"],
        "amp_sh":   ["amplit", "shear", "stress"],

        # Optional out-of-phase columns. If absent, phases default to 0.
        "phase_rad": ["phase", "radial", "load"],
        "phase_tor": ["phase", "torsional", "load"],
    }

    # -----------------------------
    # Main
    # -----------------------------
    data: Dict[str, Any] = {}

    for fp in files:
        material_name = fp.stem  # filename without extension

        wb = load_workbook(fp, data_only=True)
        if "Summary" not in wb.sheetnames:
            raise KeyError(f"[{fp.name}] Missing sheet 'Summary'.")
        ws = wb["Summary"]

        header_row, col_map = _find_header_row_and_cols(ws)

        # ------------------------------------------------------------
        # Optional material properties
        #
        # These are NOT required for the baseline workflow, but enable
        # direct use of some research methods without extra user code.
        # If they are absent we keep them as None and downstream code
        # may estimate/default them (emitting warnings).
        #
        # Expected header names (case-insensitive substring match):
        #   sigma_f, epsilon_f, b_fat, c_fat
        #   tau_f, gamma_f, b0, c0
        #   k_MGSE, a_MSWT, alpha_w, xi_ZHU, L_LI
        # ------------------------------------------------------------
        def _optional_col_idx(*needles: str) -> Optional[int]:
            row = ws[header_row]
            for c_idx, cell in enumerate(row, start=1):
                txt = norm_text(cell.value)
                for n in needles:
                    if n.lower() in txt:
                        return c_idx
            return None

        avg_loc = _get_named_cell_location(wb, "Average")
        avg_row = None
        if avg_loc is not None:
            avg_sheet, avg_cell = avg_loc
            ws_avg = wb[avg_sheet]
            cell_obj = ws_avg[avg_cell]
            val_num = to_float_or_none(cell_obj.value)
            avg_row = int(round(val_num)) if (val_num is not None and val_num > 0) else int(cell_obj.row)
        else:
            avg_row = _find_row_with_text(ws, "Average", search_rows=300)

        if avg_row is None:
            raise KeyError(f"[{fp.name}] Could not determine the 'Average' row.")

        # --- Extract material properties from Summary at avg_row ---
        E_file = to_float_or_none(ws.cell(row=avg_row, column=col_map["E_MPa"]).value)
        nu_file = to_float_or_none(ws.cell(row=avg_row, column=col_map["nu"]).value)

        if E_file is None:
            raise ValueError(f"[{fp.name}] Young's modulus at Summary row {avg_row} is missing/invalid.")
        if nu_file is None or not (0.0 < nu_file < 0.6):
            raise ValueError(f"[{fp.name}] Poisson ratio at Summary row {avg_row} is missing/invalid or not plausible.")

        # Optional extra properties (can be None if not present / invalid)
        ys = to_float_or_none(ws.cell(row=avg_row, column=col_map["Yield_strength_MPa"]).value)
        uts = to_float_or_none(ws.cell(row=avg_row, column=col_map["Tensile_strength_MPa"]).value)
        sfl = to_float_or_none(ws.cell(row=avg_row, column=col_map["Tensile_fatigue_limit_MPa"]).value)
        tfl = to_float_or_none(ws.cell(row=avg_row, column=col_map["Torsional_fatigue_limit_MPa"]).value)
        kfs = to_float_or_none(ws.cell(row=avg_row, column=col_map["k_FS"]).value)
        kfi = to_float_or_none(ws.cell(row=avg_row, column=col_map["k_FI"]).value)

        # Optional fatigue properties / extra model constants
        def _read_opt(key: str, *needles: str) -> Optional[float]:
            c = _optional_col_idx(*needles)
            return None if c is None else to_float_or_none(ws.cell(row=avg_row, column=c).value)

        # Optional maximum stresses in zero-to-maximum (R=0) fatigue tests.
        # BP estimates missing values using Eqs. (1)-(2) of Böhme et al. (2026).
        s0 = _read_opt("Sig0", "tensile fatigue strength r=0", "tensile fatigue limit r=0",
                       "tensile fatigue strength r = 0", "tensile fatigue limit r = 0")
        t0 = _read_opt("Tau0", "torsional fatigue strength r=0", "torsional fatigue limit r=0",
                       "torsional fatigue strength r = 0", "torsional fatigue limit r = 0")

        sigma_f_opt = _read_opt("sigma_f", "sigma_f", "fatigue strength coefficient")
        epsilon_f_opt = _read_opt("epsilon_f", "epsilon_f", "fatigue ductility coefficient")
        b_fat_opt = _read_opt("b_fat", "b_fat", "fatigue strength exponent")
        c_fat_opt = _read_opt("c_fat", "c_fat", "fatigue ductility exponent")

        tau_f_opt = _read_opt("tau_f", "tau_f", "shear fatigue strength")
        gamma_f_opt = _read_opt("gamma_f", "gamma_f", "shear fatigue ductility")
        b0_opt = _read_opt("b0", "b0", "shear fatigue strength exponent")
        c0_opt = _read_opt("c0", "c0", "shear fatigue ductility exponent")

        k_mgse_opt = _read_opt("k_MGSE", "k_mgse", "k mgse")
        a_mswt_opt = _read_opt("a_MSWT", "a_mswt", "a mswt")
        alpha_w_opt = _read_opt("alpha_w", "alpha_w", "alpha w")
        xi_zhu_opt = _read_opt("xi_ZHU", "xi_zhu", "xi zhu")
        L_li_opt = _read_opt("L_LI", "l_li", "li weight")

        E_file = float(E_file)
        nu_file = float(nu_file)

        if verbose:
            print(f"[{fp.name}] material='{material_name}' -> E={E_file:.3f} MPa, nu={nu_file:.6f}")

        # Dataset
        df = pd.read_excel(fp, sheet_name="Dataset", engine="openpyxl")

        col_Nf    = find_col(list(df.columns), colmap_tokens["Nf"])
        col_notch = find_col(list(df.columns), colmap_tokens["notched"])
        col_run   = find_col(list(df.columns), colmap_tokens["runout"])

        col_mr    = find_col(list(df.columns), colmap_tokens["mean_rad"])
        col_ar    = find_col(list(df.columns), colmap_tokens["amp_rad"])
        col_max   = find_col(list(df.columns), colmap_tokens["mean_ax"])
        col_aax   = find_col(list(df.columns), colmap_tokens["amp_ax"])
        col_ms    = find_col(list(df.columns), colmap_tokens["mean_sh"])
        col_as    = find_col(list(df.columns), colmap_tokens["amp_sh"])

        # Optional out-of-phase columns: phases in degrees with respect to axial.
        # If not present in a workbook, we fall back to 0° (in-phase).
        try:
            col_pr = find_col(list(df.columns), colmap_tokens["phase_rad"])
        except KeyError:
            col_pr = None
        try:
            col_pt = find_col(list(df.columns), colmap_tokens["phase_tor"])
        except KeyError:
            col_pt = None

        # Initialize container for this material (structure preserved; only adds keys inside "material")
        data[material_name] = {
            "material": {
                "E_MPa": E_file,
                "nu": nu_file,
                "Yield_strength_MPa": ys,
                "Tensile_strength_MPa": uts,
                "Tensile_fatigue_limit_MPa": sfl,
                "Torsional_fatigue_limit_MPa": tfl,
                "Tensile_fatigue_R0_MPa": s0,
                "Torsional_fatigue_R0_MPa": t0,
                "k_FS": kfs,
                "k_FI": kfi,

                # Optional fatigue properties
                "sigma_f": sigma_f_opt,
                "epsilon_f": epsilon_f_opt,
                "b_fat": b_fat_opt,
                "c_fat": c_fat_opt,
                "tau_f": tau_f_opt,
                "gamma_f": gamma_f_opt,
                "b0": b0_opt,
                "c0": c0_opt,

                # Optional model constants
                "k_MGSE": k_mgse_opt,
                "a_MSWT": a_mswt_opt,
                "alpha_w": alpha_w_opt,
                "xi_ZHU": xi_zhu_opt,
                "L_LI": L_li_opt,
            },
            "cases": {}
        }

        case_id = 0
        for idx, row in df.iterrows():
            if not is_no_value(row[col_notch]):
                continue
            if not is_no_value(row[col_run]):
                continue

            Nf = to_float_or_zero(row[col_Nf])
            if Nf <= 0:
                continue

            mr = to_float_or_zero(row[col_mr]); ar = to_float_or_zero(row[col_ar])
            ma = to_float_or_zero(row[col_max]); aa = to_float_or_zero(row[col_aax])
            ms = to_float_or_zero(row[col_ms]);  aS = to_float_or_zero(row[col_as])

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

            data[material_name]["cases"][case_id] = {
                "meta": {
                    "Nf_cycles": float(Nf),
                    "row_index": int(idx),
                    "source_file": fp.name,
                    "phase_rad_deg": float(phi_r_deg),
                    "phase_tor_deg": float(phi_t_deg),
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
                }
            }
            case_id += 1

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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
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

    Assumptions consistent with your importer:
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
            Emin = np.asarray(case[node][2]["strain"], dtype=float)  # not used directly, but kept for symmetry
            Smin = np.asarray(case[node][2]["stress"], dtype=float)  # loadstep 2: min
            Emax = np.asarray(case[node][1]["strain"], dtype=float)

            # Axial component (based on your mapping AXIAL -> sigma_yy)
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

            points.append({
                "case_id": int(case_id),
                "Nf_cycles": Nf,
                "sigma_max_MPa": smax,
                "sigma_min_MPa": smin,
                "sigma_a_MPa": float(sigma_a),
                "eps_a": float(eps_a),
            })

        if points:
            out[material_name] = {
                "material": mat_props,   # keep all material properties for convenience
                "points": points,
            }

    return out


def exclude_subdataset(
    data: Dict[str, Any],
    excluded: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Return a NEW dataset with the cases listed in `excluded` removed from `data`.

    Designed to work with sub-datasets returned by functions like
    `extract_uniaxial_tension_Rminus1`, which returns:
      excluded[material_name]["points"] = [{"case_id": int, ...}, ...]

    IMPORTANT:
    - The returned dataset keeps the SAME structure as `data`, so `extract_tensor`
      and any other utilities keep working unchanged.
    - Materials are preserved even if they end up with zero cases (you can prune later if desired).

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
