# calibrate_cp_params.py
"""Calibration utility for critical-plane / model parameters.

Originally used to calibrate:
- k_FS for Fatemi–Socie (FS)
- k_FI for Findley (FIN)

Extended (per your request) to also calibrate additional material parameters:
- L_LI     (Li 2021 CP method)
- k_MGSE   (MGSE_ZHU CP method)
- a_MSWT   (MSWT CP method)
- alpha_w  (ZhuEDP case-level model)
- xi_ZHU   (ZhuEDP case-level model)

Repository conventions (v6)
--------------------------
- Experimental data importer: fatigue.imp.experimental_data
- Plane sampling: fatigue.planes.sampling (dcmesh_xyz, gentang)
- Critical-plane evaluation: fatigue.eval.evaluator.best_cp_for_case
- CP methods: fatigue.models.cp_methods
- Case-level model: fatigue.models.energy_based_methods.zhu_edp.ZhuEDP

Calibration objective
---------------------
For each candidate parameter value(s), compute the model scalar (CP/DP) for each
selected case, fit a power law in log-space:

    CP = A * N^b

and minimize fatigue scatter defined as the STD of residuals in ln(N) using the
project's helper: fatigue.eval.fitting.fit_power_law_with_survival_std.

Datasets
--------
- FS, FIN, LI, MGSE_ZHU, MSWT:
    use uniaxial axial + uniaxial torsion points (target R=-1 with fallbacks).
    This matches the legacy calibration intent.

- ZhuEDP (alpha_w, xi_ZHU):
    uses *multiaxial* cases where both axial and shear strain amplitudes are
    non-negligible (based on the harmonic representation). If such cases are
    not present, the script falls back to all cases (with a warning).

Notes
-----
- This script does NOT write back to Excel; it prints calibrated parameters and
  returns a nested dictionary for programmatic use.
- No SciPy is used: a grid-refine search is implemented.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from fatigue.imp import experimental_data as exp
from fatigue.planes.sampling import dcmesh_xyz, gentang
from fatigue.eval.evaluator import best_cp_for_case
from fatigue.eval.fitting import fit_power_law_with_survival_std

from fatigue.models.cp_methods.fs import FatemiSocie
from fatigue.models.cp_methods.fin import Findley
from fatigue.models.cp_methods.li import LI
from fatigue.models.cp_methods.mgse_zhu import MGSE_Zhu
from fatigue.models.cp_methods.mswt import MSWT
from fatigue.models.energy_based_methods.zhu_edp import ZhuEDP


# -----------------------------
# Simple minimizers (no SciPy)
# -----------------------------

def minimize_1d_grid_refine(
    objective: Callable[[float], float],
    bounds: Tuple[float, float],
    *,
    n_grid: int = 41,
    n_refine: int = 6,
) -> Tuple[float, float]:
    """Minimize objective(x) over [lo, hi] via repeated grid search."""
    lo, hi = float(bounds[0]), float(bounds[1])
    if not (hi > lo):
        raise ValueError("bounds must satisfy hi > lo")

    x_best = np.nan
    f_best = np.inf

    for _ in range(int(n_refine)):
        xs = np.linspace(lo, hi, int(n_grid))
        fs = np.array([objective(float(x)) for x in xs], dtype=float)

        j = int(np.nanargmin(fs))
        x_best = float(xs[j])
        f_best = float(fs[j])

        step = (hi - lo) / max(1, (len(xs) - 1))
        lo = max(bounds[0], x_best - 2.0 * step)
        hi = min(bounds[1], x_best + 2.0 * step)

        if hi - lo < 1e-10:
            break

    return x_best, f_best


def minimize_2d_grid_refine(
    objective: Callable[[float, float], float],
    bounds_x: Tuple[float, float],
    bounds_y: Tuple[float, float],
    *,
    n_grid_x: int = 21,
    n_grid_y: int = 21,
    n_refine: int = 5,
) -> Tuple[float, float, float]:
    """Minimize objective(x,y) with a coarse 2D grid + iterative window shrinking."""
    x_lo, x_hi = float(bounds_x[0]), float(bounds_x[1])
    y_lo, y_hi = float(bounds_y[0]), float(bounds_y[1])
    if not (x_hi > x_lo and y_hi > y_lo):
        raise ValueError("bounds must satisfy hi > lo for both variables")

    x_best = np.nan
    y_best = np.nan
    f_best = np.inf

    for _ in range(int(n_refine)):
        xs = np.linspace(x_lo, x_hi, int(n_grid_x))
        ys = np.linspace(y_lo, y_hi, int(n_grid_y))

        for x in xs:
            for y in ys:
                f = float(objective(float(x), float(y)))
                if np.isfinite(f) and f < f_best:
                    f_best, x_best, y_best = f, float(x), float(y)

        # shrink window around best
        x_step = (x_hi - x_lo) / max(1, (len(xs) - 1))
        y_step = (y_hi - y_lo) / max(1, (len(ys) - 1))

        x_lo = max(bounds_x[0], x_best - 2.0 * x_step)
        x_hi = min(bounds_x[1], x_best + 2.0 * x_step)
        y_lo = max(bounds_y[0], y_best - 2.0 * y_step)
        y_hi = min(bounds_y[1], y_best + 2.0 * y_step)

        if (x_hi - x_lo) < 1e-10 and (y_hi - y_lo) < 1e-10:
            break

    return float(x_best), float(y_best), float(f_best)


# -----------------------------
# Dataset selection (uniaxial tension / torsion)
# -----------------------------

def _dominance_ratio(other_max: float, main_max: float, *, atol: float, rtol: float) -> bool:
    return other_max <= (rtol * abs(main_max) + atol)


def _case_get_Smax_Smin(case: Dict[str, Any], node: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    S0 = np.asarray(case[node][1]["stress"], dtype=float)
    S1 = np.asarray(case[node][2]["stress"], dtype=float)
    return S0, S1


def _uniaxial_tension_filter(
    case: Dict[str, Any],
    *,
    node: int = 1,
    r_target: Optional[float] = -1.0,
    rtol_R: float = 0.02,
    atol_R: float = 1e-6,
    rtol_uniax: float = 1e-3,
    atol_uniax: float = 1e-6,
) -> Tuple[bool, float]:
    S0, S1 = _case_get_Smax_Smin(case, node=node)

    s0 = float(S0[1, 1])
    s1 = float(S1[1, 1])
    smax = max(s0, s1)
    smin = min(s0, s1)

    if abs(smax) < atol_R:
        return False, np.nan

    R = smin / smax

    if smax <= 0.0:
        return False, R

    S_other = S0.copy()
    S_other[1, 1] = 0.0
    other_max = float(np.max(np.abs(S_other)))
    if not _dominance_ratio(other_max, smax, atol=atol_uniax, rtol=rtol_uniax):
        return False, R

    if r_target is not None:
        if not np.isclose(R, r_target, rtol=rtol_R, atol=atol_R):
            return False, R

    return True, R


def _uniaxial_torsion_filter(
    case: Dict[str, Any],
    *,
    node: int = 1,
    r_target: Optional[float] = -1.0,
    rtol_R: float = 0.02,
    atol_R: float = 1e-6,
    rtol_uniax: float = 1e-3,
    atol_uniax: float = 1e-6,
) -> Tuple[bool, float]:
    S0, S1 = _case_get_Smax_Smin(case, node=node)

    t0 = float(S0[0, 1])
    t1 = float(S1[0, 1])
    tmax = max(t0, t1)
    tmin = min(t0, t1)

    if abs(tmax) < atol_R:
        return False, np.nan

    R = tmin / tmax

    if tmax <= 0.0:
        return False, R

    S_other = S0.copy()
    S_other[0, 1] = 0.0
    S_other[1, 0] = 0.0
    other_max = float(np.max(np.abs(S_other)))
    if not _dominance_ratio(other_max, tmax, atol=atol_uniax, rtol=rtol_uniax):
        return False, R

    if r_target is not None:
        if not np.isclose(R, r_target, rtol=rtol_R, atol=atol_R):
            return False, R

    return True, R


def _select_cases_with_fallback(
    cases: Dict[int, Dict[str, Any]],
    selector: Callable[..., Tuple[bool, float]],
    *,
    selector_kwargs: Dict[str, Any],
    fallback_rtol_list: Sequence[float] = (0.02, 0.05, 0.1, 0.2, 0.5),
    allow_any_R_fallback: bool = True,
    name_for_warnings: str = "dataset",
) -> Tuple[List[int], Dict[str, Any]]:
    info: Dict[str, Any] = {"fallback_used": False, "rtol_R_used": None, "any_R_used": False}

    selected: List[int] = []
    for cid, case in cases.items():
        ok, _R = selector(case, **selector_kwargs)
        if ok:
            selected.append(int(cid))

    if selected:
        info["rtol_R_used"] = selector_kwargs.get("rtol_R", None)
        return selected, info

    base_kwargs = dict(selector_kwargs)
    for rtol_R in fallback_rtol_list:
        base_kwargs["rtol_R"] = float(rtol_R)
        selected = []
        for cid, case in cases.items():
            ok, _R = selector(case, **base_kwargs)
            if ok:
                selected.append(int(cid))
        if selected:
            warnings.warn(
                f"[{name_for_warnings}] No points found for target R with initial tolerance; using relaxed rtol_R={rtol_R}.",
                RuntimeWarning,
            )
            info.update({"fallback_used": True, "rtol_R_used": rtol_R})
            return selected, info

    if allow_any_R_fallback:
        anyR_kwargs = dict(selector_kwargs)
        anyR_kwargs["r_target"] = None
        selected = []
        for cid, case in cases.items():
            ok, _R = selector(case, **anyR_kwargs)
            if ok:
                selected.append(int(cid))
        if selected:
            warnings.warn(
                f"[{name_for_warnings}] No points found for target R even after relaxing; falling back to ANY R.",
                RuntimeWarning,
            )
            info.update({"fallback_used": True, "any_R_used": True})
            return selected, info

    warnings.warn(f"[{name_for_warnings}] No usable points found.", RuntimeWarning)
    info.update({"fallback_used": True})
    return [], info


# -----------------------------
# Calibration helpers
# -----------------------------

@dataclass(frozen=True)
class CalibrationResult1D:
    param_name: str
    p_opt: float
    scatter_STD: float
    fit: Dict[str, float]
    series: Dict[str, np.ndarray]


@dataclass(frozen=True)
class CalibrationResult2D:
    param_names: Tuple[str, str]
    p_opt: Tuple[float, float]
    scatter_STD: float
    fit: Dict[str, float]
    series: Dict[str, np.ndarray]


def _build_R_list(*, n1: int = 20, n2: int = 40) -> np.ndarray:
    xyz = dcmesh_xyz(a1=np.pi, n1=int(n1), a2=2 * np.pi, n2=int(n2))
    xyz = np.unique(np.round(xyz, 12), axis=0)
    R = gentang(xyz, return_matrices=True)
    return np.asarray(R)


def _calibrate_1d(
    *,
    cases_data: List[Dict[str, Any]],
    dp_eval: Callable[[Dict[str, Any], float], float],
    bounds: Tuple[float, float],
    stdnum: float,
    n_grid: int,
    n_refine: int,
    param_name: str,
) -> CalibrationResult1D:
    N = np.array([c["Nf"] for c in cases_data], dtype=float)

    def objective(p: float) -> float:
        DP = np.array([dp_eval(c, float(p)) for c in cases_data], dtype=float)
        fit = fit_power_law_with_survival_std(N, DP, stdnum=stdnum)
        return float(fit["STD"]) if np.isfinite(fit["STD"]) else np.inf

    p_best, std_best = minimize_1d_grid_refine(objective, bounds, n_grid=n_grid, n_refine=n_refine)

    DP_best = np.array([dp_eval(c, float(p_best)) for c in cases_data], dtype=float)
    fit = fit_power_law_with_survival_std(N, DP_best, stdnum=stdnum)

    return CalibrationResult1D(
        param_name=str(param_name),
        p_opt=float(p_best),
        scatter_STD=float(std_best),
        fit=dict(fit),
        series={"N": N, "DP": DP_best},
    )


def _calibrate_2d(
    *,
    cases_data: List[Dict[str, Any]],
    dp_eval: Callable[[Dict[str, Any], float, float], float],
    bounds_x: Tuple[float, float],
    bounds_y: Tuple[float, float],
    stdnum: float,
    n_grid_x: int,
    n_grid_y: int,
    n_refine: int,
    param_names: Tuple[str, str],
) -> CalibrationResult2D:
    N = np.array([c["Nf"] for c in cases_data], dtype=float)

    def objective(px: float, py: float) -> float:
        DP = np.array([dp_eval(c, float(px), float(py)) for c in cases_data], dtype=float)
        fit = fit_power_law_with_survival_std(N, DP, stdnum=stdnum)
        return float(fit["STD"]) if np.isfinite(fit["STD"]) else np.inf

    x_best, y_best, std_best = minimize_2d_grid_refine(
        objective,
        bounds_x,
        bounds_y,
        n_grid_x=n_grid_x,
        n_grid_y=n_grid_y,
        n_refine=n_refine,
    )

    DP_best = np.array([dp_eval(c, float(x_best), float(y_best)) for c in cases_data], dtype=float)
    fit = fit_power_law_with_survival_std(N, DP_best, stdnum=stdnum)

    return CalibrationResult2D(
        param_names=tuple(param_names),
        p_opt=(float(x_best), float(y_best)),
        scatter_STD=float(std_best),
        fit=dict(fit),
        series={"N": N, "DP": DP_best},
    )


# -----------------------------
# ZhuEDP multiaxial selection helpers
# -----------------------------

def _harmonic_amplitudes_for_zhu(h: Mapping[str, np.ndarray]) -> Tuple[float, float, float]:
    """Return (amp_axial, amp_shear, F_NP) using the same definitions as ZhuEDP."""
    E_sin = np.asarray(h.get("E_sin"), dtype=float)
    E_cos = np.asarray(h.get("E_cos"), dtype=float)

    ex_s, ex_c = float(E_sin[1, 1]), float(E_cos[1, 1])
    g_s, g_c = 2.0 * float(E_sin[0, 1]), 2.0 * float(E_cos[0, 1])
    gy_s, gy_c = g_s / math.sqrt(3.0), g_c / math.sqrt(3.0)

    amp_ax = math.sqrt(ex_s * ex_s + ex_c * ex_c)
    amp_sh = math.sqrt(gy_s * gy_s + gy_c * gy_c)

    def _phase(sin_coeff: float, cos_coeff: float) -> float:
        if abs(sin_coeff) < 1e-16 and abs(cos_coeff) < 1e-16:
            return 0.0
        return math.atan2(cos_coeff, sin_coeff)

    phi_x = _phase(ex_s, ex_c)
    phi_y = _phase(gy_s, gy_c)
    F_NP = abs(math.sin(phi_y - phi_x))

    return float(amp_ax), float(amp_sh), float(F_NP)


# -----------------------------
# Public API
# -----------------------------

def calibrate_material_params(
    *,
    data: Dict[str, Any],
    material_name: str,
    R_list: Optional[np.ndarray] = None,
    stdnum: float = 2.0,
    node: int = 1,
    # selection tolerances
    r_target: float = -1.0,
    rtol_R_init: float = 0.02,
    rtol_uniax: float = 1e-3,
    # speed
    plane_stride: int = 1,
    # search tuning
    n_grid_1d: int = 41,
    n_refine_1d: int = 6,
    n_grid_2d: int = 21,
    n_refine_2d: int = 5,
    # bounds
    bounds_k_FS: Tuple[float, float] = (0.0, 3.0),
    bounds_k_FI: Tuple[float, float] = (0.0, 3.0),
    bounds_L_LI: Tuple[float, float] = (0.0, 5.0),
    bounds_k_MGSE: Tuple[float, float] = (0.0, 5.0),
    bounds_a_MSWT: Tuple[float, float] = (0.0, 1.0),
    bounds_alpha_w: Tuple[float, float] = (0.0, 2.0),
    bounds_xi_ZHU: Tuple[float, float] = (0.0, 1.0),
    # Zhu selection
    zhu_amp_tol: float = 1e-12,
) -> Dict[str, Any]:
    """Calibrate FS/FIN + additional parameters for one material."""

    if material_name not in data:
        raise KeyError(f"Material '{material_name}' not found. Available: {list(data.keys())}")

    if R_list is None:
        R_list = _build_R_list()

    if plane_stride is not None and int(plane_stride) > 1:
        R_use = np.asarray(R_list)[:: int(plane_stride)]
    else:
        R_use = np.asarray(R_list)

    mat = data[material_name]
    mat_prop = mat.get("material", {})
    cases = mat.get("cases", {})

    Sy = mat_prop.get("Yield_strength_MPa", None)
    if Sy is None or not np.isfinite(Sy):
        raise ValueError(f"Yield strength not found/invalid for material '{material_name}'. Needed for FS.")

    selector_kwargs = dict(
        node=node,
        r_target=r_target,
        rtol_R=rtol_R_init,
        atol_R=1e-6,
        rtol_uniax=rtol_uniax,
        atol_uniax=1e-6,
    )

    tension_ids, tension_info = _select_cases_with_fallback(
        cases,
        _uniaxial_tension_filter,
        selector_kwargs=selector_kwargs,
        name_for_warnings=f"{material_name}/uniaxial_tension",
    )
    torsion_ids, torsion_info = _select_cases_with_fallback(
        cases,
        _uniaxial_torsion_filter,
        selector_kwargs=selector_kwargs,
        name_for_warnings=f"{material_name}/uniaxial_torsion",
    )

    if len(tension_ids) == 0 or len(torsion_ids) == 0:
        warnings.warn(
            f"[{material_name}] Uniaxial dataset may be incomplete: tension={len(tension_ids)} torsion={len(torsion_ids)}.",
            RuntimeWarning,
        )

    combined_ids = tension_ids + torsion_ids

    uniax_cases: List[Dict[str, Any]] = []
    for cid in combined_ids:
        c = cases[int(cid)]
        Nf = float(c["meta"]["Nf_cycles"])
        S0, E0, S1, E1, H = exp.extract_tensor(data, material_name, int(cid), node=node, include_harmonics=True)
        uniax_cases.append({"case_id": int(cid), "Nf": Nf, "S0": S0, "E0": E0, "S1": S1, "E1": E1, "H": H})

    if len(uniax_cases) < 2:
        raise RuntimeError(f"[{material_name}] Not enough uniaxial cases to calibrate (need >= 2).")

    # ------------------
    # Instantiate methods
    # ------------------
    fs_method = FatemiSocie()
    fin_method = Findley()
    li_method = LI()
    mgse_method = MGSE_Zhu()
    mswt_method = MSWT()

    # ------------------
    # 1D CP parameter evals
    # ------------------
    def fs_eval(case_dict: Dict[str, Any], k: float) -> float:
        best = best_cp_for_case(
            S0=case_dict["S0"],
            E0=case_dict["E0"],
            S1=case_dict["S1"],
            E1=case_dict["E1"],
            R_list=R_use,
            method=fs_method,
            params={"k_FS": float(k), "Sy": float(Sy)},
            S_sin=case_dict["H"].get("S_sin"),
            S_cos=case_dict["H"].get("S_cos"),
            E_sin=case_dict["H"].get("E_sin"),
            E_cos=case_dict["H"].get("E_cos"),
        )
        return float(best.by_metric)  # plane maximizing Δγ

    def fin_eval(case_dict: Dict[str, Any], k: float) -> float:
        best = best_cp_for_case(
            S0=case_dict["S0"],
            E0=case_dict["E0"],
            S1=case_dict["S1"],
            E1=case_dict["E1"],
            R_list=R_use,
            method=fin_method,
            params={"k_FI": float(k)},
            S_sin=case_dict["H"].get("S_sin"),
            S_cos=case_dict["H"].get("S_cos"),
            E_sin=case_dict["H"].get("E_sin"),
            E_cos=case_dict["H"].get("E_cos"),
        )
        return float(best.ext)  # maximize Findley scalar

    def li_eval(case_dict: Dict[str, Any], L: float) -> float:
        best = best_cp_for_case(
            S0=case_dict["S0"],
            E0=case_dict["E0"],
            S1=case_dict["S1"],
            E1=case_dict["E1"],
            R_list=R_use,
            method=li_method,
            params={"L_LI": float(L)},
            S_sin=case_dict["H"].get("S_sin"),
            S_cos=case_dict["H"].get("S_cos"),
            E_sin=case_dict["H"].get("E_sin"),
            E_cos=case_dict["H"].get("E_cos"),
        )
        # Li defines a shear-driven plane criterion; the CP method's metric
        # is designed to select that plane -> use by_metric.
        return float(best.by_metric)

    def mgse_eval(case_dict: Dict[str, Any], k: float) -> float:
        best = best_cp_for_case(
            S0=case_dict["S0"],
            E0=case_dict["E0"],
            S1=case_dict["S1"],
            E1=case_dict["E1"],
            R_list=R_use,
            method=mgse_method,
            params={"k_MGSE": float(k)},
            S_sin=case_dict["H"].get("S_sin"),
            S_cos=case_dict["H"].get("S_cos"),
            E_sin=case_dict["H"].get("E_sin"),
            E_cos=case_dict["H"].get("E_cos"),
        )
        return float(best.ext)

    def mswt_eval(case_dict: Dict[str, Any], a: float) -> float:
        best = best_cp_for_case(
            S0=case_dict["S0"],
            E0=case_dict["E0"],
            S1=case_dict["S1"],
            E1=case_dict["E1"],
            R_list=R_use,
            method=mswt_method,
            params={"a_MSWT": float(a)},
            S_sin=case_dict["H"].get("S_sin"),
            S_cos=case_dict["H"].get("S_cos"),
            E_sin=case_dict["H"].get("E_sin"),
            E_cos=case_dict["H"].get("E_cos"),
        )
        return float(best.ext)

    # ------------------
    # Calibrate 1D params on uniaxial dataset
    # ------------------
    FS_cal = _calibrate_1d(
        cases_data=uniax_cases,
        dp_eval=fs_eval,
        bounds=bounds_k_FS,
        stdnum=stdnum,
        n_grid=n_grid_1d,
        n_refine=n_refine_1d,
        param_name="k_FS",
    )
    FI_cal = _calibrate_1d(
        cases_data=uniax_cases,
        dp_eval=fin_eval,
        bounds=bounds_k_FI,
        stdnum=stdnum,
        n_grid=n_grid_1d,
        n_refine=n_refine_1d,
        param_name="k_FI",
    )
    LI_cal = _calibrate_1d(
        cases_data=uniax_cases,
        dp_eval=li_eval,
        bounds=bounds_L_LI,
        stdnum=stdnum,
        n_grid=n_grid_1d,
        n_refine=n_refine_1d,
        param_name="L_LI",
    )
    MGSE_cal = _calibrate_1d(
        cases_data=uniax_cases,
        dp_eval=mgse_eval,
        bounds=bounds_k_MGSE,
        stdnum=stdnum,
        n_grid=n_grid_1d,
        n_refine=n_refine_1d,
        param_name="k_MGSE",
    )
    MSWT_cal = _calibrate_1d(
        cases_data=uniax_cases,
        dp_eval=mswt_eval,
        bounds=bounds_a_MSWT,
        stdnum=stdnum,
        n_grid=n_grid_1d,
        n_refine=n_refine_1d,
        param_name="a_MSWT",
    )

    # ------------------
    # ZhuEDP: 2D calibration on multiaxial cases (if available)
    # ------------------
    zhu_model = ZhuEDP()

    all_cases: List[Dict[str, Any]] = []
    multiax_cases: List[Dict[str, Any]] = []

    for cid, c in cases.items():
        try:
            Nf = float(c["meta"]["Nf_cycles"])
        except Exception:
            continue

        S0, E0, S1, E1, H = exp.extract_tensor(data, material_name, int(cid), node=node, include_harmonics=True)
        entry = {"case_id": int(cid), "Nf": Nf, "S0": S0, "E0": E0, "S1": S1, "E1": E1, "H": H}
        all_cases.append(entry)

        try:
            amp_ax, amp_sh, _fnp = _harmonic_amplitudes_for_zhu(H)
            if amp_ax > zhu_amp_tol and amp_sh > zhu_amp_tol:
                multiax_cases.append(entry)
        except Exception:
            pass

    if len(all_cases) < 2:
        warnings.warn(f"[{material_name}] ZhuEDP: not enough cases to calibrate.", RuntimeWarning)
        zhu_cal: Optional[CalibrationResult2D] = None
    else:
        zhu_cases_use = multiax_cases if len(multiax_cases) >= 2 else all_cases
        if len(multiax_cases) < 2:
            warnings.warn(
                f"[{material_name}] ZhuEDP: no (or too few) multiaxial cases with both axial+shear amplitudes; using ALL cases.",
                RuntimeWarning,
            )

        # Warn if alpha_w is likely not identifiable (F_NP ~ 0 for all selected points)
        fnp_vals = []
        for e in zhu_cases_use:
            try:
                _ax, _sh, fnp = _harmonic_amplitudes_for_zhu(e["H"])
                fnp_vals.append(float(fnp))
            except Exception:
                pass
        if fnp_vals and max(fnp_vals) < 1e-3:
            warnings.warn(
                f"[{material_name}] ZhuEDP: selected cases have F_NP≈0; alpha_w may be poorly identifiable (objective weakly depends on alpha_w).",
                RuntimeWarning,
            )

        def zhu_eval(case_dict: Dict[str, Any], alpha_w: float, xi: float) -> float:
            r = zhu_model.evaluate_case(
                S0=case_dict["S0"],
                E0=case_dict["E0"],
                S1=case_dict["S1"],
                E1=case_dict["E1"],
                params={"alpha_w": float(alpha_w), "xi_ZHU": float(xi)},
                harmonics=case_dict["H"],
                R_list=None,
            )
            return float(r.values.get("", np.nan))

        zhu_cal = _calibrate_2d(
            cases_data=zhu_cases_use,
            dp_eval=zhu_eval,
            bounds_x=bounds_alpha_w,
            bounds_y=bounds_xi_ZHU,
            stdnum=stdnum,
            n_grid_x=n_grid_2d,
            n_grid_y=n_grid_2d,
            n_refine=n_refine_2d,
            param_names=("alpha_w", "xi_ZHU"),
        )

    # ------------------
    # Pack output
    # ------------------
    out: Dict[str, Any] = {
        "material": material_name,
        "selection": {
            "target_R": float(r_target),
            "tension_case_ids": tension_ids,
            "torsion_case_ids": torsion_ids,
            "tension_info": tension_info,
            "torsion_info": torsion_info,
            "n_uniax": int(len(uniax_cases)),
            "n_all": int(len(all_cases)),
            "n_zhu_multiax": int(len(multiax_cases)),
        },
        "FS": {
            "k_FS": FS_cal.p_opt,
            "scatter_STD": FS_cal.scatter_STD,
            "fit": FS_cal.fit,
        },
        "FI": {
            "k_FI": FI_cal.p_opt,
            "scatter_STD": FI_cal.scatter_STD,
            "fit": FI_cal.fit,
        },
        "LI": {
            "L_LI": LI_cal.p_opt,
            "scatter_STD": LI_cal.scatter_STD,
            "fit": LI_cal.fit,
        },
        "MGSE_ZHU": {
            "k_MGSE": MGSE_cal.p_opt,
            "scatter_STD": MGSE_cal.scatter_STD,
            "fit": MGSE_cal.fit,
        },
        "MSWT": {
            "a_MSWT": MSWT_cal.p_opt,
            "scatter_STD": MSWT_cal.scatter_STD,
            "fit": MSWT_cal.fit,
        },
    }

    if zhu_cal is not None:
        out["ZHU_EDP"] = {
            "alpha_w": zhu_cal.p_opt[0],
            "xi_ZHU": zhu_cal.p_opt[1],
            "scatter_STD": zhu_cal.scatter_STD,
            "fit": zhu_cal.fit,
        }
    else:
        out["ZHU_EDP"] = {
            "alpha_w": np.nan,
            "xi_ZHU": np.nan,
            "scatter_STD": np.nan,
            "fit": {},
        }

    return out


def calibrate_many(
    *,
    exp_folder: str | Path = "Experimental_Data",
    materials: Optional[Iterable[str]] = None,
    stdnum: float = 2.0,
    plane_stride: int = 10,
) -> Dict[str, Any]:
    """Import data and calibrate multiple materials."""
    data = exp.import_experimental_data_grouped_by_material(exp_folder, verbose=True)
    mats = list(materials) if materials is not None else list(data.keys())

    R_list = _build_R_list()
    out: Dict[str, Any] = {}
    for m in mats:
        out[m] = calibrate_material_params(
            data=data,
            material_name=m,
            R_list=R_list,
            stdnum=stdnum,
            plane_stride=plane_stride,
        )
    return out


if __name__ == "__main__":
    data = exp.import_experimental_data_grouped_by_material("Experimental_Data", verbose=True)

    mats = list(data.keys())
    if not mats:
        raise SystemExit("No materials found in Experimental_Data")

    R_list = _build_R_list()

    for material_name in mats:
        cal = calibrate_material_params(
            data=data,
            material_name=material_name,
            R_list=R_list,
            stdnum=2.0,
            plane_stride=10,
        )

        print("\n" + "=" * 72)
        print("Material:", cal["material"])
        sel = cal["selection"]
        print(f"Uniaxial cases: tension={len(sel['tension_case_ids'])} torsion={len(sel['torsion_case_ids'])} (total={sel['n_uniax']})")
        print(f"ZhuEDP cases: all={sel['n_all']} multiax(amp>tol)={sel['n_zhu_multiax']}")

        print("\nCalibrated parameters")
        print(f"  k_FS     = {cal['FS']['k_FS']:.6g}")
        print(f"  k_FI     = {cal['FI']['k_FI']:.6g}")
        print(f"  L_LI     = {cal['LI']['L_LI']:.6g}")
        print(f"  k_MGSE   = {cal['MGSE_ZHU']['k_MGSE']:.6g}")
        print(f"  a_MSWT   = {cal['MSWT']['a_MSWT']:.6g}")
        if np.isfinite(cal['ZHU_EDP'].get('alpha_w', np.nan)):
            print(f"  alpha_w  = {cal['ZHU_EDP']['alpha_w']:.6g}")
            print(f"  xi_ZHU   = {cal['ZHU_EDP']['xi_ZHU']:.6g}")
        else:
            print("  alpha_w  = NaN (not calibrated)")
            print("  xi_ZHU   = NaN (not calibrated)")

        # Scatter summary
        print("\nScatter STD (lnN residuals)")
        print(f"  FS      STD = {cal['FS']['scatter_STD']:.6g}")
        print(f"  FI      STD = {cal['FI']['scatter_STD']:.6g}")
        print(f"  LI      STD = {cal['LI']['scatter_STD']:.6g}")
        print(f"  MGSE    STD = {cal['MGSE_ZHU']['scatter_STD']:.6g}")
        print(f"  MSWT    STD = {cal['MSWT']['scatter_STD']:.6g}")
        print(f"  ZhuEDP  STD = {cal['ZHU_EDP'].get('scatter_STD', np.nan):.6g}")
