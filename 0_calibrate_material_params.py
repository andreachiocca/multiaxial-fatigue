# calibrate_cp_params.py
"""Calibration utility for critical-plane / model parameters.

Originally used to calibrate:
- k_FS for Fatemi–Socie (FS)
- k_FI for Findley (FIN)

Extended (per your request) to also calibrate additional material parameters:
- L_LI     (Li 2021 CP method)
- k_MGSE   (MGSE_ZHU CP method)
- a_MSWT   (MSWT CP method)
ZHU_EDP calibration was removed by the numerical audit; see docs/method-audit.md.

Repository conventions (v6)
--------------------------
- Experimental data importer: fatigue.imp.experimental_data
- Plane sampling: fatigue.planes.sampling (dcmesh_xyz, gentang)
- Critical-plane evaluation: fatigue.eval.evaluator.best_cp_for_case
- CP methods: fatigue.models.cp_methods

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
    if n_grid < 2 or n_refine < 1:
        raise ValueError("Calibration requires n_grid >= 2 and n_refine >= 1")
    if not (hi > lo):
        raise ValueError("bounds must satisfy hi > lo")

    x_best = np.nan
    f_best = np.inf

    for _ in range(int(n_refine)):
        xs = np.linspace(lo, hi, int(n_grid))
        fs = np.array([objective(float(x)) for x in xs], dtype=float)

        if not np.any(np.isfinite(fs)):
            raise ValueError("No finite calibration objective: check design data and master-curve slope")
        j = int(np.argmin(np.where(np.isfinite(fs), fs, np.inf)))
        x_best = float(xs[j])
        f_best = float(fs[j])

        step = (hi - lo) / max(1, (len(xs) - 1))
        lo = max(bounds[0], x_best - 2.0 * step)
        hi = min(bounds[1], x_best + 2.0 * step)

        if hi - lo < 1e-10:
            break

    return x_best, f_best



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
    # bounds
    bounds_k_FS: Tuple[float, float] = (0.0, 3.0),
    bounds_k_FI: Tuple[float, float] = (0.0, 3.0),
    bounds_L_LI: Tuple[float, float] = (0.0, 5.0),
    bounds_k_MGSE: Tuple[float, float] = (0.0, 5.0),
    bounds_a_MSWT: Tuple[float, float] = (0.0, 1.0),
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

        print("\nCalibrated parameters")
        print(f"  k_FS     = {cal['FS']['k_FS']:.6g}")
        print(f"  k_FI     = {cal['FI']['k_FI']:.6g}")
        print(f"  L_LI     = {cal['LI']['L_LI']:.6g}")
        print(f"  k_MGSE   = {cal['MGSE_ZHU']['k_MGSE']:.6g}")
        print(f"  a_MSWT   = {cal['MSWT']['a_MSWT']:.6g}")
        # Scatter summary
        print("\nScatter STD (lnN residuals)")
        print(f"  FS      STD = {cal['FS']['scatter_STD']:.6g}")
        print(f"  FI      STD = {cal['FI']['scatter_STD']:.6g}")
        print(f"  LI      STD = {cal['LI']['scatter_STD']:.6g}")
        print(f"  MGSE    STD = {cal['MGSE_ZHU']['scatter_STD']:.6g}")
        print(f"  MSWT    STD = {cal['MSWT']['scatter_STD']:.6g}")
