# run_material.py
from __future__ import annotations

import os
import numpy as np

import fatigue.imp.experimental_data as exp
import fatigue.eval.fat_cycles as fatc
from fatigue.planes.sampling import build_R_list
from fatigue.eval.calibration import collect_uniaxial_series_and_fit
from fatigue.report.export_csv import save_fatigue_csv
from fatigue.models.registry import get_models


# ----------------------------
# Configuration (edit these)
# ----------------------------
DATA_FOLDER = "Experimental_Data"
RESULTS_DIR = "Results"

MATERIAL_NAME = "AISI316L"   # e.g. "42CrMo4_QT", "AISI316L", "Al7075_T6", "Ti6Al4V_andrea_(AXIAL)"
MODEL_NAME = "FS_ext"              # e.g. "FS", "FIN", "SWT", "BP", "CAIM"


DCMESH_N1 = 20
DCMESH_N2 = 2 * DCMESH_N1
STDNUM = 0
MAKE_PLOTS = True

def _variant_from_model_name(name: str) -> tuple[str, str]:
    """
    Returns (base_model_name, variant_key).

    Convention:
      - "FIN"     -> ("FIN", "")      standard variant
      - "FIN_ext" -> ("FIN", "_ext")  extreme variant
      - "VM"      -> ("VM", "")       standard only
    """
    if name.endswith("_ext"):
        return name[:-4], "_ext"
    return name, ""


def classify_loading_type(*, H: dict, meta: dict, tol: float = 1e-6, phase_tol_deg: float = 1e-6) -> str:
    """Classify loading type from the imported stress harmonics + phase metadata.

    Categories (as requested):
      - "Tensile"
      - "Torsion"
      - "Multiaxial proportional"
      - "Multiaxial non-proportional"

    Logic:
      * Active components are determined from mean OR harmonic amplitude.
      * Uniaxial normal (axial or radial only) -> Tensile.
      * Uniaxial shear only -> Torsion.
      * Multiaxial -> proportional if all relevant phase shifts are 0°, otherwise non-proportional.
    """
    try:
        S_mean = np.asarray(H.get("S_mean"), dtype=float)
        S_sin = np.asarray(H.get("S_sin"), dtype=float)
        S_cos = np.asarray(H.get("S_cos"), dtype=float)
        amp = np.sqrt(S_sin**2 + S_cos**2)

        # Components used by your importer:
        #  radial -> sigma_xx, axial -> sigma_yy, torsion -> tau_xy
        active_rad = (abs(S_mean[0, 0]) > tol) or (abs(amp[0, 0]) > tol)
        active_ax = (abs(S_mean[1, 1]) > tol) or (abs(amp[1, 1]) > tol)
        active_tor = (abs(S_mean[0, 1]) > tol) or (abs(amp[0, 1]) > tol)

        n_active = int(active_rad) + int(active_ax) + int(active_tor)
        if n_active <= 0:
            return "UNKNOWN"

        # Uniaxial cases
        if n_active == 1:
            if active_tor:
                return "Torsion"
            # axial-only OR radial-only normal stress
            return "Tensile"

        # Multiaxial
        phi_r = float(meta.get("phase_rad_deg", 0.0) or 0.0) if active_rad else 0.0
        phi_t = float(meta.get("phase_tor_deg", 0.0) or 0.0) if active_tor else 0.0
        nonprop = (abs(phi_r) > phase_tol_deg) or (abs(phi_t) > phase_tol_deg)
        return "Multiaxial non-proportional" if nonprop else "Multiaxial proportional"
    except Exception:
        return "UNKNOWN"


def main():
    # ----------------------------
    # 1) Plane sampling
    # ----------------------------
    # Only CP models actually require R_list; direct models ignore it.
    R_list = build_R_list(n1=DCMESH_N1, n2=DCMESH_N2)
    print("R_list:", np.asarray(R_list).shape)

    # ----------------------------
    # 2) Import experimental data (same as MAIN)
    # ----------------------------
    data = exp.import_experimental_data_grouped_by_material(DATA_FOLDER, verbose=True)

    # ----------------------------
    # 3) Material properties (same keys as MAIN)
    # ----------------------------
    mat_prop = exp.get_material_properties(data, material=MATERIAL_NAME)

    Sy = mat_prop.get("Yield_strength_MPa", None)
    Su = mat_prop.get("Tensile_strength_MPa", None)
    E = mat_prop.get("E_MPa", None)
    nu = mat_prop.get("nu", None)
    Sigm1 = mat_prop.get("Tensile_fatigue_limit_MPa", None)
    Taum1 = mat_prop.get("Torsional_fatigue_limit_MPa", None)
    Sig0 = mat_prop.get("Tensile_fatigue_R0_MPa", None)
    Tau0 = mat_prop.get("Torsional_fatigue_R0_MPa", None)
    k_FS = mat_prop.get("k_FS", None)
    k_FI = mat_prop.get("k_FI", None)  # your dataset uses "k_FI"

    # Optional advanced parameters (may or may not exist in Summary sheet)
    k_MGSE = mat_prop.get("k_MGSE", None)
    a_MSWT = mat_prop.get("a_MSWT", None)
    alpha_w = mat_prop.get("alpha_w", None)
    xi_ZHU = mat_prop.get("xi_ZHU", None)

    # Optional Li-2021 weight index
    L_LI = mat_prop.get("L_LI", None)

    # Optional strain-life parameters (if you add them in Summary)
    sigma_f = mat_prop.get("sigma_f", None)
    epsilon_f = mat_prop.get("epsilon_f", None)
    b_fat = mat_prop.get("b_fat", None)
    c_fat = mat_prop.get("c_fat", None)
    tau_f = mat_prop.get("tau_f", None)
    gamma_f = mat_prop.get("gamma_f", None)
    b0 = mat_prop.get("b0", None)
    c0 = mat_prop.get("c0", None)

    params = {
        "Sy": float(Sy) if Sy is not None else None,
        "Su": float(Su) if Su is not None else None,
        "E": float(E) if E is not None else None,
        "nu": float(nu) if nu is not None else None,
        "Sigm1": float(Sigm1) if Sigm1 is not None else None,
        "Taum1": float(Taum1) if Taum1 is not None else None,
        "Sig0": float(Sig0) if Sig0 is not None else None,
        "Tau0": float(Tau0) if Tau0 is not None else None,
        "k_FS": float(k_FS) if k_FS is not None else None,
        "k_FI": float(k_FI) if k_FI is not None else None,

        # optional added constants
        "k_MGSE": float(k_MGSE) if k_MGSE is not None else None,
        "a_MSWT": float(a_MSWT) if a_MSWT is not None else None,
        "alpha_w": float(alpha_w) if alpha_w is not None else None,
        "xi_ZHU": float(xi_ZHU) if xi_ZHU is not None else None,

        "L_LI": float(L_LI) if L_LI is not None else None,

        # optional strain-life parameters
        "sigma_f": float(sigma_f) if sigma_f is not None else None,
        "epsilon_f": float(epsilon_f) if epsilon_f is not None else None,
        "b_fat": float(b_fat) if b_fat is not None else None,
        "c_fat": float(c_fat) if c_fat is not None else None,
        "tau_f": float(tau_f) if tau_f is not None else None,
        "gamma_f": float(gamma_f) if gamma_f is not None else None,
        "b0": float(b0) if b0 is not None else None,
        "c0": float(c0) if c0 is not None else None,
        # add future model params here, e.g. "E": ..., "nu": ..., "k_MY": ...
    }

    print(mat_prop)

    # ----------------------------
    # 4) Design subset (same as MAIN)
    # ----------------------------
    uniax_Rm1 = exp.extract_uniaxial_tension(data)

    # ----------------------------
    # 5) Selected model + calibration on uniaxial points
    # ----------------------------
    # Only the requested model contributes to this run's CSV. For a CP "_ext"
    # variant, calibrate its base model (which returns both variants).
    base_name, variant_key = _variant_from_model_name(MODEL_NAME)
    model = get_models([base_name])[0]

    if base_name in {"BP", "CAIM"} and (Sig0 is None or Tau0 is None):
        from fatigue.models.integral_methods.bohme_papuga import _strengths
        _, _, estimated_s0, estimated_t0 = _strengths(params, base_name)
        print(f"{base_name} R=0 strengths (MPa): Sig0={estimated_s0:.3f} "
              f"({'estimated' if Sig0 is None else 'workbook'}), "
              f"Tau0={estimated_t0:.3f} "
              f"({'estimated' if Tau0 is None else 'workbook'})")

    # Check before the potentially expensive calibration.
    missing = [k for k in model.required_params() if params.get(k, None) is None]
    if missing:
        raise ValueError(
            f"Missing required material params for {model.name}: {missing}. "
            f"Available: {params}"
        )

    cal = collect_uniaxial_series_and_fit(
        data=data,
        exp=exp,
        material_name=MATERIAL_NAME,
        uniax_cases=uniax_Rm1,
        R_list=R_list,
        models=[model],
        params=params,
        stdnum=STDNUM,
    )

    # Pull calibration for this model
    if base_name not in cal:
        raise KeyError(f"No calibration entry for model '{base_name}'. Available: {list(cal.keys())}")

    series = cal[base_name]["series"]
    
    # Your calibration returns "fit" and "fit_ext" directly (not a dict-of-variants)
    fit = cal[base_name]["fit"] if variant_key == "" else cal[base_name]["fit_ext"]

    if fit is None:
        # This can happen if:
        # - you asked for "_ext" on a non-CP model, OR
        # - there were no valid uniaxial design points for this model (non-positive/NaN values)
        raise ValueError(
            f"No available fit for model '{base_name}' with variant '{variant_key}'. "
            f"fit exists: {cal[base_name]['fit'] is not None}, "
            f"fit_ext exists: {cal[base_name]['fit_ext'] is not None}."
        )

    # Match MAIN variables
    A = fit.get("A", None)
    b = fit["b"]
    As = fit["A_surv"]
    if not (np.isfinite(As) and As > 0 and np.isfinite(b)):
        raise ValueError("Invalid master-curve fit: at least two distinct positive experimental lives are required")

    # ----------------------------
    # 5b) Decide error metric (Nf prediction vs DP difference)
    # ----------------------------
    # Power law: DP = A * N^b with b = -1/k
    # If k is too large (almost horizontal curve), predicted N becomes meaningless.
    # Switch to DP log-ratio metric:
    #   Error_ln_dp = ln(DP_point / DP_fit(Nf_exp))
    K_SLOPE_MAX = 15.0
    metric_mode = "Nf"
    k_slope = np.inf
    try:
        b_float = float(b)
        if (not np.isfinite(b_float)) or (b_float >= 0.0) or (b_float == 0.0):
            metric_mode = "DP_DIFF"
        else:
            k_slope = -1.0 / b_float
            if np.isfinite(k_slope) and k_slope > K_SLOPE_MAX:
                metric_mode = "DP_DIFF"
    except Exception:
        metric_mode = "DP_DIFF"

    print(f"Metric mode for {MODEL_NAME} / {MATERIAL_NAME}: {metric_mode} (b={b}, k={k_slope})")

    # Design points arrays
    N_d = series["N"]
    F_d = series["CP"] if variant_key == "" else series["CP_ext"]

    # Loading type for design points (uniaxial tension R=-1 by construction)
    loading_type_design = ["Tensile"] * int(len(np.atleast_1d(N_d)))

    if F_d is None or len(np.atleast_1d(F_d)) == 0:
        raise ValueError(
            f"Design series is empty for model '{base_name}' variant '{variant_key}'. "
            "Check that your uniaxial design subset exists and the model produces positive finite values."
        )

    # Predicted cycles on design points only make sense in Nf mode
    if metric_mode == "Nf":
        Npred_d = fatc.fatigue_cycles(F_d, As, b)
    else:
        Npred_d = np.full_like(np.asarray(N_d, dtype=float), np.nan, dtype=float)

    # DP_fit on design points (always computable from calibration curve)
    N_d_arr = np.asarray(N_d, dtype=float)
    DP_fit_d = np.full_like(N_d_arr, np.nan, dtype=float)
    md = np.isfinite(N_d_arr) & (N_d_arr > 0) & np.isfinite(As) & np.isfinite(b)
    DP_fit_d[md] = float(As) * (N_d_arr[md] ** float(b))
    # DP log-ratio error on design points
    F_d_arr = np.asarray(F_d, dtype=float)
    Error_ln_dp_d = np.full_like(F_d_arr, np.nan, dtype=float)
    md2 = np.isfinite(F_d_arr) & np.isfinite(DP_fit_d) & (F_d_arr > 0) & (DP_fit_d > 0)
    Error_ln_dp_d[md2] = np.log(F_d_arr[md2] / DP_fit_d[md2])

    # ----------------------------
    # 6) Exclude design subset and evaluate remaining cases (same as MAIN)
    # ----------------------------
    data_no_uniax = exp.exclude_subdataset(data, uniax_Rm1)

    values: list[float] = []
    Nexp: list[float] = []
    Npred: list[float] = []
    DP_fit: list[float] = []
    Error_ln_dp: list[float] = []
    loading_type: list[str] = []

    for case_id in data_no_uniax[MATERIAL_NAME]["cases"].keys():
        S0, E0, S1, E1, H = exp.extract_tensor(
            data_no_uniax, MATERIAL_NAME, int(case_id), node=1, include_harmonics=True
        )

        res = model.evaluate_case(
            S0=S0, E0=E0, S1=S1, E1=E1,
            params=params,
            R_list=R_list,  # CP uses it, non-CP ignores it
            harmonics=H,
        )

        if variant_key not in res.values:
            raise ValueError(
                f"Case result from model '{model.name}' missing variant '{variant_key}'. "
                f"Available: {list(res.values.keys())}"
            )

        F = float(res.values[variant_key])

        # Experimental cycles for this case
        Nf = float(data_no_uniax[MATERIAL_NAME]["cases"][int(case_id)]["meta"]["Nf_cycles"]) 

        # Loading type classification
        meta = data_no_uniax[MATERIAL_NAME]["cases"][int(case_id)].get("meta", {})
        loading_type.append(classify_loading_type(H=H, meta=meta))

        # Calibration-fit DP at experimental Nf
        dp_fit_i = np.nan
        if np.isfinite(Nf) and Nf > 0 and np.isfinite(As) and np.isfinite(b):
            dp_fit_i = float(As) * (Nf ** float(b))

        # Metric branch
        if metric_mode == "Nf" and np.isfinite(F) and F > 0:
            N_pred = fatc.fatigue_cycles(F, As, b)
        else:
            N_pred = np.nan

        # DP log-ratio error: ln(DP_point / DP_fit)
        if np.isfinite(dp_fit_i) and dp_fit_i > 0.0 and np.isfinite(F) and F > 0.0:
            err_ln_dp_i = float(np.log(float(F) / float(dp_fit_i)))
        else:
            err_ln_dp_i = np.nan

        values.append(float(F))
        Nexp.append(float(Nf))
        Npred.append(float(N_pred) if np.isfinite(N_pred) else float("nan"))
        DP_fit.append(float(dp_fit_i) if np.isfinite(dp_fit_i) else float("nan"))
        Error_ln_dp.append(float(err_ln_dp_i) if np.isfinite(err_ln_dp_i) else float("nan"))

    # ----------------------------
    # 7) Plot (same call/signature as MAIN)
    # ----------------------------
    # Plot only in Nf mode (parity plot requires meaningful N_pred)
    if metric_mode == "Nf" and MAKE_PLOTS:
        fatc.fatigue_plots(
            values, Nexp, Npred,
            A=As, B=b,
            F_design=F_d,
            N_design=N_d,
        )

    # ----------------------------
    # 8) Save CSV (same design points handling + metadata style as MAIN)
    # ----------------------------
    os.makedirs(RESULTS_DIR, exist_ok=True)

    basquin_meta = {
        "A": A,
        "b": b,
        "A_surv": As,
        "STDNUM": STDNUM,
    }
    material_meta = {
        "Sy_MPa": params.get("Sy"),
        "k_FS": params.get("k_FS"),
        "k_FI": params.get("k_FI"),
    }
    extra_meta = {
        "model_name": MODEL_NAME,
        "base_model": base_name,
        "variant": variant_key,
        "n_planes": len(R_list),
        "dcmesh_n1": DCMESH_N1,
        "dcmesh_n2": DCMESH_N2,
    }

    csv_path = save_fatigue_csv(
        out_dir=RESULTS_DIR,
        method_name=MODEL_NAME,        # keep filename compatibility
        material_name=MATERIAL_NAME,
        N_exp=Nexp,
        CP=values,
        N_pred=Npred,
        N_design=N_d,
        F_design=F_d,
        N_pred_design=Npred_d,
        loading_type=loading_type,
        loading_type_design=loading_type_design,
        metric_mode=metric_mode,
        DP_fit=DP_fit,
        Error_ln_dp=Error_ln_dp,
        DP_fit_design=DP_fit_d,
        Error_ln_dp_design=Error_ln_dp_d,
        basquin_meta=basquin_meta,
        material_meta=material_meta,
        extra_meta=extra_meta,
    )
    print("Saved:", csv_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default=MODEL_NAME)
    parser.add_argument("--material", default=MATERIAL_NAME)
    parser.add_argument("--results-dir", default=RESULTS_DIR)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    MODEL_NAME, MATERIAL_NAME, RESULTS_DIR = args.model, args.material, args.results_dir
    MAKE_PLOTS = not args.no_plots
    main()
