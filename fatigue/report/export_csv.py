# fatigue/report/export_csv.py
from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Optional

import numpy as np


def _strip_variant(method_name: str) -> str:
    base = str(method_name).strip()
    if base.endswith("_ext"):
        base = base[:-4]
    return base


def _family_from_module(mod: str) -> str:
    """Extract family from a module path like 'fatigue.models.cp_methods.fin'."""
    if not isinstance(mod, str) or not mod:
        return "UNKNOWN"
    parts = mod.split(".")
    try:
        i = parts.index("models")
        if i + 1 < len(parts):
            fam = parts[i + 1]
            return fam
    except ValueError:
        pass
    return "UNKNOWN"


def infer_method_family(method_name: str) -> str:
    """Infer the method family robustly.

    Priority:
      1) Use the central registry (most reliable, no filesystem assumptions).
      2) Fallback to filesystem search under fatigue/models using fatigue.__path__.

    Returns a folder-like name, e.g. 'cp_methods', 'invariant_methods', 'energy_based_methods'.
    """
    base = _strip_variant(method_name)

    # --- 1) Registry-based inference (works even for namespace packages) ---
    try:
        from fatigue.models import registry as reg  # type: ignore

        try:
            cp = reg._available_cp_methods()  # type: ignore[attr-defined]
            if base in cp:
                return "cp_methods"
        except Exception:
            pass

        try:
            direct = reg._available_direct_models()  # type: ignore[attr-defined]
            if base in direct:
                fam = _family_from_module(direct[base].__class__.__module__)
                return fam if fam != "UNKNOWN" else "direct_models"
        except Exception:
            pass

        # last attempt: build model and inspect module
        try:
            models = reg.get_models([base])
            if models:
                m0 = models[0]
                # CP adapters: prefer underlying plane method module
                if hasattr(m0, "_m"):
                    fam = _family_from_module(getattr(m0._m, "__class__").__module__)
                    if fam != "UNKNOWN":
                        return fam
                    return "cp_methods"
                fam = _family_from_module(m0.__class__.__module__)
                if fam != "UNKNOWN":
                    return fam
        except Exception:
            pass

    except Exception:
        # registry not importable → try filesystem
        pass

    # --- 2) Filesystem fallback (best-effort) ---
    try:
        import fatigue  # type: ignore

        roots = []
        if hasattr(fatigue, "__path__"):
            roots = list(fatigue.__path__)  # type: ignore[attr-defined]
        if roots:
            fatigue_dir = Path(roots[0]).resolve()
        elif getattr(fatigue, "__file__", None):
            fatigue_dir = Path(fatigue.__file__).resolve().parent  # type: ignore[arg-type]
        else:
            return "UNKNOWN"

        models_dir = fatigue_dir / "models"
        if not models_dir.is_dir():
            return "UNKNOWN"

        stem = base.lower()
        for fam in sorted(models_dir.iterdir()):
            if not fam.is_dir() or fam.name.startswith("__"):
                continue
            for py in fam.rglob("*.py"):
                if py.stem.lower() == stem:
                    return fam.name

        # try anywhere under models
        for py in models_dir.rglob("*.py"):
            if py.stem.lower() == stem:
                rel = py.relative_to(models_dir)
                if len(rel.parts) >= 2:
                    return rel.parts[0]
                return "models"

    except Exception:
        return "UNKNOWN"

    return "UNKNOWN"


def save_fatigue_csv(
    *,
    out_dir: str,
    method_name: str,
    material_name: str,
    N_exp,
    CP,
    N_pred,
    N_design=None,
    F_design=None,
    N_pred_design=None,
    # Optional loading-type labels (same length as N_exp / N_design)
    loading_type=None,
    loading_type_design=None,
    method_family: Optional[str] = None,
    metric_mode: str = "Nf",
    # Optional DP-based metric columns (used when metric_mode == 'DP_DIFF')
    DP_fit=None,
    Error_ln_dp=None,
    # Deprecated aliases kept for backward compatibility (ignored in favor of Error_ln_dp)
    Error_dp=None,
    DP_fit_design=None,
    Error_ln_dp_design=None,
    Error_dp_design=None,
    basquin_meta: dict | None = None,   # kept for backward compatibility (ignored)
    material_meta: dict | None = None,  # kept for backward compatibility (ignored)
    extra_meta: dict | None = None,     # kept for backward compatibility (ignored)
) -> str:
    """Export a lightweight results CSV.

    Output columns:
      Material, Method_family, Method, Metric_mode, Loading_type,
      Nf_exp, CP_value, Nf_expected, DP_fit,
      Flag, Error_ln, Error_ln_dp

    Flag convention:
      0 = design points (if provided)
      1 = other points

    Error definitions:
      - Nf mode:   Error_ln = ln(Nf_exp / Nf_expected)
      - DP_DIFF:   Error_ln_dp = ln(DP_point / DP_fit), where DP_fit is the calibration best-fit evaluated at Nf_exp.

    Notes:
      * Parameters Error_dp / Error_dp_design are deprecated aliases kept only to avoid breaking older scripts.
    """
    os.makedirs(out_dir, exist_ok=True)
    filename = f"{method_name}_{material_name}.csv"
    path = os.path.join(out_dir, filename)

    family = (method_family or infer_method_family(method_name)).strip() or "UNKNOWN"

    mode = str(metric_mode or "Nf").strip() or "Nf"

    # Convert inputs
    N_exp = np.asarray(N_exp, dtype=float).ravel()
    CP = np.asarray(CP, dtype=float).ravel()
    N_pred = np.asarray(N_pred, dtype=float).ravel()

    # Loading type arrays (strings)
    if loading_type is None:
        loading_type = ["UNKNOWN"] * int(len(N_exp))
    if isinstance(loading_type, (str, bytes)):
        loading_type = [str(loading_type)] * int(len(N_exp))
    loading_type = list(loading_type)
    if len(loading_type) != len(N_exp):
        # Best-effort: truncate/pad
        if len(loading_type) > len(N_exp):
            loading_type = loading_type[: len(N_exp)]
        else:
            loading_type = loading_type + ["UNKNOWN"] * (len(N_exp) - len(loading_type))

    have_design = (N_design is not None) and (F_design is not None) and (N_pred_design is not None)
    if have_design:
        N_design = np.asarray(N_design, dtype=float).ravel()
        F_design = np.asarray(F_design, dtype=float).ravel()
        N_pred_design = np.asarray(N_pred_design, dtype=float).ravel()

        if loading_type_design is None:
            loading_type_design = ["UNKNOWN"] * int(len(N_design))
        if isinstance(loading_type_design, (str, bytes)):
            loading_type_design = [str(loading_type_design)] * int(len(N_design))
        loading_type_design = list(loading_type_design)
        if len(loading_type_design) != len(N_design):
            if len(loading_type_design) > len(N_design):
                loading_type_design = loading_type_design[: len(N_design)]
            else:
                loading_type_design = loading_type_design + ["UNKNOWN"] * (len(N_design) - len(loading_type_design))

    # DP-fit and DP-log-ratio error arrays (optional)
    DP_fit = np.asarray(DP_fit, dtype=float).ravel() if DP_fit is not None else None
    Error_ln_dp = np.asarray(Error_ln_dp, dtype=float).ravel() if Error_ln_dp is not None else None
    DP_fit_design = np.asarray(DP_fit_design, dtype=float).ravel() if DP_fit_design is not None else None
    Error_ln_dp_design = np.asarray(Error_ln_dp_design, dtype=float).ravel() if Error_ln_dp_design is not None else None

    # If not provided, but basquin_meta contains A_surv and b, compute DP_fit and Error_dp.
    try:
        As = None
        b = None
        if isinstance(basquin_meta, dict):
            As = basquin_meta.get("A_surv", basquin_meta.get("As", None))
            b = basquin_meta.get("b", None)
        As = float(As) if As is not None else None
        b = float(b) if b is not None else None
    except Exception:
        As, b = None, None

    def _dp_fit(n: np.ndarray) -> np.ndarray:
        out = np.full_like(n, np.nan, dtype=float)
        if As is None or b is None:
            return out
        m = np.isfinite(n) & (n > 0) & np.isfinite(As) & np.isfinite(b)
        out[m] = As * (n[m] ** b)
        return out

    if DP_fit is None:
        DP_fit = _dp_fit(N_exp)
    def _safe_log_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
        num = np.asarray(num, dtype=float)
        den = np.asarray(den, dtype=float)
        out = np.full_like(num, np.nan, dtype=float)
        m = np.isfinite(num) & np.isfinite(den) & (num > 0.0) & (den > 0.0)
        out[m] = np.log(num[m] / den[m])
        return out

    if Error_ln_dp is None and DP_fit is not None:
        Error_ln_dp = _safe_log_ratio(np.asarray(CP, dtype=float), np.asarray(DP_fit, dtype=float))

    if have_design:
        if DP_fit_design is None:
            DP_fit_design = _dp_fit(N_design)
        if Error_ln_dp_design is None and DP_fit_design is not None:
            Error_ln_dp_design = _safe_log_ratio(np.asarray(F_design, dtype=float), np.asarray(DP_fit_design, dtype=float))

    def _safe_err_ln(nexp: np.ndarray, npred: np.ndarray) -> np.ndarray:
        nexp = np.asarray(nexp, dtype=float)
        npred = np.asarray(npred, dtype=float)
        out = np.full_like(nexp, np.nan, dtype=float)
        m = np.isfinite(nexp) & np.isfinite(npred) & (nexp > 0.0) & (npred > 0.0)
        out[m] = np.log(nexp[m] / npred[m])
        return out

    err_other = _safe_err_ln(N_exp, N_pred) if mode.upper() != "DP_DIFF" else np.full_like(N_exp, np.nan, dtype=float)

    if have_design:
        err_design = (
            _safe_err_ln(N_design, N_pred_design)
            if mode.upper() != "DP_DIFF"
            else np.full_like(N_design, np.nan, dtype=float)
        )
    else:
        err_design = None

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "Material",
                "Method_family",
                "Method",
                "Metric_mode",
                "Loading_type",
                "Nf_exp",
                "CP_value",
                "Nf_expected",
                "DP_fit",
                "Flag",
                "Error_ln",
                "Error_ln_dp",
            ]
        )

        if have_design:
            for i, (ne, cp, npred, err) in enumerate(zip(N_design, F_design, N_pred_design, err_design)):
                dp_fit_i = float(DP_fit_design[i]) if DP_fit_design is not None and i < len(DP_fit_design) else float("nan")
                edp_i = float(Error_ln_dp_design[i]) if Error_ln_dp_design is not None and i < len(Error_ln_dp_design) else float("nan")
                w.writerow(
                    [
                        material_name,
                        family,
                        method_name,
                        mode,
                        str(loading_type_design[i]) if loading_type_design is not None and i < len(loading_type_design) else "UNKNOWN",
                        float(ne),
                        float(cp),
                        float(npred),
                        dp_fit_i,
                        0,
                        float(err),
                        edp_i,
                    ]
                )

        for i, (ne, cp, npred, err) in enumerate(zip(N_exp, CP, N_pred, err_other)):
            dp_fit_i = float(DP_fit[i]) if DP_fit is not None and i < len(DP_fit) else float("nan")
            edp_i = float(Error_ln_dp[i]) if Error_ln_dp is not None and i < len(Error_ln_dp) else float("nan")
            w.writerow(
                [
                    material_name,
                    family,
                    method_name,
                    mode,
                    str(loading_type[i]) if i < len(loading_type) else "UNKNOWN",
                    float(ne),
                    float(cp),
                    float(npred),
                    dp_fit_i,
                    1,
                    float(err),
                    edp_i,
                ]
            )

    return path
