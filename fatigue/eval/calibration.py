# fatigue/eval/calibration.py
from __future__ import annotations

from typing import Any, Mapping
import numpy as np

from fatigue.eval.fitting import fit_power_law_with_survival_std
from fatigue.models.base import FatigueModel


def collect_uniaxial_series_and_fit(
    *,
    data,
    exp,                    # your experimental_data module
    material_name: str,
    uniax_cases: dict,      # output of exp.extract_uniaxial_tension(data)
    R_list: np.ndarray | None,
    models: list[FatigueModel],
    params: Mapping[str, float],
    stdnum: float = 2.0,
) -> dict[str, Any]:
    """
    Generic calibration for BOTH CP and non-CP models.

    Each model returns one or more variants through CaseResult.values:
      - ""     : standard value (stored as "CP")
      - "_ext" : extreme/max value (stored as "CP_ext") if available

    Returns (backward-compatible structure):
      {
        model.name: {
          "series": {
              "N": array,
              "CP": array,        # variant ""
              "CP_ext": array     # variant "_ext" (empty if not provided)
          },
          "fit": {...} ,         # fit on CP (variant "")
          "fit_ext": {...}       # fit on CP_ext (variant "_ext") if present
        }, ...
      }

    Notes:
    - For non-CP models (e.g. VM/Tresca/energy), typically only variant "" exists.
      In that case, CP_ext will be empty and fit_ext will be None.
    - For CP models wrapped by CPModelAdapter, both "" and "_ext" are usually available.
    """
    out: dict[str, Any] = {}

    if material_name not in uniax_cases:
        raise KeyError(f"material_name='{material_name}' not found in uniax_cases keys: {list(uniax_cases.keys())}")

    points = uniax_cases[material_name]["points"]

    for model in models:
        # lists for standard and ext variants
        N_std_list: list[float] = []
        CP_std_list: list[float] = []

        N_ext_list: list[float] = []
        CP_ext_list: list[float] = []

        # early check for required params
        missing = [k for k in model.required_params() if params.get(k, None) is None]
        if missing:
            raise ValueError(
                f"Missing required material params for model '{model.name}': {missing}. "
                f"Available keys: {list(params.keys())}"
            )

        for p in points:
            case_id = int(p["case_id"])
            Nf = float(p["Nf_cycles"])

            if not (np.isfinite(Nf) and Nf > 0):
                continue

            S0, E0, S1, E1, H = exp.extract_tensor(
                data, material_name, case_id, node=1, include_harmonics=True
            )

            # Model-level evaluation (CP models may need R_list; direct models ignore it)
            res = model.evaluate_case(
                S0=S0, E0=E0, S1=S1, E1=E1,
                params=params,
                R_list=R_list,
                harmonics=H,
            )

            # Standard variant ("")
            if "" in res.values:
                v = float(res.values[""])
                if np.isfinite(v) and v > 0:
                    N_std_list.append(Nf)
                    CP_std_list.append(v)

            # Extreme variant ("_ext"), optional
            if "_ext" in res.values:
                vext = float(res.values["_ext"])
                if np.isfinite(vext) and vext > 0:
                    N_ext_list.append(Nf)
                    CP_ext_list.append(vext)

        # Convert to arrays
        N_std = np.asarray(N_std_list, dtype=float)
        CP_std = np.asarray(CP_std_list, dtype=float)

        N_ext = np.asarray(N_ext_list, dtype=float)
        CP_ext = np.asarray(CP_ext_list, dtype=float)

        # Fit
        fit_std = fit_power_law_with_survival_std(N_std, CP_std, stdnum=stdnum) if len(N_std) > 0 else None
        fit_ext = fit_power_law_with_survival_std(N_ext, CP_ext, stdnum=stdnum) if len(N_ext) > 0 else None

        # Backward-compatible output structure
        out[model.name] = {
            "series": {
                "N": N_std,        # keep legacy meaning: N associated to CP (by_metric or direct)
                "CP": CP_std,      # legacy "CP"
                "CP_ext": CP_ext,  # legacy "CP_ext" (may be empty)
            },
            "fit": fit_std,        # may be None if no points
            "fit_ext": fit_ext,    # may be None if model doesn't provide "_ext"
        }

    return out
