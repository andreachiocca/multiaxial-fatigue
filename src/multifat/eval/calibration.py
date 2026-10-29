# multifat/eval/calibration.py
from __future__ import annotations

from typing import Any, Mapping
import numpy as np

from multifat.eval.fitting import fit_power_law_with_survival_std
from multifat.models.base import FatigueModel
from multifat.models.comparison import capabilities_for_model


def collect_uniaxial_series_and_fit(
    *,
    data,
    exp,  # the experimental_data module in use
    material_name: str,
    uniax_cases: dict,  # result of exp.extract_uniaxial_tension(data)
    R_list: np.ndarray | None,
    models: list[FatigueModel],
    params: Mapping[str, float],
    stdnum: float = 2.0,
) -> dict[str, Any]:
    """
    Shared calibration routine that works for CP and non-CP models alike.

    A model can expose several variants via CaseResult.values:
      - ""     : the nominal value (saved under "CP")
      - "_ext" : the extreme/maximum value (saved under "CP_ext"), when present

    The returned mapping keeps the original backward-compatible layout:
      {
        model.name: {
          "series": {
              "N": array,
              "CP": array,        # comes from variant ""
              "CP_ext": array     # comes from variant "_ext" (empty if absent)
          },
          "fit": {...} ,         # fitted on CP (variant "")
          "fit_ext": {...}       # fitted on CP_ext (variant "_ext"), if any
        }, ...
      }

    Notes:
    - Non-CP models (VM/Tresca/energy, ...) usually only provide variant "",
      so CP_ext stays empty and fit_ext is None.
    - CP models wrapped in CPModelAdapter normally provide both "" and "_ext".
    """
    results_by_model: dict[str, Any] = {}

    if material_name not in uniax_cases:
        raise KeyError(
            f"material_name='{material_name}' is absent from uniax_cases keys: {list(uniax_cases.keys())}"
        )

    load_cases = uniax_cases[material_name]["points"]

    for model in models:
        # accumulators for the nominal and extreme series
        n_std_vals: list[float] = []
        cp_std_vals: list[float] = []

        n_ext_vals: list[float] = []
        cp_ext_vals: list[float] = []
        std_case_ids: list[int] = []
        ext_case_ids: list[int] = []

        # verify up front that every needed material parameter is supplied
        absent_params = [k for k in model.required_params() if params.get(k, None) is None]
        if absent_params:
            raise ValueError(
                f"Model '{model.name}' is missing required material parameters: {absent_params}. "
                f"Keys currently available: {list(params.keys())}"
            )

        for p in load_cases:
            cid = int(p["case_id"])
            Nf = float(p["Nf_cycles"])

            if not (np.isfinite(Nf) and Nf > 0):
                continue

            s_mean, e_mean, s_amp, e_amp, harmonics = exp.extract_tensor(
                data, material_name, cid, node=1, include_harmonics=True
            )

            # Evaluate at model level (CP models may use R_list; direct models disregard it)
            evaluation = model.evaluate_case(
                S0=s_mean,
                E0=e_mean,
                S1=s_amp,
                E1=e_amp,
                params=params,
                R_list=R_list,
                harmonics=harmonics,
            )

            # nominal variant ("")
            if "" in evaluation.values:
                value_std = float(evaluation.values[""])
                if np.isfinite(value_std) and value_std > 0:
                    n_std_vals.append(Nf)
                    cp_std_vals.append(value_std)
                    std_case_ids.append(cid)

            # extreme variant ("_ext"), only if offered
            if "_ext" in evaluation.values:
                value_ext = float(evaluation.values["_ext"])
                if np.isfinite(value_ext) and value_ext > 0:
                    n_ext_vals.append(Nf)
                    cp_ext_vals.append(value_ext)
                    ext_case_ids.append(cid)

        # Turn the accumulators into arrays
        n_std_arr = np.asarray(n_std_vals, dtype=float)
        cp_std_arr = np.asarray(cp_std_vals, dtype=float)

        n_ext_arr = np.asarray(n_ext_vals, dtype=float)
        cp_ext_arr = np.asarray(cp_ext_vals, dtype=float)

        # Native life comes from the constitutive law; never fit its diagnostic index.
        native_prediction = capabilities_for_model(model).life_prediction == "native"
        # Perform the fits
        std_fit = (
            fit_power_law_with_survival_std(n_std_arr, cp_std_arr, stdnum=stdnum)
            if len(n_std_arr) > 0 and not native_prediction
            else None
        )
        ext_fit = (
            fit_power_law_with_survival_std(n_ext_arr, cp_ext_arr, stdnum=stdnum)
            if len(n_ext_arr) > 0 and not native_prediction
            else None
        )

        # Keep the output layout backward compatible
        results_by_model[model.name] = {
            "series": {
                "N": n_std_arr,  # N paired with the fitting damage (by_metric or direct)
                "CP": cp_std_arr,  # damage of the fitting curve
                "CP_ext": cp_ext_arr,  # damage of the extreme variant (possibly empty)
            },
            "fit": std_fit,  # None when there are no points
            "fit_ext": ext_fit,  # None when the model has no "_ext" variant
            # Case IDs behind the series points (upstream 0c9ea48: design provenance)
            "design_case_ids": std_case_ids,
            "design_case_ids_ext": ext_case_ids,
        }

    return results_by_model
