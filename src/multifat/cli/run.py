# run_material.py
from __future__ import annotations

import hashlib
import json
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

import multifat
import multifat.imp.experimental_data as exp
import multifat.eval.fat_cycles as fatc
from multifat.imp.schema import (
    ValidationReport,
    WorkbookError,
    material_id,
    workbook_files,
)
from multifat.planes.sampling import build_R_list
from multifat.eval.calibration import collect_uniaxial_series_and_fit
from multifat.report.export_csv import save_fatigue_csv
from multifat.report.subsets import observation_index
from multifat.models.errors import MissingParamsError, NativeLifeError
from multifat.eval.findley_calibration import (
    FindleyCalibrationUnavailable,
    amplitude_identity,
    findley_method,
    require_stress_harmonics,
)
from multifat.models.references import IMPLEMENTATION_REVISION
from multifat.models.registry import get_models
from multifat.models.comparison import capabilities_for_model


# ----------------------------
# Settings (adjust as needed)
# ----------------------------
DATA_FOLDER = "Experimental_Data"
RESULTS_DIR = "Results/current"

MATERIAL_NAME = "AISI316L"  # e.g. "42CrMo4_QT", "AISI316L", "Al7075_T6", "Ti6Al4V_andrea_(AXIAL)"
MODEL_NAME = "FS"  # defaults to the author/maintainer plane-selection rule
# Findley only (FIN, FIN_ext): method_default is the production MCE path; mcc_harmonic,
# mce_harmonic and rss run the explicit harmonic assessment, calibrated with the same option.
SHEAR_AMPLITUDE_ALGORITHM = "method_default"


DCMESH_N1 = 20
DCMESH_N2 = 2 * DCMESH_N1
STDNUM = 0
MAKE_PLOTS = True
K_SLOPE_MAX = 15.0

# Workbook material keys -> model parameter names.
PARAM_KEYS = {
    "Sy": "Yield_strength_MPa",
    "Su": "Tensile_strength_MPa",
    "E": "E_MPa",
    "nu": "nu",
    "Sigm1": "Tensile_fatigue_limit_MPa",
    "Taum1": "Torsional_fatigue_limit_MPa",
    "Sig0": "Tensile_fatigue_R0_MPa",
    "Tau0": "Torsional_fatigue_R0_MPa",
    "k_FS": "k_FS",
    "k_FI": "k_FI",
    # further optional constants
    "k_MGSE": "k_MGSE",
    "a_MSWT": "a_MSWT",
    "alpha_w": "alpha_w",
    "xi_ZHU": "xi_ZHU",
    "L_LI": "L_LI",
    # optional strain-life coefficients
    "sigma_f": "sigma_f",
    "epsilon_f": "epsilon_f",
    "b_fat": "b_fat",
    "c_fat": "c_fat",
    "tau_f": "tau_f",
    "gamma_f": "gamma_f",
    "b0": "b0",
    "c0": "c0",
}


def _split_model_variant(model_tag: str) -> tuple[str, str]:
    """Resolve a named method through the shared public variant selector."""
    from multifat.cli.method_selection import resolve_method_tag

    request = resolve_method_tag(model_tag)
    return request.base_name, request.suffix


def classify_loading_type(
    *, H: dict, meta: dict, tol: float = 1e-6, phase_tol_deg: float = 1e-6
) -> str:
    """Label the loading case from the imported stress harmonics and phase metadata.

    Labels produced:
      - "Tensile"
      - "Torsion"
      - "Multiaxial proportional"
      - "Multiaxial non-proportional"

    How it decides:
      * A component counts as engaged if either its mean or its harmonic amplitude is non-zero.
      * Normal stress alone (axial or radial) -> Tensile.
      * Shear stress alone -> Torsion.
      * Several components -> proportional when every relevant phase shift is 0°, else
        non-proportional.
    """
    try:
        mean_t = np.asarray(H.get("S_mean"), dtype=float)
        sin_t = np.asarray(H.get("S_sin"), dtype=float)
        cos_t = np.asarray(H.get("S_cos"), dtype=float)
        harm_amp = np.sqrt(sin_t**2 + cos_t**2)

        # Component slots used by the importer:
        #  radial -> sigma_xx, axial -> sigma_yy, torsion -> tau_xy
        rad_on = (abs(mean_t[0, 0]) > tol) or (abs(harm_amp[0, 0]) > tol)
        ax_on = (abs(mean_t[1, 1]) > tol) or (abs(harm_amp[1, 1]) > tol)
        tor_on = (abs(mean_t[0, 1]) > tol) or (abs(harm_amp[0, 1]) > tol)

        n_on = int(rad_on) + int(ax_on) + int(tor_on)
        if n_on <= 0:
            return "UNKNOWN"

        # Single-component cases
        if n_on == 1:
            if tor_on:
                return "Torsion"
            # axial-only OR radial-only normal stress
            return "Tensile"

        # Several components
        phi_rad = float(meta.get("phase_rad_deg", 0.0) or 0.0) if rad_on else 0.0
        phi_tor = float(meta.get("phase_tor_deg", 0.0) or 0.0) if tor_on else 0.0
        is_nonprop = (abs(phi_rad) > phase_tol_deg) or (abs(phi_tor) > phase_tol_deg)
        return "Multiaxial non-proportional" if is_nonprop else "Multiaxial proportional"
    except Exception:
        return "UNKNOWN"


def resolve_workbook(name: str, data_folder: str | Path = DATA_FOLDER) -> Path:
    """The workbook inside `data_folder` carrying material ID `name` (stem stripped of any ``6_``-style prefix)."""
    matches = [p for p in workbook_files(data_folder) if name in (material_id(p.stem), p.stem)]
    if len(matches) != 1:
        listing = ", ".join(p.name for p in matches) or "none"
        raise WorkbookError(
            f"{name}.xlsx",
            "",
            "",
            f"expected one workbook for material '{name}' in {data_folder}, " f"found {listing}",
        )
    return matches[0]


def workbook_path(material_or_path: str | Path, data_folder: str | Path = DATA_FOLDER) -> Path:
    path = Path(material_or_path)
    if path.suffix.lower() == ".xlsx":
        return path
    return resolve_workbook(str(material_or_path), data_folder)


def load_material(
    material_or_path: str | Path,
    *,
    data_folder: str | Path = DATA_FOLDER,
    verbose: bool = False,
) -> tuple[dict, str]:
    """Read one workbook, without loading the rest of the folder."""
    wb_path = workbook_path(material_or_path, data_folder)
    if not wb_path.exists():
        raise FileNotFoundError(f"Workbook not found: {wb_path}")
    return exp.load_material_workbook(wb_path, verbose=verbose), material_id(wb_path.stem)


def build_params(mat_prop: dict) -> dict:
    """Model parameters taken from the workbook material block (absent entries remain None)."""
    collected = {}
    for param_name, key_name in PARAM_KEYS.items():
        raw = mat_prop.get(key_name, None)
        collected[param_name] = float(raw) if raw is not None else None
    for extra_name, extra_val in (mat_prop.get("extra") or {}).items():
        collected.setdefault(extra_name, float(extra_val))
    return collected


def validate_workbook(fp: str | Path) -> ValidationReport:
    """Contract checks done while loading, plus dataset-level ones (design set, loading types)."""
    fp = Path(fp)
    report = ValidationReport(file=fp.name, material=material_id(fp.stem))
    try:
        workbook = exp.load_material_workbook(fp, report=report)
    except WorkbookError as e:
        report.errors.append(str(e))
        return report
    mat_key = material_id(fp.stem)
    entry = workbook[mat_key]
    if not entry["cases"]:
        from multifat.eval.endurance import is_fatigue_limit_dataset

        if is_fatigue_limit_dataset(entry):
            report.info.append(
                f"fatigue-limit dataset: {len(entry['runouts'])} run-outs judged through the fatigue "
                "index error at the R = -1 fatigue limits (no master curve)"
            )
        elif entry.get("runouts"):
            report.errors.append(
                "only run-outs present, and the Summary holds no usable tensile or torsional fatigue "
                "limits: neither a master curve nor a fatigue index can be built"
            )
        return report
    points = exp.extract_uniaxial_tension(workbook).get(mat_key, {}).get("points", [])
    life_values = {float(pt["Nf_cycles"]) for pt in points}
    report.info.append(f"design set (uniaxial tension, R = -1): {len(points)} tests")
    if len(life_values) < 2:
        report.errors.append(
            "design set needs >= 2 uniaxial tension R = -1 tests with distinct lives "
            f"(found {len(points)}); every model is calibrated on them"
        )
    tally: dict[str, int] = {}
    for record in workbook[mat_key]["cases"].values():
        ltype = classify_loading_type(H=record["harmonics"], meta=record["meta"])
        tally[ltype] = tally.get(ltype, 0) + 1
    report.info.append("loading types: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    return report


@dataclass
class Evaluation:
    """Calibration of one base model on one material, with per-case values across all its variants."""

    material: str
    base_name: str
    params: dict
    R_list: Any
    cal: dict
    values: list[dict] = field(default_factory=list)
    Nexp: list[float] = field(default_factory=list)
    loading_type: list[str] = field(default_factory=list)
    case_ids: list[int] = field(default_factory=list)
    # Identity columns per ("case" | "runout", id), see multifat.report.subsets
    observations: dict = field(default_factory=dict)
    model: Any = None
    aux: list[dict] = field(default_factory=list)
    # Explicit Findley amplitude of the run (see shear_amplitude_for); None for production
    amplitude: dict | None = None
    # Aux values of the design points (native-life fit), lined up with cal["series"]["N"]
    design_aux: list[dict] = field(default_factory=list)
    # Cases that fed the material-constant fit: reported as design points (Flag = 0)
    extra: dict = field(
        default_factory=lambda: {
            "values": [],
            "Nexp": [],
            "aux": [],
            "lt": [],
            "ids": [],
        }
    )
    # Run-outs: scored with the same constants, written out with Flag = 2, never fitted
    runout: dict = field(
        default_factory=lambda: {
            "values": [],
            "Nexp": [],
            "aux": [],
            "lt": [],
            "ids": [],
        }
    )

    @property
    def n_cases(self) -> int:
        return len(self.Nexp)

    def fit(self, variant_key: str):
        return self.cal["fit"] if variant_key == "" else self.cal["fit_ext"]

    def has_variant(self, variant_key: str) -> bool:
        return self.fit(variant_key) is not None and all(variant_key in v for v in self.values)


def check_fit(base_name: str, variant_key: str, cal: dict) -> None:
    curve = cal["fit"] if variant_key == "" else cal["fit_ext"]
    if curve is None:
        # Possible reasons:
        # - an "_ext" flavour was requested on a model that is not CP, OR
        # - no usable uniaxial design points existed for this model (non-positive or NaN values)
        raise ValueError(
            f"No available fit for model '{base_name}' with variant '{variant_key}'. "
            f"fit exists: {cal['fit'] is not None}, "
            f"fit_ext exists: {cal['fit_ext'] is not None}."
        )
    a_surv, slope = curve["A_surv"], curve["b"]
    if not (np.isfinite(a_surv) and a_surv > 0 and np.isfinite(slope)):
        raise ValueError(
            "Invalid master-curve fit: at least two distinct positive experimental lives are required"
        )


def shear_amplitude_for(base_name: str, selection: str) -> dict | None:
    """Explicit Findley amplitude of a run, or None for the production path (method_default)."""
    if selection == "method_default":
        return None
    if base_name != "FIN":
        raise ValueError(
            f"--shear-amplitude-algorithm applies to FIN and FIN_ext only, not {base_name}"
        )
    return {**amplitude_identity(selection), "selected_option": selection}


def harmonic_case_or_skip(case: dict, material: str, base_name: str, case_id: int) -> bool:
    """False, with a RuntimeWarning, when an explicit assessment cannot use this case."""
    try:
        require_stress_harmonics(case, material, case_id)
    except FindleyCalibrationUnavailable as error:
        warnings.warn(
            f"Skipping {material}/{base_name} case {case_id}: {error}",
            RuntimeWarning,
            stacklevel=2,
        )
        return False
    return True


def write_assessment_metadata(
    csv_path: str | Path, *, method: str, result_label: str, amplitude: dict
) -> Path:
    """Sidecar next to an explicit-amplitude CSV: what was selected, resolved and computed."""
    metadata = {
        "schema_version": 1,
        "method": method,
        "result_label": result_label,
        "selected_option": amplitude["selected_option"],
        "resolved_algorithm": amplitude["algorithm"],
        "definition_version": amplitude["definition_version"],
        "evaluation_path": amplitude["evaluation_path"],
        "multifat_version": multifat.__version__,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "csv_sha256": hashlib.sha256(Path(csv_path).read_bytes()).hexdigest(),
    }
    path = Path(csv_path).with_suffix(".meta.json")
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def evaluate_model(
    data: dict,
    material: str,
    base_name: str,
    params: dict,
    R_list=None,
    *,
    require_variant: str | None = None,
    stdnum: float = STDNUM,
    extra_design_ids=(),
    shear_amplitude_algorithm: str = "method_default",
) -> Evaluation:
    """Calibrate on the uniaxial R = -1 design points, then score every remaining case once.

    `extra_design_ids` lists cases that fed the material constants (e.g. the torsion tests
    behind a fitted tau_f, b0): they are reported as design points rather than predictions.
    """
    if R_list is None:
        R_list = build_R_list(n1=DCMESH_N1, n2=DCMESH_N2)
    amplitude = shear_amplitude_for(base_name, shear_amplitude_algorithm)
    mdl = findley_method(shear_amplitude_algorithm) if amplitude else get_models([base_name])[0]

    # Fail early, before the possibly costly calibration.
    absent = [k for k in mdl.required_params() if params.get(k, None) is None]
    if absent:
        raise MissingParamsError(mdl.name, absent)

    native_prediction = capabilities_for_model(mdl).life_prediction == "native"
    uniax_rminus1 = exp.extract_uniaxial_tension(data)
    if material not in uniax_rminus1 and native_prediction:
        # A native-life model needs no axial R=-1 subset or DP master curve.
        uniax_rminus1 = {material: {"points": []}}
    if material not in uniax_rminus1:
        raise ValueError(
            f"No design points for '{material}': the calibration needs uniaxial tension R = -1 tests"
        )
    calibrations = collect_uniaxial_series_and_fit(
        data=data,
        exp=exp,
        material_name=material,
        uniax_cases=uniax_rminus1,
        R_list=R_list,
        models=[mdl],
        params=params,
        stdnum=stdnum,
    )
    # Grab the calibration of this model
    if base_name not in calibrations:
        raise KeyError(
            f"No calibration entry for model '{base_name}'. Available: {list(calibrations.keys())}"
        )
    result = Evaluation(
        material,
        base_name,
        params,
        R_list,
        calibrations[base_name],
        model=mdl,
        amplitude=amplitude,
    )
    result.observations = observation_index(data[material])
    if require_variant is not None and not native_prediction:
        check_fit(base_name, require_variant, result.cal)
    if hasattr(mdl, "life_fit"):
        for pt in uniax_rminus1[material]["points"]:
            cycles = float(pt["Nf_cycles"])
            if not (np.isfinite(cycles) and cycles > 0):
                continue
            s_mean, e_mean, s_amp, e_amp, harm = exp.extract_tensor(
                data, material, int(pt["case_id"]), node=1, include_harmonics=True
            )
            case_res = mdl.evaluate_case(
                S0=s_mean,
                E0=e_mean,
                S1=s_amp,
                E1=e_amp,
                params=params,
                R_list=R_list,
                harmonics=harm,
            )
            val = float(case_res.values.get("", np.nan))
            if np.isfinite(val) and val > 0:
                result.design_aux.append(dict(case_res.aux))
    extra_set = {int(i) for i in extra_design_ids}

    # Skip the design subset and score what is left
    data_without_uniax = exp.exclude_subdataset(data, uniax_rminus1)
    for case_id in data_without_uniax[material]["cases"].keys():
        if amplitude and not harmonic_case_or_skip(
            data_without_uniax[material]["cases"][int(case_id)],
            material,
            base_name,
            int(case_id),
        ):
            continue
        s_mean, e_mean, s_amp, e_amp, harm = exp.extract_tensor(
            data_without_uniax, material, int(case_id), node=1, include_harmonics=True
        )
        case_res = mdl.evaluate_case(
            S0=s_mean,
            E0=e_mean,
            S1=s_amp,
            E1=e_amp,
            params=params,
            R_list=R_list,  # CP uses it, non-CP ignores it
            harmonics=harm,
        )
        if require_variant is not None and require_variant not in case_res.values:
            raise ValueError(
                f"Case result from model '{mdl.name}' missing variant '{require_variant}'. "
                f"Available: {list(case_res.values.keys())}"
            )
        record = data_without_uniax[material]["cases"][int(case_id)]
        ltype = classify_loading_type(H=harm, meta=record.get("meta", {}))
        if int(case_id) in extra_set:
            result.extra["values"].append(dict(case_res.values))
            result.extra["Nexp"].append(float(record["meta"]["Nf_cycles"]))
            result.extra["aux"].append(dict(case_res.aux))
            result.extra["lt"].append(ltype)
            result.extra["ids"].append(int(case_id))
            continue
        result.values.append(dict(case_res.values))
        result.aux.append(dict(case_res.aux))
        result.Nexp.append(float(record["meta"]["Nf_cycles"]))
        result.loading_type.append(ltype)
        result.case_ids.append(int(case_id))

    runout_cases = data[material].get("runouts", {})
    runout_view = {material: {"material": data[material]["material"], "cases": runout_cases}}
    for rid, record in runout_cases.items():
        s_mean, e_mean, s_amp, e_amp, harm = exp.extract_tensor(
            runout_view, material, int(rid), node=1, include_harmonics=True
        )
        case_res = mdl.evaluate_case(
            S0=s_mean,
            E0=e_mean,
            S1=s_amp,
            E1=e_amp,
            params=params,
            R_list=R_list,
            harmonics=harm,
        )
        result.runout["values"].append(dict(case_res.values))
        result.runout["aux"].append(dict(case_res.aux))
        result.runout["Nexp"].append(float(record["meta"]["Nf_cycles"]))
        result.runout["lt"].append(classify_loading_type(H=harm, meta=record.get("meta", {})))
        result.runout["ids"].append(int(rid))
    return result


def native_lives(ev: Evaluation, variant_key: str, F_d, extra_values):
    """Native-life predictions (cases, design points, run-outs), or Nones for models without them."""
    life_fn = getattr(ev.model, "native_life", None)
    if life_fn is None:
        return None, None, None
    exp_lives = np.asarray(ev.cal["series"]["N"], dtype=float)
    if hasattr(ev.model, "life_fit") and not ev.design_aux:
        warnings.warn(
            f"{ev.base_name} native life needs design data for its life coefficients",
            RuntimeWarning,
        )
        nan_value = float("nan")
        return (
            [nan_value] * len(ev.values),
            [nan_value] * len(F_d),
            [nan_value] * len(ev.runout["values"]),
        )
    try:
        life_curve = (
            ev.model.life_fit(ev.design_aux, exp_lives, ev.params) if ev.design_aux else None
        )
    except NativeLifeError as e:
        warnings.warn(f"{ev.base_name} native life not available: {e}", RuntimeWarning)
        nan_value = float("nan")
        return (
            [nan_value] * len(ev.values),
            [nan_value] * len(F_d),
            [nan_value] * len(ev.runout["values"]),
        )
    dp_aux = ev.design_aux or [{} for _ in range(len(F_d) - len(extra_values))]
    dp_aux = dp_aux + ev.extra["aux"]

    def predict(dp, aux):
        return float(life_fn({"dp": float(dp), **aux}, ev.params, life_curve))

    case_lives = [predict(v[variant_key], a) for v, a in zip(ev.values, ev.aux)]
    design_lives = [predict(dp, a) for dp, a in zip(F_d, dp_aux)]
    runout_lives = [
        predict(v[variant_key], a) for v, a in zip(ev.runout["values"], ev.runout["aux"])
    ]
    return case_lives, design_lives, runout_lives


def write_native_variant(
    ev: Evaluation,
    variant_key: str,
    results_dir: str | Path,
) -> str:
    """Export a direct fatigue-life model without inventing a DP master curve.

    The fatigue index is kept only as a diagnostic. The comparison quantity is
    the model's Nf_native, including +inf when the model predicts non-failure.
    """
    series = ev.cal["series"]
    base_name = ev.base_name
    method_name = base_name + variant_key

    design_cycles = np.asarray(series["N"], dtype=float)
    design_index = np.asarray(series["CP"] if variant_key == "" else series["CP_ext"], dtype=float)
    extra_index = [v[variant_key] for v in ev.extra["values"]]
    if extra_index:
        design_cycles = np.concatenate([design_cycles, ev.extra["Nexp"]])
        design_index = np.concatenate([design_index, extra_index])

    cases_index = [float(v[variant_key]) for v in ev.values]
    runout_index = [float(v[variant_key]) for v in ev.runout["values"]]
    native_cases, native_design, native_runouts = native_lives(
        ev, variant_key, design_index, extra_index
    )
    if native_cases is None or native_design is None or native_runouts is None:
        raise ValueError(f"{method_name} declared native life but does not provide it")

    caps = capabilities_for_model(ev.model)

    def expected_dp(cycles):
        n = np.asarray(cycles, dtype=float)
        if not caps.dp_comparison:
            return np.full(n.shape, np.nan)
        if caps.dp_reference == "model_law":
            reference = getattr(ev.model, "expected_dp_at_life", None)
            if reference is None:
                raise ValueError(f"{method_name} lacks its declared DP reference")
            return np.asarray(reference(n, ev.params), dtype=float)
        curve = ev.fit(variant_key)
        if curve is None:
            return np.full(n.shape, np.nan)
        a, b = float(curve["A_surv"]), float(curve["b"])
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return np.where(n > 0, a * np.power(n, b), np.nan)

    dp_fit_cases = expected_dp(ev.Nexp)
    dp_fit_design = expected_dp(design_cycles)
    dp_fit_runouts = expected_dp(ev.runout["Nexp"])
    obs = ev.observations
    ids_design = design_identities(ev, variant_key, len(design_cycles)) if obs else None
    os.makedirs(results_dir, exist_ok=True)
    return save_fatigue_csv(
        out_dir=str(results_dir),
        method_name=method_name,
        material_name=ev.material,
        N_exp=ev.Nexp,
        CP=cases_index,
        N_pred=native_cases,
        N_design=design_cycles,
        F_design=design_index,
        N_pred_design=native_design,
        loading_type=ev.loading_type,
        loading_type_design=["Tensile"] * len(series["N"]) + list(ev.extra["lt"]),
        metric_mode="Nf_direct",
        DP_fit=dp_fit_cases,
        DP_fit_design=dp_fit_design,
        N_native=native_cases,
        N_native_design=native_design,
        N_runout=ev.runout["Nexp"],
        CP_runout=runout_index,
        N_pred_runout=native_runouts,
        DP_fit_runout=dp_fit_runouts,
        loading_type_runout=ev.runout["lt"],
        N_native_runout=native_runouts,
        observation_ids=[obs[("case", i)] for i in ev.case_ids] if obs else None,
        observation_ids_design=ids_design,
        observation_ids_runout=([obs[("runout", i)] for i in ev.runout["ids"]] if obs else None),
    )


def design_identities(ev: Evaluation, variant_key: str, n_rows: int) -> list[dict] | None:
    """Identities of the design rows, pairing the standard series lives with this variant's DP.

    The rows are the standard design series (whose lives serve as the design `N` for every
    variant), followed by the extra design cases; None, with a warning, when the counts disagree.
    """
    case_ids = list(ev.cal.get("design_case_ids", [])) + list(ev.extra["ids"])
    if len(case_ids) != n_rows:
        warnings.warn(
            f"{ev.base_name}{variant_key}: {n_rows} design rows but {len(case_ids)} design cases; "
            "their identity columns are left empty",
            RuntimeWarning,
        )
        return None
    return [ev.observations[("case", i)] for i in case_ids]


def write_variant(
    ev: Evaluation,
    variant_key: str,
    results_dir: str | Path,
    *,
    make_plots: bool = False,
    verbose: bool = False,
) -> str:
    """Error metrics, predicted lives and the CSV for a single variant of an evaluated model."""
    base_name = ev.base_name
    model_name = base_name + variant_key
    if capabilities_for_model(ev.model).life_prediction == "native":
        return write_native_variant(ev, variant_key, results_dir)
    check_fit(base_name, variant_key, ev.cal)
    series = ev.cal["series"]
    curve = ev.fit(variant_key)

    coef_A = curve.get("A", None)
    exponent_b = curve["b"]
    surv_A = curve["A_surv"]

    # ----------------------------
    # Pick the error metric (predicted Nf vs. DP difference)
    # ----------------------------
    # Power law: DP = A * N^b with b = -1/k
    # When k gets very large (nearly flat curve), the predicted N is meaningless.
    # Then switch to the DP log-ratio metric:
    #   Error_ln_dp = ln(DP_point / DP_fit(Nf_exp))
    error_mode = "Nf"
    slope_k = np.inf
    try:
        b_float = float(exponent_b)
        if (not np.isfinite(b_float)) or (b_float >= 0.0) or (b_float == 0.0):
            error_mode = "DP_DIFF"
        else:
            slope_k = -1.0 / b_float
            if np.isfinite(slope_k) and slope_k > K_SLOPE_MAX:
                error_mode = "DP_DIFF"
    except Exception:
        error_mode = "DP_DIFF"

    if verbose:
        print(
            f"Metric mode for {model_name} / {ev.material}: {error_mode} "
            f"(b={exponent_b}, k={slope_k})"
        )

    # Arrays over the design points
    N_d = series["N"]
    F_d = series["CP"] if variant_key == "" else series["CP_ext"]

    # Loading type of the design points (uniaxial tension R=-1 by construction)
    loading_type_design = ["Tensile"] * int(len(np.atleast_1d(N_d)))
    extra_values = [v[variant_key] for v in ev.extra["values"]]
    if extra_values:
        N_d = np.concatenate([np.asarray(N_d, dtype=float), ev.extra["Nexp"]])
        F_d = np.concatenate([np.asarray(F_d, dtype=float), extra_values])
        loading_type_design = loading_type_design + list(ev.extra["lt"])

    if F_d is None or len(np.atleast_1d(F_d)) == 0:
        raise ValueError(
            f"Design series is empty for model '{base_name}' variant '{variant_key}'. "
            "Check that the uniaxial design subset exists and the model produces positive finite values."
        )

    # Predicted cycles at the design points only carry meaning in Nf mode
    if error_mode == "Nf":
        Npred_d = fatc.fatigue_cycles(F_d, surv_A, exponent_b)
    else:
        Npred_d = np.full_like(np.asarray(N_d, dtype=float), np.nan, dtype=float)

    # DP_fit at the design points (always obtainable from the calibration curve)
    N_d_arr = np.asarray(N_d, dtype=float)
    DP_fit_d = np.full_like(N_d_arr, np.nan, dtype=float)
    md = np.isfinite(N_d_arr) & (N_d_arr > 0) & np.isfinite(surv_A) & np.isfinite(exponent_b)
    DP_fit_d[md] = float(surv_A) * (N_d_arr[md] ** float(exponent_b))
    # DP log-ratio error at the design points
    F_d_arr = np.asarray(F_d, dtype=float)
    Error_ln_dp_d = np.full_like(F_d_arr, np.nan, dtype=float)
    md2 = np.isfinite(F_d_arr) & np.isfinite(DP_fit_d) & (F_d_arr > 0) & (DP_fit_d > 0)
    Error_ln_dp_d[md2] = np.log(F_d_arr[md2] / DP_fit_d[md2])

    values: list[float] = []
    Npred: list[float] = []
    DP_fit: list[float] = []
    Error_ln_dp: list[float] = []
    for res_values, Nf in zip(ev.values, ev.Nexp):
        F = float(res_values[variant_key])

        # Calibration-curve DP at the experimental Nf
        dp_fit_i = np.nan
        if np.isfinite(Nf) and Nf > 0 and np.isfinite(surv_A) and np.isfinite(exponent_b):
            dp_fit_i = float(surv_A) * (Nf ** float(exponent_b))

        # Metric branch
        if error_mode == "Nf" and np.isfinite(F) and F > 0:
            N_pred = fatc.fatigue_cycles(F, surv_A, exponent_b)
        else:
            N_pred = np.nan

        # DP log-ratio error: ln(DP_point / DP_fit)
        if np.isfinite(dp_fit_i) and dp_fit_i > 0.0 and np.isfinite(F) and F > 0.0:
            err_ln_dp_i = float(np.log(float(F) / float(dp_fit_i)))
        else:
            err_ln_dp_i = np.nan

        values.append(float(F))
        Npred.append(float(N_pred) if np.isfinite(N_pred) else float("nan"))
        DP_fit.append(float(dp_fit_i) if np.isfinite(dp_fit_i) else float("nan"))
        Error_ln_dp.append(float(err_ln_dp_i) if np.isfinite(err_ln_dp_i) else float("nan"))

    Nexp = list(ev.Nexp)

    N_native, N_native_d, N_native_r = native_lives(ev, variant_key, F_d, extra_values)

    # Run-outs: predicted life plus the master curve at the recorded cycles (a bound), no errors
    CP_r = [float(v[variant_key]) for v in ev.runout["values"]]
    Npred_r = [
        (
            float(fatc.fatigue_cycles(F, surv_A, exponent_b))
            if error_mode == "Nf" and np.isfinite(F) and F > 0
            else float("nan")
        )
        for F in CP_r
    ]
    DP_fit_r = [
        (
            float(surv_A) * (Nf ** float(exponent_b))
            if np.isfinite(Nf) and Nf > 0 and np.isfinite(surv_A) and np.isfinite(exponent_b)
            else float("nan")
        )
        for Nf in ev.runout["Nexp"]
    ]

    # Plot only in Nf mode (the parity plot needs a meaningful N_pred)
    if error_mode == "Nf" and make_plots:
        fatc.fatigue_plots(
            values,
            Nexp,
            Npred,
            A=surv_A,
            B=exponent_b,
            F_design=F_d,
            N_design=N_d,
            F_runout=CP_r,
            N_runout=ev.runout["Nexp"],
            N_pred_runout=Npred_r,
        )

    os.makedirs(results_dir, exist_ok=True)

    obs = ev.observations
    ids_design = design_identities(ev, variant_key, len(np.atleast_1d(N_d))) if obs else None

    basquin_meta = {
        "A": coef_A,
        "b": exponent_b,
        "A_surv": surv_A,
        "STDNUM": STDNUM,
    }
    material_meta = {
        "Sy_MPa": ev.params.get("Sy"),
        "k_FS": ev.params.get("k_FS"),
        "k_FI": ev.params.get("k_FI"),
    }
    extra_meta = {
        "model_name": model_name,
        "base_model": base_name,
        "variant": variant_key,
        "n_planes": len(ev.R_list),
        "dcmesh_n1": DCMESH_N1,
        "dcmesh_n2": DCMESH_N2,
    }

    amplitude = ev.amplitude
    label = f"{base_name}.{amplitude['algorithm']}{variant_key}" if amplitude else None
    csv_path = save_fatigue_csv(
        out_dir=str(results_dir),
        method_name=model_name,  # keep filename compatibility
        result_label=label,
        material_name=ev.material,
        N_exp=Nexp,
        CP=values,
        N_pred=Npred,
        N_design=N_d,
        F_design=F_d,
        N_pred_design=Npred_d,
        loading_type=ev.loading_type,
        loading_type_design=loading_type_design,
        metric_mode=error_mode,
        DP_fit=DP_fit,
        Error_ln_dp=Error_ln_dp,
        DP_fit_design=DP_fit_d,
        Error_ln_dp_design=Error_ln_dp_d,
        basquin_meta=basquin_meta,
        material_meta=material_meta,
        extra_meta=extra_meta,
        N_native=N_native,
        N_native_design=N_native_d,
        N_runout=ev.runout["Nexp"],
        CP_runout=CP_r,
        N_pred_runout=Npred_r,
        DP_fit_runout=DP_fit_r,
        loading_type_runout=ev.runout["lt"],
        N_native_runout=N_native_r,
        observation_ids=[obs[("case", i)] for i in ev.case_ids] if obs else None,
        observation_ids_design=ids_design,
        observation_ids_runout=([obs[("runout", i)] for i in ev.runout["ids"]] if obs else None),
    )
    if amplitude:
        write_assessment_metadata(
            csv_path, method=model_name, result_label=label, amplitude=amplitude
        )
    return csv_path


def main():
    # R_list is genuinely needed only by the CP models; direct models simply ignore it.
    plane_set = build_R_list(n1=DCMESH_N1, n2=DCMESH_N2)
    print("R_list:", np.asarray(plane_set).shape)

    data, material = load_material(MATERIAL_NAME, data_folder=DATA_FOLDER, verbose=True)
    mat_prop = exp.get_material_properties(data, material=material)
    params = build_params(mat_prop)
    print(mat_prop)

    # This run's CSV covers the requested model only. For a CP "_ext" flavour,
    # calibrate the base model (which hands back both variants).
    base_name, variant_key = _split_model_variant(MODEL_NAME)

    if base_name in {"BP", "CAIM"} and (params["Sig0"] is None or params["Tau0"] is None):
        from multifat.models.integral_methods.bohme_papuga import _strengths

        _, _, estimated_s0, estimated_t0 = _strengths(params, base_name)
        print(
            f"{base_name} R=0 strengths (MPa): Sig0={estimated_s0:.3f} "
            f"({'estimated' if params['Sig0'] is None else 'workbook'}), "
            f"Tau0={estimated_t0:.3f} "
            f"({'estimated' if params['Tau0'] is None else 'workbook'})"
        )

    ev = evaluate_model(
        data,
        material,
        base_name,
        params,
        plane_set,
        require_variant=variant_key,
        shear_amplitude_algorithm=SHEAR_AMPLITUDE_ALGORITHM,
    )
    csv_path = write_variant(ev, variant_key, RESULTS_DIR, make_plots=MAKE_PLOTS, verbose=True)
    print("Saved:", csv_path)


def cli(argv=None) -> None:
    """Command-line entry point; the flags take precedence over the module-level settings."""
    import argparse

    global MODEL_NAME, MATERIAL_NAME, RESULTS_DIR, MAKE_PLOTS, SHEAR_AMPLITUDE_ALGORITHM
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--model", default=MODEL_NAME)
    arg_parser.add_argument("--material", default=MATERIAL_NAME)
    arg_parser.add_argument("--results-dir", default=RESULTS_DIR)
    arg_parser.add_argument("--no-plots", action="store_true")
    arg_parser.add_argument(
        "--no-overlay",
        action="store_true",
        help="use the workbooks as supplied, without the overlays",
    )
    arg_parser.add_argument(
        "--shear-amplitude-algorithm",
        default=SHEAR_AMPLITUDE_ALGORITHM,
        choices=("method_default", "mcc_harmonic", "mce_harmonic", "rss"),
        help="Findley (FIN, FIN_ext) only; needs a calibration made with the same option",
    )
    opts = arg_parser.parse_args(argv)
    exp.OVERLAYS_ENABLED = not opts.no_overlay
    SHEAR_AMPLITUDE_ALGORITHM = opts.shear_amplitude_algorithm
    MODEL_NAME, MATERIAL_NAME, RESULTS_DIR = opts.model, opts.material, opts.results_dir
    MAKE_PLOTS = not opts.no_plots
    main()


if __name__ == "__main__":
    cli()
