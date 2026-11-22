"""Single-command pipeline: validate -> calibrate -> run every model -> compare, workbook by workbook.

    multifat-pipeline Experimental_Data/NEW.xlsx
    multifat-pipeline --all --jobs 4

Workbooks are only read, never edited. Outputs (default ``Results/current``):
``<METHOD>_<MATERIAL>.csv`` per model and variant, ``<MATERIAL>_validation.txt``,
``<MATERIAL>_calibration.json``, ``<MATERIAL>_run_report.tsv`` and ``Comparison/``.
A model that cannot run on a given material is reported as skipped while the others carry on.
The calibration JSON is always written: the constants handed over by the pipeline (if any) plus a
``provenance`` block (workbook hash, failure and run-out counts, code revision, numerical settings).
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import platform
import subprocess
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from multifat.cli import run as runmod
import numpy as np

import multifat.imp.experimental_data as exp
from multifat.eval.fitting import stdnum_for_survival
from multifat.eval.shear_amplitude_option import resolved_algorithm
from multifat.imp.schema import WorkbookError, material_id, workbook_files
from multifat.models.errors import (
    DesignSetError,
    LoadingScopeError,
    MaterialScopeError,
    MissingParamsError,
)

DATA_FOLDER = "Experimental_Data"
RESULTS_DIR = "Results/current"


@dataclass
class Outcome:
    model: str
    status: str  # ok | skipped | failed
    reason: str = ""
    n_cases: int = 0
    seconds: float = 0.0
    calibrated: str = ""  # constants taken from the calibration rather than the workbook
    notes: str = ""  # distinct warnings the model raised (e.g. fallback estimates used)
    missing: tuple = ()  # constants the skipped model was missing
    mode: str = "life"  # "fatigue_index" on a fatigue-limit dataset, "out_of_scope" for a skip
    # because the loading lies outside the criterion's scope
    fit: dict | None = None  # master curve of the variant (issue #74); None for native life laws


# Optional constants that a calibrated value may fill; required constants are always filled.
# Fitted Basquin/strain-life constants are never handed to models that read them only
# optionally (GSA and LIU2021 would silently abandon their documented fallbacks).
CALIBRATED_OPTIONAL = {"L_LI", "k_MGSE", "a_MSWT"}
# Optional constants and numerical settings that a constants file may set: the TVEIT2025 fatigue
# limits its A, C, K and L were calibrated with, and its time-integration settings.
SUPPLIED_OPTIONAL = {
    "Sigm1_TVEIT",
    "Taum1_TVEIT",
    "npts_TVEIT",
    "min_blocks_TVEIT",
    "max_blocks_TVEIT",
    "tol_alpha_TVEIT",
    "tol_damage_TVEIT",
}


@dataclass(frozen=True)
class Job:
    """Inputs for evaluating one base model in the pipeline worker."""

    workbook: str
    base_name: str
    results_dir: str
    overrides: dict
    design_ids: dict
    why: dict
    make_plots: bool
    calibration: dict
    suffixes: tuple[str, ...]
    survival_probability: float | None
    design_workbook: str | None = None
    # Master-curve model (--sn-model); an attribute only, so the legacy 11-item tuple form of a
    # Job is unchanged
    sn_model: str = "basquin"

    def __getitem__(self, index: int | slice) -> Any:
        return self._as_tuple()[index]

    def __len__(self) -> int:
        return 11

    def __iter__(self):
        return iter(self._as_tuple())

    def _as_tuple(self) -> tuple[Any, ...]:
        return (
            self.workbook,
            self.base_name,
            self.results_dir,
            self.overrides,
            self.design_ids,
            self.why,
            self.make_plots,
            self.calibration,
            self.suffixes,
            self.survival_probability,
            self.design_workbook,
        )


def every_model() -> list[str]:
    from multifat.models.registry import model_names

    return model_names()


_every_model = every_model


def _evaluate_base(job: Job | tuple) -> list[Outcome]:
    """Run one base model once, writing its std and `_ext` CSVs; warnings become notes."""
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        produced = _evaluate_base_core(job)
    note_text = " | ".join(dict.fromkeys(str(w.message).split(". ")[0] for w in captured))
    for outcome in produced:
        outcome.notes = " | ".join(text for text in (outcome.notes, note_text) if text)
    return produced


def _infeasible_pointer(calibration: dict) -> str:
    """The one-line pointer to a failed plane calibration, appended to notes and reasons."""
    block = calibration.get("block")
    diagnostic = calibration.get("diagnostic")
    if block:
        return f"calibration infeasible at {block}" + (f" ({diagnostic})" if diagnostic else "")
    return f"calibration infeasible ({calibration.get('reason') or 'no finite objective'})"


def _annotate_calibration(outcome: Outcome, model, calibration: dict) -> Outcome:
    """Name the lost constants on the models that ran on their default or were skipped for them.

    A result that adopted a plane constant calibrated on tests outside R = -1 is marked
    ``ok (R fallback)``, with the description of the tests in its notes (issue #77).
    """
    from multifat.cli.calibrate import CALIBRATED

    adopted = set(filter(None, outcome.calibrated.split(",")))
    if (
        outcome.status == "ok"
        and calibration.get("r_fallback")
        and adopted & set(CALIBRATED.values())
    ):
        outcome.status = "ok (R fallback)"
        names = ", ".join(sorted(adopted & set(CALIBRATED.values())))
        flag = f"{names} {calibration['r_fallback']}, not on R = -1 tests only"
        outcome.notes = " | ".join(text for text in (outcome.notes, flag) if text)
    lost = set(calibration.get("lost", ()))
    if not lost:
        return outcome
    pointer = _infeasible_pointer(calibration)
    if outcome.status == "ok":
        defaulted = sorted(k for k in lost if k in getattr(model, "optional_params", {}))
        if defaulted:
            outcome.status = "ok (default constant)"
            outcome.notes = " | ".join(
                text
                for text in (
                    outcome.notes,
                    *(
                        f"{k}: {pointer}; default {model.optional_params[k]} used"
                        for k in defaulted
                    ),
                )
                if text
            )
    elif outcome.status == "skipped" and set(outcome.missing) & lost:
        outcome.reason = "; ".join(
            [outcome.reason, *(f"{k}: {pointer}" for k in sorted(set(outcome.missing) & lost))]
        )
    return outcome


def _evaluate_base_core(job: Job | tuple) -> list[Outcome]:
    if isinstance(job, tuple):
        job = Job(*job)
    workbook = job.workbook
    base_name = job.base_name
    dest = job.results_dir
    overrides = job.overrides
    design_ids = job.design_ids
    why = job.why
    make_plots = job.make_plots
    calibration = job.calibration
    suffixes = job.suffixes
    survival_probability = job.survival_probability
    design_workbook = job.design_workbook
    sn_model = job.sn_model
    start = time.perf_counter()
    adopted: list[str] = []
    try:
        from multifat.models.registry import get_models, stress_only_models

        flagged: dict = {}
        model = get_models([base_name])[0]
        suppliable = _constants_calibration_can_supply(model)
        design_data = None
        if design_workbook is None:
            data, material = runmod.load_material(workbook)
        else:  # cross-batch transfer: constants and design set from design_workbook (#71)
            data, material, design_data = runmod.load_transfer(workbook, design_workbook)
        params = runmod.build_params(runmod.exp.get_material_properties(data, material))
        # A Summary entry the importer flagged (e.g. Sigm1 averaged over stress ratios) is left
        # unused; models that need it are skipped with the importer's own explanation.
        flagged = data[material].get("property_issues", {})
        for key in flagged:
            params[key] = None
        if (params.get("E") is None or params.get("nu") is None) and (
            base_name not in stress_only_models()
        ):
            text = "needs E and nu (absent from the Summary): strains are undefined"
            return [
                _annotate_calibration(
                    Outcome(base_name + suffix, "skipped", text, 0, time.perf_counter() - start),
                    model,
                    calibration,
                )
                for suffix in suffixes
            ]
        for key, value in overrides.items():
            if key in suppliable and params.get(key) is None:
                params[key] = value
                adopted.append(key)
        # Tests underlying a fitted constant become design points of this model, not predictions.
        design_points = sorted({i for key in adopted for i in design_ids.get(key, [])})
        evaluation = runmod.evaluate_model(
            data,
            material,
            base_name,
            params,
            extra_design_ids=design_points,
            survival_probability=survival_probability,
            sn_model=sn_model,
            design_data=design_data,
        )
    except MissingParamsError as exc:
        explanations = sorted({why[key] for key in exc.missing if key in why})
        text = "; ".join(
            [str(exc)]
            + [f"not identifiable: {hint}" for hint in explanations]
            + [f"{key}: {flagged[key]}" for key in exc.missing if key in flagged]
        )
        return [
            _annotate_calibration(
                Outcome(
                    base_name + suffix,
                    "skipped",
                    text,
                    0,
                    time.perf_counter() - start,
                    missing=tuple(exc.missing),
                ),
                model,
                calibration,
            )
            for suffix in suffixes
        ]
    except DesignSetError as exc:
        return [
            Outcome(base_name + suffix, "skipped", str(exc), 0, time.perf_counter() - start)
            for suffix in suffixes
        ]
    except (LoadingScopeError, MaterialScopeError) as exc:
        return [
            Outcome(
                base_name + suffix,
                "skipped",
                str(exc),
                0,
                time.perf_counter() - start,
                mode="out_of_scope",
            )
            for suffix in suffixes
        ]
    except Exception as exc:  # carry on with the remaining models
        text = f"{type(exc).__name__}: {exc} | {traceback.format_exc(limit=-2).strip()}"
        return [
            Outcome(base_name + suffix, "failed", text, 0, time.perf_counter() - start)
            for suffix in suffixes
        ]
    results = []
    fitted_label = ",".join(sorted(adopted))
    from multifat.models.comparison import capabilities_for_model

    native = capabilities_for_model(model).life_prediction == "native"
    for suffix in suffixes:
        label = base_name + suffix
        if suffix == "_ext" and not any("_ext" in v for v in evaluation.values):
            results.append(Outcome(label, "skipped", "model has no _ext variant"))
            continue
        curve = None if native else evaluation.fit(suffix)
        if curve is not None:
            curve = {**curve, "design_loading": evaluation.design_loading}
        try:
            if not native:
                runmod.check_fit(base_name, suffix, evaluation.cal)
        except ValueError as exc:
            results.append(
                Outcome(
                    label,
                    "skipped",
                    f"design-set fit unusable: {exc}",
                    calibrated=fitted_label,
                    fit=curve,
                )
            )
            continue
        try:
            runmod.write_variant(evaluation, suffix, dest, make_plots=make_plots)
            results.append(
                Outcome(
                    label,
                    "ok",
                    "",
                    evaluation.n_cases,
                    time.perf_counter() - start,
                    fitted_label,
                    notes=_life_notes(model, evaluation, suffix, adopted, calibration),
                    fit=curve,
                )
            )
        except (LoadingScopeError, MaterialScopeError) as exc:  # e.g. a run-out out of scope
            results.append(
                Outcome(label, "skipped", str(exc), calibrated=fitted_label, mode="out_of_scope")
            )
        except Exception as exc:
            results.append(
                Outcome(
                    label,
                    "failed",
                    f"{type(exc).__name__}: {exc}",
                    calibrated=fitted_label,
                )
            )
    return [_annotate_calibration(outcome, model, calibration) for outcome in results]


def _life_notes(model, evaluation, suffix: str, adopted, calibration: dict) -> str:
    """Run-report notes on what an "ok" result cannot support: no life, or a bound-limited constant."""
    from multifat.models.comparison import capabilities_for_model

    notes = []
    if capabilities_for_model(model).life_prediction != "native":
        mode, slope_k = runmod.life_mode(evaluation.fit(suffix))
        curve = evaluation.fit(suffix)
        if mode == "Nf" and curve.get("sn_model", "basquin") != "basquin":
            dp = np.array([v.get(suffix, np.nan) for v in evaluation.values], dtype=float)
            n_inf = int(np.sum(np.isinf(runmod.curve_life_at(curve, dp))))
            if n_inf:
                notes.append(
                    f"{n_inf} of {dp.size} assessed failures predicted infinite: DP below the "
                    f"{curve['sn_model']} fatigue limit; their life errors are undefined"
                )
        if mode == "DP_DIFF":
            shape = (
                f"k = -1/b = {slope_k:.3g} > {runmod.K_SLOPE_MAX:g}"
                if math.isfinite(slope_k)
                else "b >= 0"
            )
            notes.append(f"no life predictions: master curve too flat ({shape}); DP errors only")
    for name in sorted(set(adopted) & set(calibration.get("at_bound", {}))):
        low, high = calibration["at_bound"][name]
        notes.append(f"{name} at the search bound of [{low:g}, {high:g}]: not a calibrated optimum")
    return " | ".join(notes)


def _constants_calibration_can_supply(model) -> set[str]:
    return set(model.required_params()) | (
        set(getattr(model, "optional_params", {})) & (CALIBRATED_OPTIONAL | SUPPLIED_OPTIONAL)
    )


def _absent_suppliable_constants(
    workbook: Path, names: list[str], data: dict | None = None
) -> set[str]:
    """Constants a selected model needs, calibration can supply, and the workbook does not have."""
    from multifat.cli.calibrate import CALIBRATED
    from multifat.eval.material_fits import BASQUIN_KEYS, OP_CURVE_KEYS, PLASTIC_KEYS
    from multifat.models.registry import get_models, stress_only_models

    if data is None:
        data, material = runmod.load_material(workbook)
    else:
        material = material_id(workbook.stem)
    params = runmod.build_params(runmod.exp.get_material_properties(data, material))
    supplied = set(CALIBRATED.values()) | set(BASQUIN_KEYS) | set(PLASTIC_KEYS) | set(OP_CURVE_KEYS)
    absent = {name for name in supplied if params.get(name) is None}
    if params.get("E") is None or params.get("nu") is None:  # strain models are skipped anyway
        names = [name for name in names if name in stress_only_models()]
    return set().union(
        set(),
        *(absent & _constants_calibration_can_supply(m) for m in get_models(names)),
    )


def _repo_snapshot() -> dict:
    def run_git(*flags):
        proc = subprocess.run(["git", *flags], capture_output=True, text=True, timeout=10)
        return proc.stdout.strip() if proc.returncode == 0 else None

    try:
        return {
            "commit": run_git("rev-parse", "HEAD"),
            "dirty": bool(run_git("status", "--porcelain")),
        }
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None}


def load_constants(path: Path, material: str, n_failures: int | None = None) -> tuple[dict, dict]:
    """Supplied constants from an upstream calibration record or from the sidecar written here.

    An upstream record (``Calibrations/<material>/<id>.json``, A. Chiocca, dc4fe6d) keeps its
    constants under ``parameters``; its a_MSWT uses the weighting of the former workbooks and is
    converted with a' = a/(2 - a) (Jiang 2000 Eq. 6, as in 23e1dfe). The sidecar written here keeps
    them in ``params``. A published-constants file keeps constants printed in a paper under
    ``published_constants``, with ``material_id`` and ``reference``; a hand-entered file keeps
    constants calibrated outside the pipeline under ``hand_entered_constants``, with
    ``material_id`` and the ``provenance`` of the entry (#60). Returns (constants,
    description of the source for the new sidecar).
    """
    path = Path(path)
    blob = path.read_bytes()
    record = json.loads(blob.decode("utf-8"))
    origin = {"file": path.as_posix(), "sha256": hashlib.sha256(blob).hexdigest()}
    if "parameters" in record and "schema_version" in record:
        declared_material = record.get("material_id")
        constants = {k: float(v) for k, v in record["parameters"].items() if v is not None}
        origin.update(
            kind="upstream calibration record",
            calibration_id=record.get("calibration_id"),
        )
        from multifat.models.references import IMPLEMENTATION_REVISION

        warnings_list: list[str] = []
        recorded_failures = (record.get("counts") or {}).get("failure")
        if (
            n_failures is not None
            and recorded_failures is not None
            and recorded_failures != n_failures
        ):
            warnings_list.append(
                f"record made on {recorded_failures} failures, the workbook has {n_failures}: "
                "its constants were not calibrated on these data"
            )
        recorded_revision = record.get("implementation_revision")
        if recorded_revision is not None and recorded_revision != IMPLEMENTATION_REVISION:
            warnings_list.append(
                f"made with implementation revision {recorded_revision}, current is "
                f"{IMPLEMENTATION_REVISION}: its constants belong to earlier formulations"
            )
        if warnings_list:
            origin["warning"] = "; ".join(warnings_list)
        if "a_MSWT" in constants:
            weight = constants["a_MSWT"]
            constants["a_MSWT"] = weight / (2.0 - weight)
            origin["note"] = f"a_MSWT {weight} converted to the printed weighting a/(2 - a)"
    elif "hand_entered_constants" in record:
        # Constants entered by hand from an external or article-style calibration (TVEIT2025
        # stage 1, #60): not a pipeline calibration, so no fingerprint or revision check; the
        # record names its material, and its provenance text travels into the run sidecar.
        declared_material = record.get("material_id")
        constants = {k: float(v) for k, v in record["hand_entered_constants"].items()}
        origin.update(kind="hand-entered constants", provenance=record.get("provenance"))
    elif "published_constants" in record:
        # Constants printed in a publication (e.g. the SSF surface of Anes et al. 2014): not a
        # calibration of these data, so no fingerprint; the record names its source and material.
        declared_material = record.get("material_id")
        constants = {k: float(v) for k, v in record["published_constants"].items()}
        origin.update(kind="published constants", reference=record.get("reference"))
    elif "params" in record:
        declared_material = (
            material_id(Path(record.get("provenance", {}).get("workbook", "")).stem) or None
        )
        constants = {k: float(v) for k, v in record["params"].items() if v is not None}
        origin.update(kind="multifat calibration sidecar")
    else:
        raise ValueError(f"{path}: neither a calibration record nor a calibration sidecar")
    if declared_material is not None and declared_material != material:
        raise ValueError(
            f"{path} holds constants for material '{declared_material}', not '{material}'"
        )
    from multifat.models.registry import get_models

    for method_name, setting in (record.get("amplitude_settings") or {}).items():
        # The pipeline evaluates every method with its own amplitude (method_default; Findley
        # MCE, A. Chiocca 36ff0a1): a constant fitted with another amplitude would be misapplied.
        algorithm = (setting or {}).get("algorithm")
        expected = resolved_algorithm(get_models([method_name])[0])
        if algorithm not in (None, expected):
            raise ValueError(
                f"{path}: the {method_name} constant was calibrated with the shear amplitude "
                f"'{algorithm}', but the pipeline evaluates {method_name} with {expected}"
            )
    return constants, origin


def is_published_constants(path: Path, material: str) -> bool:
    """True for a published-constants record (``published_constants``) of `material`.

    The one kind of supplied record a cross-batch transfer accepts (issue #71): constants printed
    in a paper for the design batch, such as the SSF surface of 42CrMo4_Reis. Calibration records
    and sidecars are refused there, since the transfer calibrates the design batch itself.
    """
    try:
        record = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return "published_constants" in record and record.get("material_id") == material


def supplied_design_ids(path: Path, data_block: dict, supplied: dict) -> tuple[dict, str]:
    """Design tests of supplied constants, from a sidecar made on the same failures.

    A sidecar written by this pipeline records the tests behind each fitted constant. They are
    reused only when its ``failure_fingerprint`` equals that of the current data, since case
    indices are not stable across workbook revisions. Upstream records carry no case indices.
    Returns ({constant: case ids}, note for the sidecar or "").
    """
    from multifat.report.subsets import dataset_fingerprint

    record = json.loads(Path(path).read_text(encoding="utf-8"))
    recorded = record.get("design_case_ids") or {}
    wanted = {k: [int(i) for i in v] for k, v in recorded.items() if k in supplied}
    unknown = "the tests behind the supplied constants are unknown; they stay Flag = 1"
    if "params" not in record:
        return {}, unknown
    fingerprint = (record.get("provenance") or {}).get("failure_fingerprint")
    if fingerprint != dataset_fingerprint(data_block):
        return {}, "sidecar made on other failures: " + unknown
    return wanted, ""


def provenance(
    workbook: Path,
    data: dict | None = None,
    survival_probability: float = 0.5,
    sn_model: str = "basquin",
) -> dict:
    """What produced a run: the input file, the observations, the code and the numerical settings."""
    import numpy as np

    import multifat
    from multifat.models.references import IMPLEMENTATION_REVISION
    from multifat.report.subsets import dataset_fingerprint

    if data is None:
        data, material = runmod.load_material(workbook)
    else:
        material = material_id(workbook.stem)
    return {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "workbook": workbook.as_posix(),
        "workbook_sha256": hashlib.sha256(workbook.read_bytes()).hexdigest(),
        "failures": len(data[material]["cases"]),
        "runouts": len(data[material].get("runouts", {})),
        # Overlay identity (name, sha256, operations) applied to the workbook, or None
        "assumptions_applied": data[material].get("overlay"),
        # Failure-only hash of the properties and stress histories (multifat.report.subsets); it
        # is not comparable with the input_fingerprint of upstream records, which hashes other
        # fields
        "failure_fingerprint": dataset_fingerprint(data[material]),
        "multifat_version": multifat.__version__,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "git": _repo_snapshot(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "settings": {
            "plane_grid": [runmod.DCMESH_N1, runmod.DCMESH_N2],
            "refined_plane_search": True,
            "master_curve_stdnum": stdnum_for_survival(survival_probability),
            "master_curve_survival_probability": survival_probability,
            "master_curve_model": sn_model,
            "flat_curve_k_max": runmod.K_SLOPE_MAX,
        },
    }


def _calibration_state(sidecar: dict, needed: set[str]) -> dict:
    """What an infeasible calibration lost, for the models that depended on the lost constants.

    Also carries ``r_fallback``, the description of a calibration made on tests outside
    R = -1 (issue #77), so that the models adopting those constants are flagged.
    """
    r_fallback = (sidecar.get("selection") or {}).get("r_fallback", "")
    extra = {"r_fallback": r_fallback} if r_fallback else {}
    if sidecar.get("status") != "infeasible":
        state = {"at_bound": sidecar.get("at_bound", {})} if sidecar.get("at_bound") else {}
        return {**state, **extra}
    from multifat.cli.calibrate import CALIBRATED

    info = sidecar.get("infeasible") or {}
    lost = sorted(set(CALIBRATED.values()) & set(needed) - set(sidecar.get("params", {})))
    return {
        "status": "infeasible",
        "block": info.get("block"),
        "diagnostic": info.get("diagnostic") or "",
        "reason": sidecar.get("reason", ""),
        "lost": lost,
        "at_bound": sidecar.get("at_bound", {}),
        **extra,
    }


def run_workbook(
    workbook: Path,
    *,
    models: list[str] | None = None,
    all_variants: bool = False,
    results_dir: Path = Path(RESULTS_DIR),
    jobs: int = 1,
    calibrate: bool = True,
    plots: bool = False,
    constants: Path | None = None,
    survival_probability: float = 0.5,
    sn_model: str = "basquin",
    calibration_archive: Path | None = None,
    design_workbook: Path | None = None,
) -> list[Outcome]:
    """Run the selected models on one workbook.

    With `design_workbook` (cross-batch transfer, issue #71) the failures of `workbook` are
    predicted from another batch: the Summary properties, the calibrated constants, the design
    set and the master curves all come from `design_workbook`, and the outputs are written under
    ``<material>_from_<design material>``. The design-side tests appear as design rows (Flag 0);
    every failure of `workbook` is assessed (Flag 1).
    """
    stdnum_for_survival(survival_probability)  # reject invalid values before any work
    results_dir.mkdir(parents=True, exist_ok=True)
    material = material_id(workbook.stem)
    target_workbook = workbook
    constants_material = material
    if design_workbook is not None:
        constants_material = material_id(design_workbook.stem)
        if constants is not None and not is_published_constants(constants, constants_material):
            raise ValueError(
                "a transfer takes its constants from the design workbook; the only supplied "
                "constants it accepts are a published-constants record of the design material "
                f"'{constants_material}'"
            )
        material = runmod.transfer_material_id(material, constants_material)
        # the design side must hold a valid design set; the target only needs to be readable
        workbook = design_workbook

    report = runmod.validate_workbook(workbook)
    (results_dir / f"{material}_validation.txt").write_text(report.text(), encoding="utf-8")
    if not report.ok:
        raise WorkbookError(workbook.name, "", "", "; ".join(report.errors))

    from multifat.cli.method_selection import select_method_variants

    requested = select_method_variants(models, all_variants=all_variants)
    selected = list(dict.fromkeys(request.base_name for request in requested))
    if design_workbook is None:
        data, _ = runmod.load_material(workbook)
    else:
        data, _, _ = runmod.load_transfer(target_workbook, design_workbook)
    # Explicit workflow and observation roles from the shared importer (issue #18)
    from multifat.imp.classification import classify_dataset

    classification = classify_dataset(data[material], material)
    (results_dir / f"{material}_classification.json").write_text(
        json.dumps(classification, indent=1, default=float), encoding="utf-8"
    )
    if classification["workflow"] not in ("fatigue_life", "fatigue_limit"):
        raise WorkbookError(workbook.name, "", "", "; ".join(classification["reasons"]))
    if classification["workflow"] == "fatigue_limit" and design_workbook is not None:
        raise WorkbookError(workbook.name, "", "", "a transfer needs fatigue-life datasets")
    if classification["workflow"] == "fatigue_limit":
        results, sidecar = _evaluate_fatigue_limit(workbook, material, requested, results_dir, jobs)
        sidecar["provenance"] = provenance(workbook, data, survival_probability, sn_model)
        sidecar["status"] = "not_applicable"
        sidecar["reason"] = "fatigue-limit dataset: constants from the R = -1 limits"
        if constants:
            sidecar["constants_source"] = {
                "file": Path(constants).as_posix(),
                "note": "not used: fatigue-limit dataset (constants are taken from the limits)",
            }
        (results_dir / f"{material}_calibration.json").write_text(
            json.dumps(sidecar, indent=1, default=float), encoding="utf-8"
        )
        _dump_run_report(results_dir, material, results)
        return results
    n_fail = len(data[material]["cases"])
    supplied, source = (
        load_constants(constants, constants_material, n_fail) if constants else ({}, None)
    )
    if source is not None and source.get("kind") == "multifat calibration sidecar":
        from multifat.cli.records import check_record_compatibility, record_warnings

        record = json.loads(Path(constants).read_text(encoding="utf-8"))
        problems = check_record_compatibility(record, data[material])
        if problems:
            raise ValueError(f"{constants}: not reusable here: {'; '.join(problems)}")
        notes = record_warnings(record, data[material])
        if record.get("calibration_id"):
            source["calibration_id"] = record["calibration_id"]
        if notes:
            source["warning"] = "; ".join(
                [source["warning"], *notes] if "warning" in source else notes
            )
    needed = (
        (
            _absent_suppliable_constants(
                workbook, selected, data if design_workbook is None else None
            )
            - set(supplied)
        )
        if calibrate
        else set()
    )
    sidecar: dict = {"params": {}}
    if needed:
        from multifat.cli.calibrate import CALIBRATED, calibrate_workbook

        sidecar = calibrate_workbook(
            workbook,
            planes=bool(needed & set(CALIBRATED.values())),
            # k_MGSE comes from Zhu (2018) Eq. (12) on the Basquin fits whenever they are
            # identifiable, so it must not depend on which other models were requested.
            constants=bool(needed - set(CALIBRATED.values())) or "k_MGSE" in needed,
        )
    elif calibrate:
        sidecar["status"] = "not_applicable"
        sidecar["reason"] = (
            "nothing to calibrate: every constant these models need is in the workbook "
            "or in the supplied record"
        )
    else:
        sidecar["status"] = "not_applicable"
        sidecar["reason"] = "calibration not requested (--no-calibrate)"
    if supplied:
        sidecar["params"] = {**sidecar.get("params", {}), **supplied}
        reused, note = supplied_design_ids(constants, data[material], supplied)
        sidecar["design_case_ids"] = {**sidecar.get("design_case_ids", {}), **reused}
        if note:
            source = {**source, "design_rows": note}
        sidecar["constants_source"] = source
    if design_workbook is None:
        sidecar["provenance"] = provenance(workbook, data, survival_probability, sn_model)
    else:  # the constants describe the design batch; the assessed tests are the target's
        sidecar["provenance"] = provenance(workbook, None, survival_probability, sn_model)
        sidecar["transfer"] = {
            "design_material": material_id(design_workbook.stem),
            "target_material": material_id(target_workbook.stem),
            "target_provenance": provenance(target_workbook, None, survival_probability, sn_model),
        }
    from multifat.cli.records import CALIBRATION_PROCEDURE, archive_record

    sidecar["calibration_procedure"] = CALIBRATION_PROCEDURE
    if sidecar.get("params"):
        # A persistent, never-overwritten copy; the results refer to it by calibration_id
        archive = calibration_archive or results_dir / "calibrations"
        sidecar["calibration_id"], _, _ = archive_record(sidecar, archive, material)
    (results_dir / f"{material}_calibration.json").write_text(
        json.dumps(sidecar, indent=1, default=float), encoding="utf-8"
    )
    overrides = sidecar.get("params", {})
    design_ids = sidecar.get("design_case_ids", {})
    why = sidecar.get("unidentified", {})
    calibration = _calibration_state(sidecar, needed)

    suffixes_by_base: dict[str, list[str]] = {}
    for request in requested:
        suffixes_by_base.setdefault(request.base_name, []).append(request.suffix)
    queue = [
        Job(
            workbook=str(target_workbook),
            base_name=base_name,
            results_dir=str(results_dir),
            overrides=overrides,
            design_ids=design_ids,
            why=why,
            make_plots=plots,
            calibration=calibration,
            suffixes=tuple(suffixes),
            survival_probability=survival_probability,
            design_workbook=None if design_workbook is None else str(design_workbook),
            sn_model=sn_model,
        )
        for base_name, suffixes in suffixes_by_base.items()
    ]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as pool:
            results = [item for chunk in pool.map(_evaluate_base, queue) for item in chunk]
    else:
        results = [item for job in queue for item in _evaluate_base(job)]

    _dump_run_report(results_dir, material, results, sidecar.get("calibration_id", ""))
    _dump_master_curves(results_dir, material, results)
    if classification["roles"]["limit_state"] and design_workbook is None:
        _evaluate_declared_limits(workbook, material, requested, results_dir, jobs)
    return results


def _evaluate_declared_limits(workbook: Path, material: str, requests, dest: Path, jobs: int):
    """Fatigue index of the declared fatigue limits of a mixed life workbook (d113, issue #18).

    The life outcomes stay the workbook's result (run report, eligibility); the limit states are
    reported apart, in ``FatigueLimit/`` with their own run report and endurance constants.
    """
    results, sidecar = _evaluate_fatigue_limit(workbook, material, requests, dest, jobs)
    _dump_run_report(dest / "FatigueLimit", material, results)
    (dest / "FatigueLimit" / f"{material}_fatigue_limit.json").write_text(
        json.dumps(sidecar, indent=1, default=float), encoding="utf-8"
    )
    return results


MASTER_CURVE_COLUMNS = [
    "model", "status", "design_loading", "n", "decades", "A", "b", "b_se", "R2_log",
    "STD_ln_DP", "STD", "identifiable", "reason",
]  # fmt: skip


def _dump_master_curves(dest: Path, material: str, results: list[Outcome]) -> None:
    """Fit of the master curve behind every result: ``<material>_master_curves.csv`` (#74).

    One row per model variant that fitted a master curve, also when it was skipped because the
    curve is not identifiable, so that the reader can judge every curve the results rest on.
    """
    import csv

    rows = [o for o in results if o.fit]
    with (dest / f"{material}_master_curves.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(MASTER_CURVE_COLUMNS)
        for outcome in rows:
            fit = outcome.fit
            writer.writerow(
                [outcome.model, outcome.status, fit.get("design_loading", "")]
                + [
                    f"{fit[k]:.6g}" if isinstance(fit.get(k), float) else fit.get(k, "")
                    for k in MASTER_CURVE_COLUMNS[3:]
                ]
            )


def _dump_run_report(
    dest: Path, material: str, results: list[Outcome], calibration_id: str = ""
) -> None:
    sep = "\t"
    # calibration_id: the archived calibration record the constants come from (issue #15)
    head = ["model", "status", "reason", "calibrated", "n_cases", "seconds", "notes"]
    head.append("calibration_id")
    lines = [sep.join(head)]
    for outcome in results:
        row = [
            outcome.model,
            outcome.status,
            " ".join(outcome.reason.split()),
            outcome.calibrated,
        ]
        row += [
            str(outcome.n_cases),
            f"{outcome.seconds:.1f}",
            " ".join(outcome.notes.split()),
            calibration_id,
        ]
        lines.append(sep.join(row))
    (dest / f"{material}_run_report.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


FI_COLUMNS = [
    "Material", "Method", "Loading_type", "Specimen_id", "N_ref", "DP", "DP_ref_axial",
    "DP_ref_torsion", "FI_axial", "FI_torsion", "FI_load", "Implementation_revision",
    "Method_reference",
]  # fmt: skip


def _fatigue_limit_model(job: tuple) -> tuple[Outcome, dict]:
    """Fatigue index error of one base model on a fatigue-limit dataset (decision d057)."""
    import csv

    from multifat.eval.endurance import (
        ENDURANCE_CONSTANTS,
        calibrate_endurance_constant,
        case_from_block,
        fatigue_index,
        reference_cases,
    )
    from multifat.models.references import IMPLEMENTATION_REVISION, method_reference
    from multifat.models.registry import get_models, stress_only_models

    workbook, base_name, suffix, dest = job
    model_tag = base_name + suffix
    start = time.perf_counter()
    adopted: dict = {}
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        try:
            (model,) = get_models([base_name])
            data, material = runmod.load_material(workbook)
            block = data[material]
            params = runmod.build_params(runmod.exp.get_material_properties(data, material))
            if (params.get("E") is None or params.get("nu") is None) and (
                base_name not in stress_only_models()
            ):
                text = "needs E and nu (absent from the Summary): strains are undefined"
                return (
                    Outcome(model_tag, "skipped", text, 0, time.perf_counter() - start),
                    adopted,
                )
            ratios = runmod.build_R_list(n1=runmod.DCMESH_N1, n2=runmod.DCMESH_N2)
            refs = reference_cases(params["Sigm1"], params["Taum1"], params["E"], params["nu"])
            for param_name in sorted(model.required_params()):
                if params.get(param_name) is None and param_name in ENDURANCE_CONSTANTS:
                    value = calibrate_endurance_constant(
                        model,
                        param_name,
                        params,
                        refs,
                        ENDURANCE_CONSTANTS[param_name],
                        ratios,
                        variant=suffix,
                    )
                    if value is None:
                        text = f"{param_name}: no value in {ENDURANCE_CONSTANTS[param_name]} equalises the "
                        text += "axial and torsional fatigue limits"
                        return (
                            Outcome(
                                model_tag,
                                "skipped",
                                text,
                                0,
                                time.perf_counter() - start,
                            ),
                            adopted,
                        )
                    params[param_name] = adopted[param_name] = value
            absent = [name for name in model.required_params() if params.get(name) is None]
            if absent:
                raise MissingParamsError(model.name, absent)
            from multifat.imp.classification import limit_state_runouts

            # the limit states only: declared interrupted tests are not limits (d113, #18)
            limits = list(limit_state_runouts(block, material).values())
            rows = fatigue_index(
                model,
                params,
                [case_from_block(c) for c in limits],
                refs,
                ratios,
                variant=suffix,
            )
            citation = method_reference(base_name)
            target = Path(dest) / "FatigueLimit" / f"{model_tag}_{material}_FI.csv"
            target.parent.mkdir(parents=True, exist_ok=True)
            with open(target, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(FI_COLUMNS)
                for case, row in zip(limits, rows):
                    loading = runmod.classify_loading_type(H=case["harmonics"], meta=case["meta"])
                    writer.writerow(
                        [material, model_tag, loading, case["meta"].get("specimen_id", ""),
                         case["meta"]["Nf_cycles"], row["DP"], row["DP_ref_axial"], row["DP_ref_torsion"],
                         row["FI_axial"], row["FI_torsion"], row["FI_load"], IMPLEMENTATION_REVISION,
                         citation["url"] or citation["citation"]]
                    )  # fmt: skip
            outcome = Outcome(
                model_tag,
                "ok",
                "",
                len(rows),
                time.perf_counter() - start,
                ",".join(adopted),
            )
            outcome.mode = "fatigue_index"
        except MissingParamsError as exc:
            outcome = Outcome(
                model_tag,
                "skipped",
                str(exc),
                0,
                time.perf_counter() - start,
                missing=tuple(exc.missing),
            )
        except (LoadingScopeError, MaterialScopeError) as exc:
            outcome = Outcome(model_tag, "skipped", str(exc), 0, time.perf_counter() - start)
            outcome.mode = "out_of_scope"
        except Exception as exc:  # carry on with the remaining models
            outcome = Outcome(
                model_tag,
                "failed",
                f"{type(exc).__name__}: {exc}",
                0,
                time.perf_counter() - start,
            )
    outcome.notes = " | ".join(
        dict.fromkeys(str(warning.message).split(". ")[0] for warning in captured)
    )
    return outcome, adopted


def _evaluate_fatigue_limit(workbook: Path, material: str, requests, dest: Path, jobs: int):
    """Every model on a fatigue-limit dataset; returns the outcomes and the sidecar payload."""
    import pandas as pd

    queue = [(str(workbook), request.base_name, request.suffix, str(dest)) for request in requests]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as pool:
            pairs = list(pool.map(_fatigue_limit_model, queue))
    else:
        pairs = [_fatigue_limit_model(job) for job in queue]
    frames = [
        pd.read_csv(dest / "FatigueLimit" / f"{outcome.model}_{material}_FI.csv")
        for outcome, _ in pairs
        if outcome.status == "ok"
    ]
    if frames:
        from multifat.eval.endurance import fatigue_index_summary

        summary = fatigue_index_summary(pd.concat(frames), material)
        summary.to_csv(dest / "FatigueLimit" / f"{material}_fatigue_index_summary.csv", index=False)
    sidecar = {
        "mode": "fatigue-limit dataset: fatigue index error (decision d057)",
        "endurance_constants": {outcome.model: adopted for outcome, adopted in pairs if adopted},
    }
    return [outcome for outcome, _ in pairs], sidecar


def _calibratable_keys() -> set[str]:
    """Constants the pipeline can calibrate or fit when a workbook does not provide them."""
    from multifat.cli.calibrate import CALIBRATED
    from multifat.eval.material_fits import BASQUIN_KEYS, OP_CURVE_KEYS, PLASTIC_KEYS

    return set(CALIBRATED.values()) | set(BASQUIN_KEYS) | set(PLASTIC_KEYS) | set(OP_CURVE_KEYS)


def _eligibility_status(outcome: Outcome, calibrate: bool) -> str:
    if outcome.status.startswith("ok") and outcome.mode == "fatigue_index":
        return "fatigue_index"
    if outcome.status.startswith("ok"):
        return "available" if outcome.n_cases else "prediction_only"
    if outcome.mode == "out_of_scope":
        return "out_of_scope"
    if (
        outcome.status == "skipped"
        and outcome.missing
        and not calibrate
        and set(outcome.missing) <= _calibratable_keys()
    ):
        return "calibration_required"
    return "insufficient_inputs" if outcome.status == "skipped" else "failed"


def write_eligibility(runs: dict, out_dir: Path, *, calibrate: bool, models=None) -> None:
    """Method x material status matrix, per-method counts and reasons (upstream 904c619).

    `runs` maps a material to its base-model outcomes, to ``{"outcomes": ..., "calibration":
    {"status", "reason"}}`` for a calibrated run, or to (status, reason) for a rejected
    workbook. The status records what this run produced; it is not an accuracy ranking. The
    calibration status of each material is written next to the method rows and in
    ``calibration_status.csv``.
    """
    import pandas as pd

    rows = []
    calibration_rows = []
    for material, record in runs.items():
        if isinstance(record, dict):
            calibration = record.get("calibration") or {}
            outcomes = record.get("outcomes", [])
        elif isinstance(record, tuple):
            calibration = {"status": record[0], "reason": record[1]}
            outcomes = []
        else:
            calibration = {}
            outcomes = record
        calibration_rows.append(
            dict(
                material=material,
                status=calibration.get("status", ""),
                reason=calibration.get("reason", ""),
            )
        )
        if isinstance(record, tuple):
            for method in models or _every_model():
                rows.append(
                    dict(
                        method=method,
                        material=material,
                        status=record[0],
                        reason=record[1],
                        calibration_status=calibration["status"],
                    )
                )
            continue
        for outcome in outcomes:
            if outcome.model.endswith("_ext") and any(
                item.model == outcome.model.removesuffix("_ext") for item in outcomes
            ):
                continue
            state = _eligibility_status(outcome, calibrate)
            count = (
                0 if outcome.mode == "fatigue_index" else outcome.n_cases
            )  # limit states are not failures
            rows.append(
                dict(
                    method=outcome.model,
                    material=material,
                    status=state,
                    reason=outcome.reason,
                    n_cases=count,
                    calibration_status=calibration.get("status", ""),
                )
            )
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    (out_dir / "eligibility.json").write_text(
        json.dumps(rows, indent=1, default=str), encoding="utf-8"
    )
    pd.DataFrame(calibration_rows).to_csv(out_dir / "calibration_status.csv", index=False)
    frame.pivot(index="method", columns="material", values="status").to_csv(
        out_dir / "eligibility_matrix.csv"
    )
    frame["available"] = frame["status"] == "available"
    frame.groupby("method").agg(
        materials=("available", "sum"),
        failure_observations=("n_cases", lambda series: int(series.fillna(0).sum())),
        excluded_combinations=("available", lambda series: int((~series).sum())),
    ).to_csv(out_dir / "method_availability.csv")


def _rejection_record(workbook: Path, failure: Exception) -> tuple[str, str]:
    report = runmod.validate_workbook(workbook)
    only_notched = any(s.startswith("rows dropped (notched)") for s in report.info) and any(
        "kept: 0 failures, 0 run-outs" in s for s in report.info
    )
    return ("out_of_scope" if only_notched else "insufficient_inputs"), str(failure)


def _calibration_summary(results_dir: Path, material: str) -> dict:
    """The status and reason the run wrote to ``<MATERIAL>_calibration.json``."""
    path = Path(results_dir) / f"{material}_calibration.json"
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {"status": record.get("status", ""), "reason": record.get("reason", "")}


def main(cli_args=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workbooks", nargs="*", type=Path, help="workbook paths (.xlsx)")
    parser.add_argument(
        "--material", action="append", default=[], help="name inside Experimental_Data"
    )
    parser.add_argument("--all", action="store_true", help="every workbook in Experimental_Data")
    parser.add_argument(
        "--models",
        default="all",
        help="comma-separated method names, including explicit _ext, or 'all'",
    )
    parser.add_argument(
        "--all-variants",
        action="store_true",
        help="also evaluate each selected critical-plane method's available _ext variant",
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--results-dir", type=Path, default=Path(RESULTS_DIR))
    parser.add_argument("--no-calibrate", action="store_true")
    parser.add_argument("--no-compare", action="store_true")
    parser.add_argument("--plots", action="store_true")
    parser.add_argument(
        "--sn-model",
        default="basquin",
        choices=runmod.SN_MODELS,
        help="master curve converting the damage parameter to life; native life laws and the "
        "calibration of plane constants are not affected",
    )
    parser.add_argument(
        "--calibration-archive",
        type=Path,
        default=None,
        help="folder of the persistent calibration records (default: <results-dir>/calibrations)",
    )
    parser.add_argument(
        "--calibration-id",
        default=None,
        help="reuse the archived calibration record ID of the (single) workbook's material",
    )
    parser.add_argument(
        "--survival-probability",
        type=float,
        default=0.5,
        help="survival probability P_s of the master curves, 0.5 <= P_s < 1 (0.5: median); "
        "native life laws are not shifted",
    )
    parser.add_argument(
        "--constants",
        type=Path,
        default=None,
        help="constants from a Calibrations/ record or a calibration sidecar (one workbook)",
    )
    parser.add_argument(
        "--design-material",
        default=None,
        help="cross-batch transfer (issue #71): predict every given workbook from this material, "
        "whose Summary, calibration and design set are used; outputs <material>_from_<design>",
    )
    parser.add_argument(
        "--no-overlay",
        action="store_true",
        help="use the workbooks as supplied, without the overlays",
    )
    args = parser.parse_args(cli_args)
    exp.OVERLAYS_ENABLED = not args.no_overlay

    workbooks = list(args.workbooks)
    try:
        workbooks += [runmod.resolve_workbook(name, DATA_FOLDER) for name in args.material]
    except WorkbookError as exc:
        parser.error(str(exc))
    if args.all:
        workbooks += workbook_files(DATA_FOLDER)
    if not workbooks:
        parser.error("supply a workbook path, --material NAME or --all")
    if args.constants is not None and len(workbooks) != 1 and args.design_material is None:
        parser.error("--constants accepts exactly one workbook (constants belong to one material)")
    if args.calibration_id is not None:
        if args.constants is not None or len(workbooks) != 1:
            parser.error("--calibration-id needs exactly one workbook and no --constants")
        from multifat.cli.records import find_record

        archive = args.calibration_archive or args.results_dir / "calibrations"
        try:
            args.constants = find_record(
                archive, material_id(workbooks[0].stem), args.calibration_id
            )
        except FileNotFoundError as exc:
            parser.error(str(exc))
    design_workbook = None
    if args.design_material is not None:
        if args.calibration_id is not None:
            parser.error("--design-material takes the constants from the design workbook only")
        try:
            design_workbook = runmod.resolve_workbook(args.design_material, DATA_FOLDER)
        except WorkbookError as exc:
            parser.error(str(exc))
        if args.constants is not None and not is_published_constants(
            args.constants, material_id(design_workbook.stem)
        ):
            parser.error(
                "with --design-material, --constants must be a published-constants record of "
                "the design material"
            )
        workbooks = [w for w in workbooks if Path(w).resolve() != design_workbook.resolve()]
        if not workbooks:
            parser.error("--design-material needs at least one other workbook to predict")
    if args.results_dir.resolve() == Path.cwd().resolve():
        parser.error("--results-dir must be a subdirectory, not the working directory")
    selected = None if args.models == "all" else [m.strip() for m in args.models.split(",")]
    grouped: dict = {}
    for workbook in workbooks:
        grouped.setdefault(material_id(Path(workbook).stem), []).append(Path(workbook).name)
    duplicates = {material: files for material, files in grouped.items() if len(files) > 1}
    if duplicates:
        raise SystemExit(
            "duplicate material IDs: "
            + "; ".join(
                f"'{material}' in {', '.join(files)}" for material, files in duplicates.items()
            )
        )

    exit_code = 0
    run_map: dict = {}
    for workbook in workbooks:
        print(f"== {workbook}" + (f" predicted from {design_workbook}" if design_workbook else ""))
        run_id = material_id(workbook.stem)
        if design_workbook is not None:
            run_id = runmod.transfer_material_id(run_id, material_id(design_workbook.stem))
        try:
            outcomes = run_workbook(
                workbook,
                models=selected,
                all_variants=args.all_variants,
                results_dir=args.results_dir,
                jobs=args.jobs,
                calibrate=not args.no_calibrate,
                plots=args.plots,
                constants=args.constants,
                survival_probability=args.survival_probability,
                sn_model=args.sn_model,
                calibration_archive=args.calibration_archive,
                design_workbook=design_workbook,
            )
        except WorkbookError as exc:
            print(f"   workbook rejected: {exc}")
            run_map[run_id] = _rejection_record(workbook, exc)
            exit_code = 1
            continue
        run_map[run_id] = {
            "outcomes": outcomes,
            "calibration": _calibration_summary(args.results_dir, run_id),
        }
        for state in ("ok", "skipped", "failed"):
            subset = [
                outcome
                for outcome in outcomes
                if outcome.status == state or (state == "ok" and outcome.status.startswith("ok"))
            ]
            print(f"   {state}: {len(subset)}")
            if state != "ok":
                for outcome in subset:
                    print(f"     {outcome.model}: {outcome.reason.splitlines()[0][:160]}")
        exit_code |= any(outcome.status == "failed" for outcome in outcomes)

    if len(workbooks) > 1:  # a single-material run keeps the matrix of an earlier --all
        write_eligibility(
            run_map,
            args.results_dir / "Comparison",
            calibrate=not args.no_calibrate,
            models=selected,
        )
    if not args.no_compare:
        from multifat.cli import compare

        compare.main(
            [
                "--results_dir",
                str(args.results_dir),
                "--out_dir",
                str(args.results_dir / "Comparison"),
            ]
        )
    raise SystemExit(int(exit_code))


if __name__ == "__main__":
    main()
