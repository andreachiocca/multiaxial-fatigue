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
import platform
import subprocess
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from multifat.cli import run as runmod
import multifat.imp.experimental_data as exp
from multifat.imp.schema import WorkbookError, material_id, workbook_files
from multifat.models.errors import MissingParamsError

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
    mode: str = "life"  # "fatigue_index" on a fatigue-limit dataset


# Optional constants that a calibrated value may fill; required constants are always filled.
# Fitted Basquin/strain-life constants are never handed to models that read them only
# optionally (GSA and LIU2021 would silently abandon their documented fallbacks).
CALIBRATED_OPTIONAL = {"L_LI", "k_MGSE", "a_MSWT"}


def _every_model() -> list[str]:
    from multifat.models.registry import _available_cp_methods, _available_direct_models

    return list(_available_cp_methods()) + list(_available_direct_models())


def _evaluate_base(job: tuple) -> list[Outcome]:
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
    """Name the lost constants on the models that ran on their default or were skipped for them."""
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


def _evaluate_base_core(job: tuple) -> list[Outcome]:
    workbook, base_name, dest, overrides, design_ids, why, make_plots, calibration, suffixes = job
    start = time.perf_counter()
    adopted: list[str] = []
    try:
        from multifat.models.registry import get_models, stress_only_models

        flagged: dict = {}
        model = get_models([base_name])[0]
        suppliable = _constants_calibration_can_supply(model)
        data, material = runmod.load_material(workbook)
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
            data, material, base_name, params, extra_design_ids=design_points
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
    except Exception as exc:  # carry on with the remaining models
        text = f"{type(exc).__name__}: {exc} | {traceback.format_exc(limit=-2).strip()}"
        return [
            Outcome(base_name + suffix, "failed", text, 0, time.perf_counter() - start)
            for suffix in suffixes
        ]
    results = []
    fitted_label = ",".join(sorted(adopted))
    for suffix in suffixes:
        label = base_name + suffix
        if suffix == "_ext" and not any("_ext" in v for v in evaluation.values):
            results.append(Outcome(label, "skipped", "model has no _ext variant"))
            continue
        try:
            runmod.check_fit(base_name, suffix, evaluation.cal)
        except ValueError as exc:
            results.append(
                Outcome(
                    label,
                    "skipped",
                    f"design-set fit unusable: {exc}",
                    calibrated=fitted_label,
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
                )
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


def _constants_calibration_can_supply(model) -> set[str]:
    return set(model.required_params()) | (
        set(getattr(model, "optional_params", {})) & CALIBRATED_OPTIONAL
    )


def _absent_suppliable_constants(
    workbook: Path, names: list[str], data: dict | None = None
) -> set[str]:
    """Constants a selected model needs, calibration can supply, and the workbook does not have."""
    from multifat.cli.calibrate import CALIBRATED
    from multifat.eval.material_fits import BASQUIN_KEYS, PLASTIC_KEYS
    from multifat.models.registry import get_models, stress_only_models

    if data is None:
        data, material = runmod.load_material(workbook)
    else:
        material = material_id(workbook.stem)
    params = runmod.build_params(runmod.exp.get_material_properties(data, material))
    supplied = set(CALIBRATED.values()) | set(BASQUIN_KEYS) | set(PLASTIC_KEYS)
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
    them in ``params``. Returns (constants, description of the source for the new sidecar).
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
        recorded_failures = (record.get("counts") or {}).get("failure")
        if (
            n_failures is not None
            and recorded_failures is not None
            and recorded_failures != n_failures
        ):
            origin["warning"] = (
                f"record made on {recorded_failures} failures, the workbook has {n_failures}: "
                "its constants were not calibrated on these data"
            )
        if "a_MSWT" in constants:
            weight = constants["a_MSWT"]
            constants["a_MSWT"] = weight / (2.0 - weight)
            origin["note"] = f"a_MSWT {weight} converted to the printed weighting a/(2 - a)"
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
    algorithm = ((record.get("amplitude_settings") or {}).get("FIN") or {}).get("algorithm")
    if algorithm not in (None, "mce_harmonic"):
        # Production Findley evaluates MCE (A. Chiocca, 36ff0a1); an MCC k_FI would be misapplied.
        raise ValueError(
            f"{path}: k_FI was calibrated with the Findley amplitude '{algorithm}', but the "
            "pipeline evaluates Findley with mce_harmonic"
        )
    return constants, origin


def provenance(workbook: Path, data: dict | None = None) -> dict:
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
            "master_curve_stdnum": runmod.STDNUM,
            "flat_curve_k_max": runmod.K_SLOPE_MAX,
        },
    }


def _calibration_state(sidecar: dict, needed: set[str]) -> dict:
    """What an infeasible calibration lost, for the models that depended on the lost constants."""
    if sidecar.get("status") != "infeasible":
        return {}
    from multifat.cli.calibrate import CALIBRATED

    info = sidecar.get("infeasible") or {}
    lost = sorted(set(CALIBRATED.values()) & set(needed) - set(sidecar.get("params", {})))
    return {
        "status": "infeasible",
        "block": info.get("block"),
        "diagnostic": info.get("diagnostic") or "",
        "reason": sidecar.get("reason", ""),
        "lost": lost,
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
) -> list[Outcome]:
    results_dir.mkdir(parents=True, exist_ok=True)
    material = material_id(workbook.stem)

    report = runmod.validate_workbook(workbook)
    (results_dir / f"{material}_validation.txt").write_text(report.text(), encoding="utf-8")
    if not report.ok:
        raise WorkbookError(workbook.name, "", "", "; ".join(report.errors))

    from multifat.cli.method_selection import select_method_variants

    requested = select_method_variants(models, all_variants=all_variants)
    selected = list(dict.fromkeys(request.base_name for request in requested))
    data, _ = runmod.load_material(workbook)
    if not data[material]["cases"]:  # validated: a fatigue-limit dataset
        results, sidecar = _evaluate_fatigue_limit(workbook, material, requested, results_dir, jobs)
        sidecar["provenance"] = provenance(workbook, data)
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
    supplied, source = load_constants(constants, material, n_fail) if constants else ({}, None)
    needed = (
        (_absent_suppliable_constants(workbook, selected, data) - set(supplied))
        if calibrate
        else set()
    )
    sidecar: dict = {"params": {}}
    if needed:
        from multifat.cli.calibrate import CALIBRATED, calibrate_workbook

        sidecar = calibrate_workbook(
            workbook,
            planes=bool(needed & set(CALIBRATED.values())),
            constants=bool(needed - set(CALIBRATED.values())),
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
        sidecar["constants_source"] = source
    sidecar["provenance"] = provenance(workbook, data)
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
        (
            str(workbook),
            base_name,
            str(results_dir),
            overrides,
            design_ids,
            why,
            plots,
            calibration,
            tuple(suffixes),
        )
        for base_name, suffixes in suffixes_by_base.items()
    ]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as pool:
            results = [item for chunk in pool.map(_evaluate_base, queue) for item in chunk]
    else:
        results = [item for job in queue for item in _evaluate_base(job)]

    _dump_run_report(results_dir, material, results)
    return results


def _dump_run_report(dest: Path, material: str, results: list[Outcome]) -> None:
    sep = "\t"
    head = ["model", "status", "reason", "calibrated", "n_cases", "seconds", "notes"]
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
            limits = list(block["runouts"].values())
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
    from multifat.eval.material_fits import BASQUIN_KEYS, PLASTIC_KEYS

    return set(CALIBRATED.values()) | set(BASQUIN_KEYS) | set(PLASTIC_KEYS)


def _eligibility_status(outcome: Outcome, calibrate: bool) -> str:
    if outcome.status.startswith("ok") and outcome.mode == "fatigue_index":
        return "fatigue_index"
    if outcome.status.startswith("ok"):
        return "available" if outcome.n_cases else "prediction_only"
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
        "--constants",
        type=Path,
        default=None,
        help="constants from a Calibrations/ record or a calibration sidecar (one workbook)",
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
    if args.constants is not None and len(workbooks) != 1:
        parser.error("--constants accepts exactly one workbook (constants belong to one material)")
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
        print(f"== {workbook}")
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
            )
        except WorkbookError as exc:
            print(f"   workbook rejected: {exc}")
            run_map[material_id(workbook.stem)] = _rejection_record(workbook, exc)
            exit_code = 1
            continue
        run_map[material_id(workbook.stem)] = {
            "outcomes": outcomes,
            "calibration": _calibration_summary(args.results_dir, material_id(workbook.stem)),
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
