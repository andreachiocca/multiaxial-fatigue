"""Observation identities, dataset fingerprints and matched failure subsets.

Ported from upstream ``fatigue/report/subsets.py`` and ``calibration_store.input_fingerprint``
(A. Chiocca, 904c619). A case is identified by its specimen and stress history; an observation
also by its cycle count, status and replicate occurrence. A dataset fingerprint hashes the
material properties and the failure histories only, so row order and run-outs do not change it.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

IDENTITY_COLUMNS = [
    "Test_status",
    "Specimen_id",
    "Case_id",
    "Observation_id",
    "Dataset_fingerprint",
]
KEYS = ["Material", "Dataset_fingerprint", "Observation_id"]
METRICS = {"life": "Error_ln", "damage": "Error_log10_dp"}


def json_value(value):
    """JSON-safe copy: arrays to lists, non-finite floats to None."""
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _sha(content) -> str:
    return hashlib.sha256(json.dumps(json_value(content), sort_keys=True).encode()).hexdigest()


def _history(case) -> dict:
    harm = case.get("harmonics") or {}
    return {
        "stress": {step: case[1][step]["stress"] for step in (1, 2)},
        "harmonics": {k: v for k, v in harm.items() if k.startswith("S_")},
    }


def dataset_fingerprint(block: dict, runouts: bool = False) -> str:
    """Hash of the material properties and the failure histories (lives and stress histories).

    With ``runouts`` the run-out histories are hashed too: a fatigue-limit dataset holds only
    run-outs (its limit states), so without them its fingerprint would track the properties only.
    The default keeps every recorded fingerprint unchanged.
    """
    props = {k: v for k, v in block["material"].items() if k != "extra"}
    props.update(block["material"].get("extra") or {})

    # Failures at a phase that is not a multiple of 90 deg hold trigonometric values whose last
    # bits differ between platforms (first in 18_Al7075_T651_Wang): such a dataset hashes its
    # failure histories rounded to 1e-6 MPa, like the run-outs. Every other recorded
    # fingerprint is unchanged.
    decimals = RUNOUT_DECIMALS if _has_oblique_phase(block["cases"]) else None

    def histories(cases) -> list[str]:
        def history(c):
            h = _history(c)
            return h if decimals is None else _rounded(h, decimals)

        return sorted(
            json.dumps(json_value({"N": c["meta"]["Nf_cycles"], **history(c)}), sort_keys=True)
            for c in cases.values()
        )

    content = {"properties": props, "failures": histories(block["cases"])}
    if runouts:
        # rounded to 1e-6 MPa: out-of-phase harmonics hold trigonometric values (cos 90 deg =
        # 6e-17) whose last bits differ between platforms' math libraries
        content["runouts"] = sorted(
            json.dumps(
                json_value({"N": c["meta"]["Nf_cycles"], **_rounded(_history(c), RUNOUT_DECIMALS)}),
                sort_keys=True,
            )
            for c in (block.get("runouts") or {}).values()
        )
    return _sha(content)


RUNOUT_DECIMALS = 6


def _has_oblique_phase(cases) -> bool:
    """True when a case has a load phase that is not a multiple of 90 deg."""
    for c in cases.values():
        for key, value in c.get("meta", {}).items():
            if "phase" not in key.lower():
                continue
            try:
                phase = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(phase) and phase % 90.0 != 0.0:
                return True
    return False


def _rounded(value, decimals: int):
    """Copy with every array rounded (and -0.0 turned into 0.0)."""
    if isinstance(value, dict):
        return {k: _rounded(v, decimals) for k, v in value.items()}
    if isinstance(value, (np.ndarray, list, tuple)):
        return np.round(np.asarray(value, dtype=float), decimals) + 0.0
    return value


def observation_index(block: dict) -> dict[tuple[str, int], dict]:
    """Identity columns of every failure (``("case", id)``) and run-out (``("runout", id)``)."""
    fingerprint = dataset_fingerprint(block)
    occurrences: dict[str, int] = defaultdict(int)
    out = {}
    groups = (("case", block["cases"]), ("runout", block.get("runouts", {})))
    for kind, cases in groups:
        for cid, case in cases.items():
            meta = case["meta"]
            specimen = str(meta.get("specimen_id", ""))
            case_key = _sha({"specimen": specimen, **_history(case)})
            status = meta.get("test_status", "failure" if kind == "case" else "runout")
            observation = f"{case_key}/{meta['Nf_cycles']}/{status}"
            occurrences[observation] += 1
            obs_id = hashlib.sha256(
                f"{observation}/{occurrences[observation]}".encode()
            ).hexdigest()
            out[(kind, int(cid))] = {
                "Test_status": status,
                "Specimen_id": specimen,
                "Case_id": case_key,
                "Observation_id": obs_id,
                "Dataset_fingerprint": fingerprint,
            }
    return out


def eligible_failures(df: pd.DataFrame, metric: str = "life") -> pd.DataFrame:
    values = pd.to_numeric(df[METRICS[metric]], errors="coerce")
    nf = pd.to_numeric(df["Nf_exp"], errors="coerce")
    mask = (df["Test_status"] == "failure") & np.isfinite(values) & np.isfinite(nf) & (nf > 0)
    if metric == "life":
        npred = pd.to_numeric(df["Nf_expected"], errors="coerce")
        mask &= np.isfinite(npred) & (npred > 0)
    return df.loc[mask].copy()


def common_subset(df: pd.DataFrame, methods, metric: str = "life") -> pd.DataFrame:
    """Rows of `methods` on the failure observations that every method assesses with `metric`."""
    methods = list(dict.fromkeys(methods))
    if len(methods) < 2:
        raise ValueError("Select at least two methods for a common comparison")
    if any(k not in df for k in KEYS):
        raise ValueError("Common comparisons need results written with observation identities")
    selected = eligible_failures(df[df["Method"].isin(methods)], metric)
    if selected[KEYS].isna().any().any() or (selected[KEYS].astype(str) == "").any().any():
        raise ValueError("Unresolved observation identity; regenerate the results")
    if selected.duplicated(["Method", *KEYS]).any():
        raise ValueError("Duplicate method-observation results; keep one result per observation")
    sets = [set(map(tuple, selected.loc[selected.Method == m, KEYS].to_numpy())) for m in methods]
    common = set.intersection(*sets)
    subset = selected.loc[[tuple(r) in common for r in selected[KEYS].to_numpy()]].copy()
    # The same identity must not hide different lives or design/prediction roles.
    for column in ["Nf_exp", "Flag", "Case_id"]:
        if column in subset and (subset.groupby(KEYS)[column].nunique() > 1).any():
            hint = (
                ": the methods use different design sets; compare predictions only with "
                "--flag_filter other"
                if column == "Flag"
                else ""
            )
            raise ValueError(f"Conflicting {column} for matched observations{hint}")
    return subset


def write_subset_reports(df: pd.DataFrame, directory, methods=None):
    """Per-method availability and pairwise/group overlaps, for life and damage metrics."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    methods = list(methods or sorted(df["Method"].unique()))
    counts, overlaps = [], []
    for metric in METRICS:
        eligible = eligible_failures(df, metric)
        for method in methods:
            part = eligible[eligible.Method == method]
            counts.append(
                {
                    "method": method,
                    "metric": metric,
                    "materials": part.Material.nunique(),
                    "observations": len(part),
                    "cases": part.Case_id.nunique() if "Case_id" in part else None,
                }
            )
        groups = list(combinations(methods, 2))
        if len(methods) > 2:
            groups.append(tuple(methods))
        for group in groups:
            try:
                common = common_subset(df, group, metric)
                first = common[common.Method == group[0]]
                overlaps.append(
                    {
                        "methods": group,
                        "metric": metric,
                        "materials": sorted(first.Material.unique()),
                        "n_observations": len(first),
                        "n_cases": first.Case_id.nunique(),
                        "observation_keys": first[KEYS].to_dict("records"),
                    }
                )
            except ValueError as exc:
                overlaps.append(
                    {"methods": group, "metric": metric, "unavailable_reason": str(exc)}
                )
    pd.DataFrame(counts).to_csv(directory / "individual_availability.csv", index=False)
    text = json.dumps(json_value(overlaps), indent=2, allow_nan=False) + "\n"
    (directory / "common_subsets.json").write_text(text, encoding="utf-8")
    return counts, overlaps
