"""Assess the nodes of an FE import: node histories, unit-case superposition and a result table.

``node_histories`` turns each location of a ``multifat.fe`` history file into a ``TimeHistory``:
the load steps of a ``time_sample`` file are the samples (axis: the file's time, else the solver step
numbers as a cycle axis), and a ``unit_case`` file (linear elastic only) is superposed with a table
of load functions, ``S(t) = sum_k f_k(t) S_k`` and likewise for the strain. ``assess_nodes`` runs
``assess_history`` on every selected node with the setup (criterion, constants, Basquin curve) of the
node's material and returns one row per node, recording a refusal as a status instead of stopping.
Exact for proportional node histories; a non-proportional one is refused unless the caller passes
``allow_nonproportional_approximation=True``. This is software integration, not a validation of any
fatigue criterion.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Literal, Mapping, Sequence

import numpy as np

from multifat.fe.store import HistoryFile
from multifat.fe.tensors import elastic_strain

from .accumulate import SNCurve
from .assess import assess_history
from .eligibility import IneligibleHistoryError
from .signal import TimeHistory

STEP_AXIS_UNIT = "step"


@dataclass(frozen=True)
class MaterialSetup:
    """Criterion, its constants and the Basquin curve used for the nodes of one material.

    ``E`` and ``nu`` (MPa, dimensionless) are optional: when both are given, a node history of a
    linear-elastic file without strain gets the Hooke strain at assessment (recorded in the
    provenance), for criteria that read strain.
    """

    model: Any
    params: Mapping[str, float]
    curve: SNCurve
    E: float | None = None
    nu: float | None = None

    def __post_init__(self):
        if (self.E is None) != (self.nu is None):
            raise ValueError("give both E and nu, or neither")


@dataclass(frozen=True)
class NodeRef:
    """A location of a history file: its position, node label, instance and material."""

    position: int
    location_id: int
    instance: str
    material: str


def read_load_functions(
    path: str | Path, steps: Sequence[int]
) -> tuple[np.ndarray, np.ndarray, Literal["time", "cycle"]]:
    """Read the load functions of the unit cases of a ``unit_case`` file.

    The CSV has one axis column, ``time`` or ``index``, and one column ``f_<step>`` per unit case,
    matched to the file's step numbers by name (never by position); a missing or extra function
    column, a second axis column or a non-increasing axis is refused.

    Returns
    -------
    tuple
        The axis ``(Nt,)``, the functions ``F`` of shape ``(Nt, K)`` in the order of ``steps`` and
        the axis kind (``"time"`` or ``"cycle"`` for an index).
    """
    import pandas as pd

    path = Path(path)
    frame = pd.read_csv(path)
    columns = list(frame.columns)
    axes = [c for c in ("time", "index") if c in columns]
    if len(axes) != 1:
        raise ValueError(
            f"{path.name}: give exactly one axis column, 'time' or 'index' (found {axes or 'none'})"
        )
    wanted = [f"f_{int(s)}" for s in steps]
    others = [c for c in columns if c != axes[0]]
    missing = [c for c in wanted if c not in others]
    extra = [c for c in others if c not in wanted]
    if missing or extra:
        raise ValueError(
            f"{path.name}: the load-function columns must be exactly {wanted}; "
            f"missing {missing}, unexpected {extra}"
        )
    axis = frame[axes[0]].to_numpy(float)
    F = frame[wanted].to_numpy(float)
    if not (np.all(np.isfinite(axis)) and np.all(np.isfinite(F))):
        raise ValueError(f"{path.name}: the axis and the load functions must be finite")
    if len(axis) < 2 or np.any(np.diff(axis) <= 0.0):
        raise ValueError(f"{path.name}: the axis must have two or more strictly increasing values")
    return axis, F, "time" if axes[0] == "time" else "cycle"


def node_histories(
    path: str | Path,
    *,
    positions: Sequence[int] | None = None,
    load_functions: str | Path | tuple | None = None,
    chunk_locations: int = 4096,
) -> Iterator[tuple[NodeRef, TimeHistory]]:
    """Yield ``(NodeRef, TimeHistory)`` for the locations of a history file, chunk by chunk.

    Parameters
    ----------
    path : str or Path
        A ``multifat.fe`` history file (needs ``h5py``).
    positions : sequence of int, optional
        File positions to build, in this order; default all locations.
    load_functions : str, Path or tuple, optional
        Required for ``step_role: unit_case`` files (and refused for others): the CSV read by
        ``read_load_functions``, or its result ``(axis, F, axis_kind)``.
    chunk_locations : int
        Locations read from the file at a time.

    Notes
    -----
    The axis of a ``time_sample`` file is the load-case time if the file has one, else the solver
    step numbers (a cycle axis, unit ``"step"``). A unit-case superposition needs a linear-elastic
    file. A file with a single load step cannot give a history and is refused.
    """
    path = Path(path)
    with HistoryFile(path) as f:
        meta = f.meta
        unit_case = meta.step_role == "unit_case"
        steps = f.load_cases.steps
        functions = None
        provenance_extra: dict[str, Any] = {}
        if unit_case:
            if meta.analysis_type != "linear_elastic":
                raise ValueError(
                    f"{path.name}: unit-case superposition needs a linear-elastic analysis, the "
                    f"file declares {meta.analysis_type!r}"
                )
            if load_functions is None:
                raise ValueError(f"{path.name} holds unit cases: pass load_functions")
            if steps is None:
                raise ValueError(f"{path.name}: a unit-case file needs its solver step numbers")
            if isinstance(load_functions, tuple):
                functions = load_functions
                provenance_extra["load_functions"] = {"given": "arrays"}
            else:
                functions = read_load_functions(load_functions, steps)
                provenance_extra["load_functions"] = {"file": Path(load_functions).name}
            axis, F, axis_kind = functions
            if np.asarray(F).shape[1:] != (len(steps),):
                raise ValueError(
                    f"load functions have {np.asarray(F).shape[1:]} columns for {len(steps)} unit cases"
                )
            time_unit = meta.source_units.get("time") if axis_kind == "time" else STEP_AXIS_UNIT
        else:
            if load_functions is not None:
                raise ValueError(
                    f"{path.name} is a time-sample file: load functions apply only to unit cases"
                )
            if f.n_load_cases < 2:
                raise ValueError(
                    f"{path.name}: one load step is not a history (two or more needed)"
                )
            if f.load_cases.time is not None:
                axis, axis_kind = f.load_cases.time, "time"
                time_unit = meta.source_units.get("time")
            else:
                axis = (
                    steps.astype(float)
                    if steps is not None
                    else np.arange(f.n_load_cases, dtype=float)
                )
                axis_kind, time_unit = "cycle", STEP_AXIS_UNIT
        base = {
            "source": "fe_node",
            "file": path.name,
            "source_files": [{"name": s.name, "sha256": s.sha256} for s in meta.source_files],
            "solver": meta.solver,
            "analysis_type": meta.analysis_type,
            "step_role": meta.step_role,
            "cycle": meta.cycle,
            "strain_kind": meta.strain_kind,
            **provenance_extra,
        }
        order = np.arange(f.n_locations) if positions is None else np.asarray(positions, dtype=int)
        for start in range(0, len(order), max(1, chunk_locations)):
            chosen = order[start : start + max(1, chunk_locations)]
            block = f.read(chosen)
            for i, position in enumerate(chosen):
                S = block.S_hist[i]
                E = None if block.E_hist is None else block.E_hist[i]
                if functions is not None:
                    S = np.einsum("tk,kij->tij", F, S)
                    E = None if E is None else np.einsum("tk,kij->tij", F, E)
                ref = NodeRef(
                    int(position),
                    int(block.location_id[i]),
                    "" if block.instance is None else str(block.instance[i]),
                    "" if block.material is None else str(block.material[i]),
                )
                provenance = {
                    **base,
                    "position": ref.position,
                    "location_id": ref.location_id,
                    "instance": ref.instance,
                    "material": ref.material,
                    "material_id": ref.material,
                }
                yield ref, TimeHistory.from_arrays(
                    axis, S, E, stress_unit="MPa", axis_kind=axis_kind, time_unit=time_unit,
                    provenance=provenance,
                )  # fmt: skip


_EMPTY = dict(
    damage_per_block=None,
    blocks_to_failure=None,
    cycles_to_failure=None,
    n_events=None,
    approximation=None,
)


def assess_nodes(
    path: str | Path,
    materials: MaterialSetup | Mapping[str, MaterialSetup],
    *,
    R_list,
    positions: Sequence[int] | None = None,
    load_functions: str | Path | tuple | None = None,
    variant: str = "",
    on_refusal: Literal["record", "raise"] = "record",
    chunk_locations: int = 4096,
    **assess_options: Any,
) -> list[dict]:
    """Damage and life of every selected node of a history file, one row per node.

    Parameters
    ----------
    path : str or Path
        A ``multifat.fe`` history file (for example from ``import_fe_csv``).
    materials : MaterialSetup or mapping
        One setup for every node, or ``material label -> MaterialSetup`` (the label is the file's
        ``material`` value, the empty string when the file has none).
    R_list : ndarray
        Candidate plane frames, passed to ``assess_history``.
    positions, load_functions, chunk_locations
        See ``node_histories``.
    variant : {"", "_ext"}
        Damage-parameter variant of the criterion.
    on_refusal : {"record", "raise"}
        ``"record"``: a node that is refused (ineligible history, no setup for its material) gets
        a ``status`` and the run continues; ``"raise"``: the first refusal propagates.
    **assess_options
        Passed to ``assess_history`` (``counting``, ``channel``, ``residue``, ``endurance``,
        ``allow_nonproportional_approximation``, ...). Unexpected errors are never caught.

    Returns
    -------
    list of dict
        Keys ``position``, ``location_id``, ``instance``, ``material``, ``status`` (``"ok"``,
        ``"ineligible: <reason>"`` or ``"no_material_setup"``), ``damage_per_block``,
        ``blocks_to_failure``, ``cycles_to_failure``, ``n_events``, ``approximation``, ``notes``.
        The numeric entries are ``None`` for a refused node.
    """
    if on_refusal not in ("record", "raise"):
        raise ValueError(f"on_refusal must be 'record' or 'raise', got {on_refusal!r}")
    rows = []
    for ref, history in node_histories(
        path, positions=positions, load_functions=load_functions, chunk_locations=chunk_locations
    ):
        row = {
            "position": ref.position,
            "location_id": ref.location_id,
            "instance": ref.instance,
            "material": ref.material,
        }
        setup = materials if isinstance(materials, MaterialSetup) else materials.get(ref.material)
        if setup is None:
            if on_refusal == "raise":
                raise ValueError(f"no MaterialSetup for material {ref.material!r}")
            rows.append({**row, "status": "no_material_setup", **_EMPTY, "notes": ""})
            continue
        if (
            history.E_hist is None
            and setup.E is not None
            and history.provenance.get("analysis_type") == "linear_elastic"
        ):
            derived = elastic_strain(history.S_hist, setup.E, setup.nu)
            provenance = {
                **history.provenance,
                "strain_kind": "derived_linear_elastic_at_assessment",
            }
            history = dataclasses.replace(history, E_hist=derived, provenance=provenance)
        try:
            result = assess_history(
                history, setup.model, setup.params, setup.curve, R_list=R_list, variant=variant,
                **assess_options,
            )  # fmt: skip
        except IneligibleHistoryError as exc:
            if on_refusal == "raise":
                raise
            rows.append({**row, "status": f"ineligible: {exc}", **_EMPTY, "notes": ""})
            continue
        rows.append(
            {
                **row,
                "status": "ok",
                "damage_per_block": result.damage_per_block,
                "blocks_to_failure": result.blocks_to_failure,
                "cycles_to_failure": result.cycles_to_failure,
                "n_events": len(result.events),
                "approximation": result.approximation,
                "notes": "; ".join(result.notes),
            }
        )
    return rows
