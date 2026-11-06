"""Assessment of a general load history: events, damage parameters per event, accumulation, record.

The chain is: check eligibility; take a scalar channel; count it (rainflow) or take the longest
chord; build events (pairs of samples with a count); evaluate the criterion on each distinct event
through the two-state bridge of ``multifat.fe``; turn each damage parameter into a life with the
caller's Basquin curve; sum the damage of one block. The block is assumed to repeat until failure,
so the default residue policy is ``"rearrange"`` (every cycle closes, exact for strict repetition);
``"half"`` is the standard rule for a single, non-repeating history and under-counts the cycle that
spans consecutive blocks, which the result notes.

Exact for proportional histories (a straight stress path: the counted cycles are two-state paths);
a non-proportional history is refused unless ``allow_nonproportional_approximation=True``, which is
recorded. For criteria that read strain the strain path is checked too. ``event_damage_parameters``
gives the damage parameter of each event without a curve (also for native-life methods, which cannot
accumulate through a Basquin curve). The result carries everything needed to reproduce it: the
settings, the curve and the SHA-256 of the sampled history (``TimeHistory.record`` stores the arrays).
Nothing here validates a fatigue criterion against experiments.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Mapping

import numpy as np

from multifat.fe.bridge import to_case_inputs
from multifat.fe.ranges import ExplicitPair
from multifat.models.errors import LoadingScopeError
from multifat.models.registry import stress_only_models

from .accumulate import Endurance, SNCurve, miner
from .channels import resolve_channel, weighted_path
from .eligibility import (
    IneligibleHistoryError,
    check_accumulation_eligibility,
    check_event_eligibility,
    check_history_kind,
)
from .events import (
    check_event_count,
    events_from_cycles,
    events_from_longest_chord,
    merge_identical,
)
from .proportionality import Proportionality, proportionality
from .rainflow import rainflow
from .signal import TimeHistory

HALF_RESIDUE_NOTE = (
    "residue='half' counts the cycle that spans consecutive blocks as half cycles; for a block that "
    "repeats until failure use residue='rearrange' (the damage per block is then typically larger)"
)


def _jsonable(value: Any) -> Any:
    """JSON-safe copy: numpy to Python, non-finite floats to the strings ``"inf"``, ``"-inf"``, ``"nan"``."""
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return "nan" if np.isnan(value) else ("inf" if value > 0 else "-inf")
    return value


@dataclass(frozen=True)
class HistoryAssessment:
    """Result of ``assess_history``; ``blocks_to_failure`` and the lives are infinite at zero damage."""

    method: str
    variant: str
    curve: dict
    settings: dict
    history_sha256: str
    history_provenance: dict
    events: list[dict]
    damage_per_block: float
    blocks_to_failure: float
    cycles_per_block: float
    cycles_to_failure: float
    approximation: str | None
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON-safe dictionary (infinite values are the string ``"inf"``)."""
        return _jsonable({name: getattr(self, name) for name in self.__dataclass_fields__})

    def to_json(self, path: str | Path) -> Path:
        """Write ``to_dict`` as JSON; returns the path."""
        path = Path(path)
        path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True, allow_nan=False), encoding="utf-8"
        )
        return path


def _verdicts(history: TimeHistory, name: str, tol: float, allow: bool) -> tuple[dict, str | None]:
    """Proportionality of the stress path, and of the strain path for criteria that read strain."""
    record, notes = {}, []
    stress = proportionality(history.S_hist, tol=tol)
    record["stress"] = {
        "rank_ratio": stress.rank_ratio,
        "tol": stress.tol,
        "proportional": stress.proportional,
    }
    note = check_history_kind(stress, allow_approximation=allow, what="stress")
    notes.append(note) if note else None
    if history.E_hist is not None and name not in stress_only_models():
        strain = proportionality(history.E_hist, tol=tol)
        record["strain"] = {
            "rank_ratio": strain.rank_ratio,
            "tol": strain.tol,
            "proportional": strain.proportional,
        }
        note = check_history_kind(strain, allow_approximation=allow, what="strain")
        notes.append(note) if note else None
        # Rank one in each tensor space does not guarantee a shared load function.
        # Normalize each path separately so stress units cannot hide cycling strain.
        paths = [weighted_path(values - values[:1]) for values in (history.S_hist, history.E_hist)]
        amplitudes = [float(np.max(np.abs(p), initial=0.0)) for p in paths]
        normalised = [p / a / np.linalg.norm(p / a) if a else p for p, a in zip(paths, amplitudes)]
        singular = np.linalg.svd(np.concatenate(normalised, axis=1), compute_uv=False)
        rank_ratio = 0.0 if singular[0] == 0.0 else float(singular[1] / singular[0])
        # Counting uses a stress channel: a constant stress cannot parameterize varying strain.
        if amplitudes[0] == 0.0 and amplitudes[1] != 0.0:
            rank_ratio = 1.0
        joint = Proportionality(rank_ratio, rank_ratio <= tol, float(tol))
        record["joint"] = {
            "rank_ratio": joint.rank_ratio,
            "tol": joint.tol,
            "proportional": joint.proportional,
        }
        note = check_history_kind(joint, allow_approximation=allow, what="joint stress/strain")
        notes.append(note) if note else None
    record.update(
        record["stress"]
    )  # flat keys of the stress verdict, as recorded since the first version
    return record, ("; ".join(notes) if notes else None)


def _evaluate(model, params, R_list, block, meta, pair, variant: str, name: str) -> float:
    try:
        (case,) = to_case_inputs(block, ExplicitPair(*pair), meta=meta, method=model)
        result = model.evaluate_case(
            S0=case.S0, E0=case.E0, S1=case.S1, E1=case.E1, params=params, R_list=R_list,
            harmonics=case.harmonics,
        )  # fmt: skip
    except LoadingScopeError as exc:
        raise IneligibleHistoryError(
            f"{name}: the event between samples {tuple(pair)} lies outside the loading scope of the "
            f"criterion: {exc}"
        ) from exc
    if variant not in result.values:
        raise IneligibleHistoryError(
            f"{name} has no variant {variant!r} (available: {sorted(result.values)})"
        )
    value = float(result.values[variant])
    if not np.isfinite(value):
        raise ValueError(
            f"{name}: non-finite damage parameter for the event between samples {tuple(pair)}"
        )
    return value


def _table(
    history, model, params, *, R_list, variant, counting, channel, gate, residue, chord_metric, chord_on,
    proportionality_tol, allow_nonproportional_approximation, max_events,
):  # fmt: skip
    """Events with their damage parameters, and the settings and notes common to both entry points."""
    name = getattr(model, "name", type(model).__name__)
    check_event_eligibility(model, has_strain=history.E_hist is not None)
    verdict, approximation = _verdicts(
        history, name, proportionality_tol, allow_nonproportional_approximation
    )

    notes = []
    if counting == "rainflow":
        scalar = resolve_channel(channel, history)
        events = events_from_cycles(history, rainflow(scalar, gate=gate, residue=residue))
        channel_name = (
            channel if isinstance(channel, str) else getattr(channel, "__name__", repr(channel))
        )
        chord = None
        if residue == "half":
            notes.append(HALF_RESIDUE_NOTE)
    elif counting == "longest_chord":
        events = events_from_longest_chord(history, metric=chord_metric, on=chord_on)
        channel_name, residue = "not applicable", "not applicable"
        chord = {"metric": chord_metric, "on": chord_on}
    else:
        raise ValueError(f"counting must be 'rainflow' or 'longest_chord', got {counting!r}")
    events = merge_identical(history, events)
    check_event_count(events, max_events)

    block, meta = history.to_block(), history.meta()
    dp = np.array(
        [_evaluate(model, params, R_list, block, meta, e.pair, variant, name) for e in events]
    )
    if not events:  # nothing to evaluate: still validate the variant on the end-to-end pair
        try:
            _evaluate(model, params, R_list, block, meta, (0, history.n_samples - 1), variant, name)
        except IneligibleHistoryError as exc:
            if "variant" in str(exc):
                raise

    grid = (
        None
        if R_list is None
        else np.ascontiguousarray(np.asarray(R_list, dtype=float), dtype="<f8")
    )
    settings = {
        "counting": counting,
        "channel": channel_name,
        "gate": gate,
        "residue": residue,
        "chord": chord,
        "variant": variant,
        "plane_grid": {
            "n_planes": None if grid is None else int(len(grid)),
            "sha256": None if grid is None else hashlib.sha256(grid.tobytes()).hexdigest(),
        },
        "proportionality": verdict,
        "allow_nonproportional_approximation": allow_nonproportional_approximation,
        "max_events": max_events,
    }
    return name, events, dp, settings, approximation, notes


def event_damage_parameters(
    history: TimeHistory,
    model,
    params: Mapping[str, float],
    *,
    R_list,
    variant: str = "",
    counting: Literal["rainflow", "longest_chord"] = "rainflow",
    channel: Any = "dominant",
    gate: float = 0.0,
    residue: Literal["half", "rearrange"] = "rearrange",
    chord_metric: str = "deviatoric_norm",
    chord_on: str = "stress",
    proportionality_tol: float = 1e-6,
    allow_nonproportional_approximation: bool = False,
    max_events: int = 2000,
) -> list[dict]:
    """Damage parameter of each distinct event of a history, without a life curve.

    Available for every method that can evaluate the history, including native-life methods (which
    cannot accumulate through a Basquin curve of their damage parameter). Arguments as in
    ``assess_history``. Returns one dict per event: ``pair``, ``time``, ``count``, ``span``, ``dp``.
    """
    _, events, dp, _, _, _ = _table(
        history, model, params, R_list=R_list, variant=variant, counting=counting, channel=channel,
        gate=gate, residue=residue, chord_metric=chord_metric, chord_on=chord_on,
        proportionality_tol=proportionality_tol,
        allow_nonproportional_approximation=allow_nonproportional_approximation, max_events=max_events,
    )  # fmt: skip
    return [
        {
            "pair": list(e.pair),
            "time": [float(history.time[e.pair[0]]), float(history.time[e.pair[1]])],
            "count": e.count,
            "span": e.span,
            "dp": float(dp[k]),
        }
        for k, e in enumerate(events)
    ]


def assess_history(
    history: TimeHistory,
    model,
    params: Mapping[str, float],
    curve: SNCurve,
    *,
    R_list,
    variant: str = "",
    counting: Literal["rainflow", "longest_chord"] = "rainflow",
    channel: Any = "dominant",
    gate: float = 0.0,
    residue: Literal["half", "rearrange"] = "rearrange",
    chord_metric: str = "deviatoric_norm",
    chord_on: str = "stress",
    endurance: Endurance | None = None,
    proportionality_tol: float = 1e-6,
    allow_nonproportional_approximation: bool = False,
    max_events: int = 2000,
) -> HistoryAssessment:
    """Damage and life of a block of load history under one criterion.

    Parameters
    ----------
    history : TimeHistory
        One block, assumed to repeat until failure.
    model : FatigueModel
        A case-level model (``multifat.models.registry.get_models``).
    params : mapping
        Material constants of the criterion.
    curve : SNCurve
        Basquin curve of the same method and variant, calibrated on the same plane set as ``R_list``.
    R_list : ndarray
        Candidate plane frames, as in the existing pipeline (its size and hash are recorded).
    variant : {"", "_ext"}
        Which damage parameter of the model is used.
    counting : {"rainflow", "longest_chord"}
        Rainflow on a scalar channel, or the single longest-chord event.
    channel, gate
        Scalar channel (see ``channels.resolve_channel``) and hysteresis gate of the rainflow
        counting (ignored by ``longest_chord``).
    residue : {"rearrange", "half"}
        ``"rearrange"`` (default): the block repeats back to back and every cycle closes, exact for
        strict repetition. ``"half"``: the standard rule for a single history; it counts the cycle
        that spans consecutive blocks as halves and so under-counts a repeating block (noted).
    chord_metric, chord_on
        Distance and tensor of the longest chord (see ``multifat.fe.ranges.LongestChord``).
    endurance : Endurance, optional
        Treatment of events below a damage-parameter limit; default extrapolates the curve.
    proportionality_tol : float
        Tolerance of the proportionality verdict (second over first singular value of the path).
    allow_nonproportional_approximation : bool
        Accept a non-proportional history as an approximation; recorded in the result.
    max_events : int
        Guard on the number of distinct events.

    Raises
    ------
    IneligibleHistoryError
        With the reason: native-life or nonlocal method, missing strain, curve mismatch, a
        non-proportional stress (or strain) path, a loading outside the criterion's scope, a missing
        variant.
    """
    endurance = endurance or Endurance()
    check_accumulation_eligibility(
        model, curve, variant=variant, has_strain=history.E_hist is not None
    )
    name, events, dp, settings, approximation, notes = _table(
        history, model, params, R_list=R_list, variant=variant, counting=counting, channel=channel,
        gate=gate, residue=residue, chord_metric=chord_metric, chord_on=chord_on,
        proportionality_tol=proportionality_tol,
        allow_nonproportional_approximation=allow_nonproportional_approximation, max_events=max_events,
    )  # fmt: skip
    counts = np.array([e.count for e in events])
    summed = miner(dp, counts, curve, endurance)
    damage = summed.total
    cycles_per_block = float(counts.sum())
    blocks = 1.0 / damage if damage > 0.0 else float("inf")
    short = int(np.sum(np.isfinite(summed.n_f) & (summed.n_f < 1.0)))
    if short:
        notes.append(
            f"N_f < 1 for {short} event(s): the damage parameter exceeds the curve's A, outside the "
            "range the curve describes"
        )
    rows = [
        {
            "pair": list(e.pair),
            "time": [float(history.time[e.pair[0]]), float(history.time[e.pair[1]])],
            "count": e.count,
            "span": e.span,
            "dp": float(dp[k]),
            "n_f": float(summed.n_f[k]),
            "damage": float(summed.damage[k]),
        }
        for k, e in enumerate(events)
    ]
    settings["endurance"] = {"mode": endurance.mode, "dp_limit": endurance.dp_limit}
    return HistoryAssessment(
        method=name,
        variant=variant,
        curve={
            "A": curve.A,
            "b": curve.b,
            "method": curve.method,
            "variant": curve.variant,
            "design_loading": curve.design_loading,
            "source": curve.source,
        },  # fmt: skip
        settings=settings,
        history_sha256=history.sha256(),
        history_provenance=dict(history.provenance),
        events=rows,
        damage_per_block=damage,
        blocks_to_failure=blocks,
        cycles_per_block=cycles_per_block,
        cycles_to_failure=cycles_per_block * blocks if damage > 0.0 else float("inf"),
        approximation=approximation,
        notes=notes,
    )
