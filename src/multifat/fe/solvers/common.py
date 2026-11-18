"""Shared validation and atomic CSV output for solver adapters."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import tempfile

from ..csvio import COLUMN_ALIASES, CsvHeader, resolve_columns
from ..csvwriter import FeCsvWriter
from ..model import sha256_of


@dataclass
class Frame:
    """One result state with original identity, time and canonical rows."""

    identity: tuple
    time: float | None
    rows: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    averaging: str | None = None


def paths_of(paths):
    if isinstance(paths, (str, Path)):
        paths = [paths]
    result = [Path(p) for p in paths]
    if not result:
        raise ValueError("at least one source file is required")
    return result


def number(text):
    value = float(text.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("solver numbers must be finite")
    return value


def label(text):
    if not text.isdigit() or int(text) < 1:
        raise ValueError(f"positive integer node/element label required: {text!r}")
    return int(text)


def column(name):
    """Resolve solver component names through the exchange-format alias table."""
    token = name.split(".")[-1].lower()
    if token in ("coord1", "coor1", "coord2", "coor2", "coord3", "coor3"):
        return {"1": "x", "2": "y", "3": "z"}[token[-1]]
    if token in ("element", "elementlabel"):
        return "element"
    if token in ("nodelabel", "node_number"):
        token = "node"
    try:
        return COLUMN_ALIASES[token]
    except KeyError:
        raise ValueError(f"unsupported solver column {name!r}") from None


def collapse(rows, *, element_nodal=False):
    """Collapse agreeing element-node values without averaging."""
    unique = {}
    for row in rows:
        row = {k: v for k, v in row.items() if k != "element"}
        key = tuple(row[k] for k in ("instance", "node_id", "material"))
        if key in unique:
            previous = unique[key]
            same = previous.keys() == row.keys() and all(
                (
                    math.isclose(v, previous[k], rel_tol=1e-10, abs_tol=1e-14)
                    if isinstance(v, (int, float))
                    else v == previous[k]
                )
                for k, v in row.items()
            )
            if not element_nodal or not same:
                raise ValueError(
                    f"non-unique node {key}; export Unique Nodal values within material, "
                    "or export separate element/material selections with labels"
                )
        else:
            unique[key] = row
    return list(unique.values())


def write_frames(frames, out_csv, *, solver, sources, metadata):
    """Validate a complete grid and publish CSV only after successful conversion."""
    metadata = {k: v for k, v in metadata.items() if v is not None}
    header = {
        "format": "multifat-fe-csv 1.0",
        "solver": solver,
        "solver_version": "unknown",
        **metadata,
    }
    # Validate declarations before opening an output, including missing units/frame.
    CsvHeader.from_mapping({"averaging": "none", **header})
    source_text = [f"{p} sha256={sha256_of(p)}" for p in sources]
    notes, columns, locations, identities = [], None, None, set()
    previous_time, has_time, count, averaging = None, None, 0, header.get("averaging")
    previous_identity = None
    out_csv = Path(out_csv)
    if out_csv.resolve() in {p.resolve() for p in sources}:
        raise ValueError("output must differ from every source")
    # Spool rows to disk: memory is bounded by a frame and its location grid.
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as spool:
        for frame in frames:
            if frame.identity in identities:
                raise ValueError("result frames must ascend without repetition")
            if identities:
                last = previous_identity
                if (
                    isinstance(frame.identity[0], int)
                    and isinstance(last[0], int)
                    and frame.identity < last
                ):
                    raise ValueError("numeric solver steps must ascend")
                if frame.identity[0] == last[0] and frame.identity[1:] <= last[1:]:
                    raise ValueError("result frames must ascend")
                if frame.identity[0] != last[0] and any(
                    i[0] == frame.identity[0] for i in identities
                ):
                    raise ValueError("solver steps must ascend without revisiting a step")
            identities.add(frame.identity)
            previous_identity = frame.identity
            count += 1
            if has_time is None:
                has_time = frame.time is not None
            if has_time != (frame.time is not None):
                raise ValueError("time must be present in every frame or none")
            if has_time:
                if not math.isfinite(frame.time) or (
                    previous_time is not None and frame.time <= previous_time
                ):
                    raise ValueError(
                        "time must strictly ascend; step-time resets need a total-time export"
                    )
                previous_time = frame.time
            if frame.averaging:
                if averaging is not None and averaging != frame.averaging:
                    raise ValueError(
                        "averaging argument or report metadata conflicts between frames"
                    )
                averaging = frame.averaging
            if not frame.rows:
                raise ValueError("no nodal stress rows in frame")
            grid = {}
            for row in frame.rows:
                row = dict(row, step=count)
                if has_time:
                    row["time"] = frame.time
                current_columns = tuple(row)
                if columns is None:
                    columns = current_columns
                    mapping = resolve_columns(columns)
                    mapping.check_header(
                        CsvHeader.from_mapping({**header, "averaging": averaging or "none"})
                    )
                if set(columns) != set(current_columns):
                    raise ValueError("column set changed between frames or nodes")
                key = tuple(row[k] for k in ("instance", "node_id", "material"))
                if key in grid:
                    raise ValueError("duplicate location in frame")
                xyz = tuple(row[k] for k in ("x", "y", "z"))
                if not all(
                    math.isfinite(float(v))
                    for k, v in row.items()
                    if k not in ("instance", "material")
                ):
                    raise ValueError("solver numbers must be finite")
                grid[key] = xyz
                spool.write(json.dumps(row) + "\n")
            if locations is None:
                locations = grid
            elif locations != grid:
                raise ValueError("location grid or undeformed coordinates changed between frames")
            notes.extend(frame.notes)
            notes.append(f"CSV step {count}: solver state {frame.identity}, time={frame.time}")
        if not count:
            raise ValueError("no stress frames found")
        if averaging is None:
            raise ValueError("averaging declaration is required")
        if averaging == "nodal_average_global" and (
            len({k[0] for k in locations}) > 1 or len({k[2] for k in locations}) > 1
        ):
            raise ValueError("global averaging across materials/instances is refused")
        header.update(
            averaging=averaging,
            n_nodes=len(locations),
            n_steps=count,
            source="; ".join(source_text + [n for n in notes if n.lower().startswith("odb:")]),
            notes="; ".join(dict.fromkeys(notes)),
        )
        validated = CsvHeader.from_mapping(header)
        handle = tempfile.NamedTemporaryFile(dir=out_csv.parent, suffix=".csv", delete=False)
        temporary = Path(handle.name)
        handle.close()
        try:
            spool.seek(0)
            raw = {k: v for k, v in asdict(validated).items() if v is not None and k != "extra"}
            raw.update(validated.extra)
            with FeCsvWriter(temporary, raw, columns=columns) as writer:
                for line in spool:
                    writer.write_row(json.loads(line))
            temporary.replace(out_csv)
        finally:
            temporary.unlink(missing_ok=True)
    return validated
