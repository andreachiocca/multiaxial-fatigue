"""Convert Abaqus field output reports to the FE CSV exchange format.

Repository choices:
Only complete global Cartesian nodal tensor tables are supported. Ambiguous section points,
local frames and non-agreeing element-node values are refused rather than averaged.
"""

from __future__ import annotations

from pathlib import Path
import re

from .common import Frame, collapse, column, label, number, paths_of, write_frames


def _coordinates(inp):
    coordinates, reading = {}, False
    if inp is None:
        return coordinates
    with open(inp, encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                command = line.split(",")[0].upper()
                if command in (
                    "*PART",
                    "*INSTANCE",
                    "*SYSTEM",
                    "*INCLUDE",
                    "*NGEN",
                    "*NFILL",
                    "*NCOPY",
                    "*NMAP",
                ):
                    raise ValueError(
                        "INP coordinates require flat explicit global *NODE data; part/instance transforms and generated nodes are unsupported"
                    )
                reading = command == "*NODE"
                if reading and re.search(r"SYSTEM\s*=\s*[CS]", line, re.I):
                    raise ValueError("INP *NODE coordinates must be global Cartesian")
                continue
            if reading:
                cells = [v.strip() for v in line.split(",")]
                if len(cells) != 4:
                    raise ValueError("*NODE needs label and three undeformed coordinates")
                node = label(cells[0])
                if node in coordinates:
                    raise ValueError("duplicate INP node label")
                coordinates[node] = tuple(number(v) for v in cells[1:])
    return coordinates


def _frames(paths, coordinates, material, instance):
    time_basis, coordinate_instance = None, None
    for path in paths:
        frame, step, odb, names, position = None, None, None, None, None
        part, table, table_started = instance or "", False, False
        schema = None
        with open(path, encoding="utf-8-sig") as handle:
            for raw in handle:
                line = raw.strip()
                lower = line.lower()
                if lower.startswith(
                    ("field output report", "odb:", "step:", "frame:", "field output reported")
                ):
                    table = False
                if lower.startswith("odb:"):
                    odb = line
                if lower.startswith("step:"):
                    step = line.split(":", 1)[1].strip()
                if lower.startswith("frame:"):
                    if frame is not None:
                        frame.rows = collapse(frame.rows, element_nodal=position == "element")
                        yield frame
                    match = re.search(
                        r"Increment\s+(\d+).*?(Step Time|Total Time)\s*=\s*(\S+)", line, re.I
                    )
                    if match is None or step is None:
                        raise ValueError(
                            "frame requires Step, Increment and Step Time or Total Time"
                        )
                    if time_basis is not None and time_basis != match[2].lower():
                        raise ValueError("time basis changed between Step Time and Total Time")
                    time_basis = match[2].lower()
                    frame = Frame((step, int(match[1])), number(match[3]))
                    frame.notes.extend(v for v in (odb, line) if v)
                    names, schema, table, table_started, position = None, None, False, False, None
                if re.search(
                    r"\blocal\b|coordinate system\s*[:=]\s*[1-9]|\bsection point\b", lower
                ):
                    raise ValueError(
                        "local frames and section points are unsupported; export global Cartesian solid nodal tensors"
                    )
                if "material:" in lower:
                    stated = line.split(":", 1)[1].strip()
                    if stated != material:
                        raise ValueError("material argument conflicts with report material")
                if frame is None:
                    continue
                if "element nodal" in lower or "element nodes" in lower:
                    position = "element"
                elif "unique nodal" in lower:
                    position = "unique"
                if any(term in lower for term in ("averag", "computation algorithm", "loc ")):
                    frame.notes.append(line)
                    if "not averaged across region boundaries" in lower or "odb_regions" in lower:
                        declared = "nodal_average_within_material"
                    elif "averaged at nodes" in lower or "unique nodal" in lower:
                        declared = frame.averaging or "nodal_average_global"
                    else:
                        declared = None
                    if declared:
                        frame.averaging = declared
                if "unique nodal" in lower and frame.averaging is None:
                    frame.averaging = "nodal_average_global"
                if "for part:" in lower or "for instance:" in lower:
                    stated = line.split(":", 1)[1].strip()
                    if instance is not None and instance != stated:
                        raise ValueError("instance argument conflicts with report")
                    part = stated
                if re.match(r"(?:Element|Node)\s+(?:Label\s+)?(?:Node|S\.|E\.|EE\.|PE\.)", line):
                    tokens = re.sub(r"\b(Node|Element)\s+Label\b", r"\1", line).split()
                    names = [column(token) for token in tokens]
                    if len(names) != len(set(names)):
                        raise ValueError("duplicate solver columns")
                    if schema is not None and schema != names:
                        raise ValueError("column set changed within frame")
                    schema = names
                    table, table_started = True, False
                    continue
                if table and names:
                    cells = line.split()
                    if cells and re.fullmatch(r"\d+", cells[0]):
                        if len(cells) != len(names):
                            raise ValueError("ragged nodal table: column count differs")
                        row = {
                            name: (
                                label(value) if name in ("node_id", "element") else number(value)
                            )
                            for name, value in zip(names, cells)
                        }
                        if position is None:
                            raise ValueError(
                                "report must state Unique Nodal or Element Nodal position"
                            )
                        if not all(k in row for k in ("x", "y", "z")):
                            if coordinate_instance is not None and coordinate_instance != part:
                                raise ValueError(
                                    "a flat INP cannot supply coordinates for multiple instances; export COORD columns"
                                )
                            coordinate_instance = part
                            try:
                                row.update(zip(("x", "y", "z"), coordinates[row["node_id"]]))
                            except KeyError:
                                raise ValueError(
                                    "missing undeformed coordinates: provide a flat INP or COORD columns"
                                ) from None
                        elif (
                            row["node_id"] in coordinates
                            and tuple(row[k] for k in ("x", "y", "z"))
                            != coordinates[row["node_id"]]
                        ):
                            raise ValueError("report coordinates conflict with INP")
                        row.update(instance=part, material=material)
                        frame.rows.append(row)
                        table_started = True
                    elif table_started and line:
                        if lower.startswith(("minimum", "maximum", "total")):
                            table = False
                        elif not set(line) <= {"-", " "}:
                            raise ValueError("unexpected text inside nodal table")
        if frame is not None:
            frame.rows = collapse(frame.rows, element_nodal=position == "element")
            yield frame


def abaqus_report_to_csv(
    report_paths, out_csv, *, inp=None, material="", instance=None, **metadata
):
    """Convert complete Abaqus field reports, with explicit CSV units and conventions.

    Parameters
    ----------
    report_paths : path or sequence of paths
        Field Output Reports in increasing result-state order.
    out_csv : path
        Destination, replaced atomically on success.
    inp : path, optional
        Flat global *NODE coordinates if not present as COORD columns.
    material, instance : str, optional
        Material label and an instance assertion (otherwise parsed from the report).
    **metadata
        CSV header declarations: stress_unit, length_unit, analysis_type, step_role,
        coordinate_system and, where relevant, strain_convention, strain_unit, cycle, time_unit.

    Returns
    -------
    CsvHeader
        The validated header written to disk.
    """
    paths = paths_of(report_paths)
    sources = paths + ([Path(inp)] if inp is not None else [])
    return write_frames(
        _frames(paths, _coordinates(inp), material, instance),
        out_csv,
        solver="Abaqus",
        sources=sources,
        metadata=metadata,
    )
