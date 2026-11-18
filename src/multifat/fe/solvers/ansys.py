"""Convert complete Ansys PRNSOL/NLIST and Workbench tensor tables to FE CSV.

Repository choices:
Only global Cartesian nodal output with explicit averaging is accepted. Shell section points,
element tables and scalar Workbench results cannot supply a unique complete nodal tensor.
"""

from __future__ import annotations

import re

from .common import Frame, column, label, number, paths_of, write_frames


def _nodes(path):
    result, reading = {}, False
    if path is None:
        return result
    with open(path, encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            match = re.search(r"DSYS\s*=\s*(\d+)", line, re.I)
            if match and int(match[1]) != 0:
                raise ValueError("NLIST DSYS must be 0 (global Cartesian)")
            if re.match(r"NODE\s+X\s+Y\s+Z\b", line, re.I):
                reading = True
                continue
            cells = line.split()
            if reading and cells and cells[0].isdigit():
                if len(cells) not in (4, 7):
                    raise ValueError("ragged NLIST coordinates")
                node = label(cells[0])
                if node in result:
                    raise ValueError("duplicate NLIST node")
                result[node] = tuple(number(v) for v in cells[1:4])
            elif reading and line and not set(line) <= {"-", " "}:
                reading = False
    return result


def _finish(frame, coordinates, material, instance):
    joined = {}
    for row in frame.rows:
        node = row["node_id"]
        if node not in joined:
            joined[node] = dict(node_id=node, instance=instance, material=material)
        overlap = set(joined[node]) & (set(row) - {"node_id"})
        if overlap:
            raise ValueError(
                f"duplicate Ansys node/component: {node}, {overlap}; frames must ascend and component tables must be disjoint"
            )
        joined[node].update(row)
    for node, row in joined.items():
        if not all(k in row for k in ("x", "y", "z")):
            if node not in coordinates:
                raise ValueError(
                    "missing undeformed coordinates: provide global NLIST or X/Y/Z columns"
                )
            row.update(zip(("x", "y", "z"), coordinates[node]))
    frame.rows = list(joined.values())
    return frame


def _workbench(path, coordinates, material, instance, metadata, step, time):
    with open(path, encoding="utf-8-sig") as handle:
        first = handle.readline().strip()
        delimiter = "\t" if "\t" in first else ","
        tokens = first.split(delimiter)
        names = []
        for token in tokens:
            match = re.fullmatch(r"(.+?)\s*\[([^]]+)\]", token.strip())
            name, unit = (match[1].strip(), match[2]) if match else (token.strip(), None)
            name = column("node" if name.lower() == "node number" else name)
            required_unit = (
                metadata.get("stress_unit")
                if name.startswith("S_")
                else (
                    metadata.get("length_unit")
                    if name in ("x", "y", "z")
                    else (
                        metadata.get("strain_unit")
                        if name.startswith(("E_", "Eel_", "Epl_"))
                        else None
                    )
                )
            )
            if unit and unit != required_unit:
                raise ValueError(
                    f"Workbench column unit {unit!r} conflicts with declared unit {required_unit!r}"
                )
            names.append(name)
        if len(set(names)) != len(names):
            raise ValueError("duplicate Workbench columns")
        if step is None:
            raise ValueError("Workbench result table needs explicit step")
        frame = Frame((label(str(step)), 1), None if time is None else number(str(time)))
        frame.notes.append(f"Workbench table: {path.name}")
        for line in handle:
            if not line.strip():
                continue
            cells = [v.strip() for v in line.strip().split(delimiter)]
            if len(cells) != len(names):
                raise ValueError("ragged Workbench tensor table")
            frame.rows.append(
                {
                    name: label(v) if name == "node_id" else number(v)
                    for name, v in zip(names, cells)
                }
            )
        return _finish(frame, coordinates, material, instance)


def _frames(paths, coordinates, material, instance, metadata, step, time):
    for path in paths:
        with open(path, encoding="utf-8-sig") as handle:
            first = handle.readline()
        if first.lower().startswith("node number"):
            yield _workbench(path, coordinates, material, instance, metadata, step, time)
            continue
        frame, names, reading = None, None, False
        with open(path, encoding="utf-8-sig") as handle:
            for raw in handle:
                line = raw.strip()
                upper = line.upper()
                if "LOCAL COORDINATE" in upper or re.search(r"RSYS\s*[,=]\s*(?!0\b)\S+", upper):
                    raise ValueError(
                        "Ansys results must use global Cartesian RSYS,0, not local coordinates"
                    )
                if any(
                    term in upper
                    for term in (
                        "TOP/BOTTOM",
                        "ELEMENT NODAL",
                        "ELEMENT=",
                        "SECTION POINT",
                        "CYLINDRICAL",
                        "SPHERICAL",
                    )
                ):
                    raise ValueError(
                        "shell/element section-point tables are unsupported; export unique nodal solid tensors"
                    )
                match = re.search(r"LOAD STEP\s*=\s*(\d+)\s+SUBSTEP\s*=\s*(\d+)", upper)
                if match:
                    identity = (int(match[1]), int(match[2]))
                    if frame is not None and frame.identity != identity:
                        yield _finish(frame, coordinates, material, instance)
                        frame = None
                    if frame is None:
                        frame = Frame(identity, None)
                    frame.notes.append(line)
                    reading = False
                    continue
                match = re.match(r"TIME\s*=\s*(\S+)", upper)
                if match and frame is not None:
                    stated_time = number(match[1])
                    if frame.time is not None and frame.time != stated_time:
                        raise ValueError("time conflicts within result set")
                    frame.time = stated_time
                    frame.notes.append(line)
                    continue
                if re.match(r"NODE\s+(?:SX|EPELX|EPPLX|EPTOX)\b", upper):
                    if frame is None:
                        if step is None:
                            raise ValueError("listing lacks load step; provide step")
                        frame = Frame((label(str(step)), 1), time)
                    names = [column(token) for token in line.split()]
                    if len(set(names)) != len(names):
                        raise ValueError("duplicate Ansys columns")
                    reading = True
                    continue
                if reading and line:
                    cells = line.split()
                    if cells[0].isdigit():
                        if len(cells) != len(names):
                            raise ValueError("ragged Ansys tensor table")
                        frame.rows.append(
                            {
                                name: label(v) if name == "node_id" else number(v)
                                for name, v in zip(names, cells)
                            }
                        )
                    elif set(line) <= {"-", " "}:
                        continue
                    elif upper.startswith(
                        ("MINIMUM", "MAXIMUM", "PRINT", "*****", "*** ANSYS", "*** MAPDL")
                    ):
                        reading = False
                    else:
                        raise ValueError("unexpected text in Ansys tensor table")
                elif frame is not None and ("COORDINATE" in upper or "AVERAG" in upper):
                    frame.notes.append(line)
        if frame is not None:
            yield _finish(frame, coordinates, material, instance)


def ansys_listing_to_csv(
    report_paths, out_csv, *, nodes=None, material="", instance="", step=None, time=None, **metadata
):
    """Convert complete global nodal Ansys text tensors with explicit CSV declarations.

    Parameters
    ----------
    report_paths : path or sequence of paths
        PRNSOL listings or complete tab/comma-delimited Workbench tensor tables.
    out_csv : path
        Atomic CSV destination.
    nodes : path, optional
        Global Cartesian NLIST undeformed coordinates.
    material, instance : str, optional
        Labels for this material-specific selection.
    step, time : optional
        Required state identity and optional time when absent from a table.
    **metadata
        CSV header declarations, including units, global frame, averaging, analysis and step role.

    Returns
    -------
    CsvHeader
        The header written to disk.
    """
    paths = paths_of(report_paths)
    from pathlib import Path

    sources = paths + ([Path(nodes)] if nodes is not None else [])
    return write_frames(
        _frames(paths, _nodes(nodes), material, instance, metadata, step, time),
        out_csv,
        solver="Ansys text",
        sources=sources,
        metadata=metadata,
    )
