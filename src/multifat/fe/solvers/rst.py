"""Optional native Ansys solid stress conversion using ansys-mapdl-reader 0.56.

Repository choices:
Only real static/transient, small-displacement, single-material solid stress results with default
element options are accepted. The reader rotates element stresses globally before nodal averaging;
its strain rotation is not equivalent, so this adapter exports stress only. Native ODB needs Abaqus
and is outside this adapter's scope.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .common import Frame, column, label, paths_of, write_frames


def _read_binary(path):
    try:
        from ansys.mapdl.reader import read_binary
    except ImportError as exc:
        raise ImportError("native RST conversion needs the optional multifat[ansys] extra") from exc
    return read_binary(path)


def ansys_rst_to_csv(result_path, out_csv, *, nodes=None, material="", instance="", **metadata):
    """Write global averaged solid stresses from an RST to exchange CSV.

    Parameters
    ----------
    result_path, out_csv : path
        Source native result and atomic CSV destination.
    nodes : sequence of int, optional
        Explicit node selection. Missing or non-finite stress is refused, including midside nodes
        without stored results. No node is silently discarded.
    material, instance : str, optional
        Labels for the single-material model.
    **metadata
        CSV units, analysis_type, step_role, coordinate_system and time_unit declarations.
        Only analysis_type='linear_elastic' is supported; no native strains are exported.

    Returns
    -------
    CsvHeader
        The validated header written to disk.
    """
    path = paths_of(result_path)[0]
    result = _read_binary(path)
    info = result._resultheader
    if result._is_cyclic:
        raise ValueError("cyclic symmetry RST is unsupported")
    if info["kan"] not in (0, 4):
        raise ValueError("only real static or transient RST states are supported")
    if info["cpxrst"]:
        raise ValueError("complex RST data is unsupported")
    if info["nlgeom"]:
        raise ValueError("geometrically nonlinear RST frames are unsupported")
    if metadata.get("analysis_type") != "linear_elastic":
        raise ValueError("native stress-only conversion requires analysis_type='linear_elastic'")
    if len(set(result.mesh.material_type)) != 1:
        raise ValueError("native nodal averaging across multiple materials is refused")
    if (
        not set(result.mesh.ekey[:, 1]).issubset({182, 183, 185, 186, 187})
        or result.mesh.key_option
    ):
        raise ValueError(
            "only solid elements SOLID185/186/187 and plane elements PLANE182/183 are supported"
        )
    units = info["units"]
    if units != "User Defined":
        if units != "SI" or (metadata.get("stress_unit"), metadata.get("length_unit")) != (
            "Pa",
            "m",
        ):
            raise ValueError(
                f"RST unit system {units!r} conflicts with declarations or is unsupported"
            )
    if metadata.get("averaging", "nodal_average_global") != "nodal_average_global":
        raise ValueError("native reader uses nodal_average_global")
    coordinates = {
        int(n): tuple(float(v) for v in x[:3]) for n, x in zip(result.mesh.nnum, result.mesh.nodes)
    }
    selected = None if nodes is None else [label(str(n)) for n in nodes]
    if selected is not None and (not selected or len(set(selected)) != len(selected)):
        raise ValueError("node selection must be nonempty and unique")

    def frames():
        for i in range(result.nsets):
            ids, stress = result.nodal_stress(i)
            if stress.shape != (len(ids), 6) or len(set(ids)) != len(ids):
                raise ValueError("RST nodal stress needs unique labels and six components")
            lookup = {int(n): j for j, n in enumerate(ids)}
            requested = list(lookup) if selected is None else selected
            if any(n not in lookup or n not in coordinates for n in requested):
                raise ValueError("selected node lacks coordinates or stress")
            # Preserve reader order, even for an unordered caller selection.
            keep = set(requested)
            frame = Frame(
                tuple(int(v) for v in info["ls_table"][i, :2]), float(result.time_values[i])
            )
            frame.notes.append(
                f"reader=ansys-mapdl-reader 0.56; global Cartesian stress; selected {len(keep)} of {len(ids)} nodes; material number {next(iter(set(result.mesh.material_type)))}"
            )
            for n in ids:
                node = int(n)
                if node not in keep:
                    continue
                values = stress[lookup[node]]
                if not np.isfinite(values).all():
                    raise ValueError(
                        "non-finite native stress; explicitly select nodes with stored corner stresses"
                    )
                row = dict(node_id=node, instance=instance, material=material)
                row.update(zip(("x", "y", "z"), coordinates[node]))
                row.update(
                    (column(k), float(v))
                    for k, v in zip(("SX", "SY", "SZ", "SXY", "SYZ", "SXZ"), values)
                )
                frame.rows.append(row)
            yield frame

    return write_frames(
        frames(),
        out_csv,
        solver="Ansys RST",
        sources=[Path(path)],
        metadata={
            "solver_version": info["verstring"],
            "averaging": "nodal_average_global",
            **metadata,
        },
    )
