"""VTK XML importer (``.vtu`` and ``.pvd`` time series) for local tensor histories.

Supported: unstructured grids with ASCII or inline base64 ``binary`` data arrays (uncompressed,
``UInt32`` or ``UInt64`` header, little endian), point data or cell data (cell values sit at the
cell centroid). Appended data and compressed data are rejected with a diagnostic, never guessed.
Component order, shear convention, stress unit, length unit and the source frame are arguments:
nothing is inferred. Only the standard library and numpy are needed; writing the universal file
needs ``h5py``.
"""

from __future__ import annotations

import base64
import binascii
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Sequence

import numpy as np

from ._atomic import atomic_destination
from .model import HistoryMeta, LoadCases, SourceFile, sha256_of
from .tensors import (
    components_to_sym3,
    cylindrical_rotations,
    rotate_to_global,
    to_mm,
    to_mpa,
)

_DTYPES = {
    "Float32": "<f4", "Float64": "<f8", "Int8": "i1", "UInt8": "u1", "Int16": "<i2",
    "UInt16": "<u2", "Int32": "<i4", "UInt32": "<u4", "Int64": "<i8", "UInt64": "<u8",
}  # fmt: skip
_HEADER_BYTES = {"UInt32": 4, "UInt64": 8}


@dataclass(frozen=True)
class CylindricalFrame:
    """Source tensors are in ``(r, theta, z)`` about ``axis`` through ``origin`` (source length unit)."""

    axis: tuple[float, float, float] = (0.0, 0.0, 1.0)
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True)
class RotationArray:
    """Source tensors are in per-point frames given by a 9-component data array.

    The array holds the row-major rotation matrix whose columns are the local axes in global
    coordinates.
    """

    name: str


@dataclass(frozen=True)
class TensorSpec:
    """Where a tensor field lives in the file and how to read it.

    ``unit`` is the stress unit and is required for stress (nothing is assumed); it is ignored for
    strain, which is dimensionless. ``shear`` only matters for strain.
    """

    array: str
    order: Sequence[str]
    unit: str | None = None
    shear: Literal["tensor", "engineering"] = "tensor"
    location: Literal["point", "cell"] = "point"


@dataclass
class VtuData:
    """Contents of one ``.vtu`` file."""

    points: np.ndarray
    point_data: dict[str, np.ndarray]
    cell_data: dict[str, np.ndarray]
    cell_centroids: np.ndarray | None


@dataclass
class StepData:
    """One load step converted to the repository conventions (MPa, mm, global Cartesian)."""

    location_id: np.ndarray
    coordinates: np.ndarray
    S: np.ndarray
    E: np.ndarray | None


def _binary_payload(text: str, header_bytes: int) -> bytes | None:
    """Data bytes of an uncompressed inline base64 array, or ``None`` if the layout is unknown.

    VTK writes the byte-count header and the data as one base64 stream; some writers encode the
    header as a block of its own. Both are accepted, and each is confirmed by the byte count.
    """
    try:
        whole = base64.b64decode(text)
    except binascii.Error:
        whole = b""
    if (
        len(whole) >= header_bytes
        and int.from_bytes(whole[:header_bytes], "little") == len(whole) - header_bytes
    ):
        return whole[header_bytes:]
    n_chars = -(-header_bytes // 3) * 4
    try:
        size = int.from_bytes(base64.b64decode(text[:n_chars]), "little")
        raw = base64.b64decode(text[n_chars:])
    except binascii.Error:
        return None
    return raw if len(raw) == size else None


def _decode(element: ET.Element, header_type: str) -> np.ndarray:
    fmt = element.get("format")
    kind = element.get("type")
    if kind not in _DTYPES:
        raise ValueError(f"unsupported DataArray type {kind!r}")
    dtype = np.dtype(_DTYPES[kind])
    text = (element.text or "").strip()
    if fmt == "ascii":
        values = np.array(text.split(), dtype=dtype) if text else np.empty(0, dtype=dtype)
    elif fmt == "binary":
        raw = _binary_payload(text, _HEADER_BYTES[header_type])
        if raw is None:
            raise ValueError(
                f"binary DataArray {element.get('Name')!r}: the byte count in the header does not "
                "match the data (compressed or malformed data are not supported)"
            )
        values = np.frombuffer(raw, dtype=dtype)
    else:
        raise ValueError(
            f"unsupported DataArray format {fmt!r} (appended and compressed data are not supported)"
        )
    n_comp = int(element.get("NumberOfComponents", "1"))
    return values.reshape(-1, n_comp) if n_comp > 1 else values


def read_vtu(path: str | Path) -> VtuData:
    """Parse one ``.vtu`` file; unsupported encodings raise ``ValueError`` naming the feature."""
    root = ET.parse(path).getroot()
    if root.tag != "VTKFile" or root.get("type") != "UnstructuredGrid":
        raise ValueError(f"{path}: not a VTK XML UnstructuredGrid file")
    if root.get("compressor"):
        raise ValueError(f"{path}: compressed VTK data ({root.get('compressor')}) is not supported")
    if root.get("byte_order", "LittleEndian") != "LittleEndian":
        raise ValueError(f"{path}: only little-endian files are supported")
    if root.find("AppendedData") is not None:
        raise ValueError(f"{path}: appended data is not supported")
    header_type = root.get("header_type", "UInt32")
    if header_type not in _HEADER_BYTES:
        raise ValueError(f"{path}: unsupported header_type {header_type!r}")
    pieces = root.findall("UnstructuredGrid/Piece")
    if len(pieces) != 1:
        raise ValueError(f"{path}: expected one Piece, found {len(pieces)}")
    piece = pieces[0]
    n_points = int(piece.get("NumberOfPoints"))
    n_cells = int(piece.get("NumberOfCells", "0"))

    def arrays(tag: str) -> dict[str, np.ndarray]:
        node = piece.find(tag)
        if node is None:
            return {}
        return {a.get("Name"): _decode(a, header_type) for a in node.findall("DataArray")}

    points = _decode(piece.find("Points/DataArray"), header_type).astype(float)
    if points.shape != (n_points, 3):
        raise ValueError(f"{path}: Points has shape {points.shape}, expected ({n_points}, 3)")
    centroids = None
    cells = piece.find("Cells")
    if cells is not None and n_cells:
        named = {a.get("Name"): _decode(a, header_type) for a in cells.findall("DataArray")}
        conn, offsets = named["connectivity"].astype(int), named["offsets"].astype(int)
        starts = np.concatenate([[0], offsets[:-1]])
        centroids = (
            np.add.reduceat(points[conn], starts, axis=0)
            / np.diff(np.concatenate([starts, [offsets[-1]]]))[:, None]
        )
    return VtuData(points, arrays("PointData"), arrays("CellData"), centroids)


def read_pvd(path: str | Path) -> list[tuple[float, Path]]:
    """Time series of a ``.pvd`` collection as ``(timestep, file)``, in timestep order."""
    path = Path(path)
    root = ET.parse(path).getroot()
    entries = list(root.iter("DataSet"))
    if not entries:
        raise ValueError(f"{path}: no DataSet entries")
    if any(d.get("part", "0") != "0" for d in entries):
        raise ValueError(f"{path}: multi-part collections (part != 0) are not supported")
    steps = [(float(d.get("timestep", "0")), path.parent / d.get("file")) for d in entries]
    if len({t for t, _ in steps}) != len(steps):
        raise ValueError(f"{path}: repeated timestep: each load case needs its own timestep")
    return sorted(steps, key=lambda s: s[0])


def _field(vtu: VtuData, tensor_spec: TensorSpec) -> np.ndarray:
    data = vtu.point_data if tensor_spec.location == "point" else vtu.cell_data
    if tensor_spec.array not in data:
        raise KeyError(
            f"array {tensor_spec.array!r} not found in {tensor_spec.location} data; available: {sorted(data)}"
        )
    a = data[tensor_spec.array]
    return a.reshape(len(a), -1).astype(float)


def read_step(
    vtu: VtuData,
    stress: TensorSpec,
    strain: TensorSpec | None,
    *,
    frame: Literal["cartesian"] | CylindricalFrame | RotationArray = "cartesian",
    length_unit: str,
    id_array: str | None = None,
) -> StepData:
    """Convert one parsed file to MPa, mm, tensor strain and global Cartesian axes."""
    if strain is not None and strain.location != stress.location:
        raise ValueError("stress and strain must both be point data or both cell data")
    S = components_to_sym3(_field(vtu, stress), stress.order, kind="stress")
    if stress.unit is None:
        raise ValueError("the stress unit must be given (TensorSpec.unit): none is assumed")
    S = to_mpa(S, stress.unit)
    E = None
    if strain is not None:
        E = components_to_sym3(_field(vtu, strain), strain.order, kind="strain", shear=strain.shear)
    m = S.shape[0]
    if stress.location == "point":
        xyz = to_mm(vtu.points, length_unit)
    elif vtu.cell_centroids is None:
        raise ValueError("cell data needs a Cells section to place the cell centroids")
    else:
        xyz = to_mm(vtu.cell_centroids, length_unit)
    if xyz.shape[0] != m:
        raise ValueError(f"{stress.location} count {m} does not match {xyz.shape[0]} locations")
    if id_array is None:
        ids = np.arange(m)
    else:
        ids = _field(vtu, TensorSpec(id_array, (), location=stress.location)).ravel()
        if np.any(ids != np.round(ids)):
            raise ValueError(f"ID array {id_array!r} is not integer-valued")
        ids = ids.astype(np.int64)
    if len(np.unique(ids)) != m:
        raise ValueError("duplicate location IDs")
    if frame != "cartesian":
        if isinstance(frame, CylindricalFrame):
            R = cylindrical_rotations(xyz, frame.axis, to_mm(np.asarray(frame.origin), length_unit))
        elif isinstance(frame, RotationArray):
            R = _field(vtu, TensorSpec(frame.name, (), location=stress.location)).reshape(m, 3, 3)
            if not np.allclose(R @ np.swapaxes(R, -1, -2), np.eye(3), atol=1e-6):
                raise ValueError(f"rotation array {frame.name!r} is not orthonormal")
            if np.any(np.linalg.det(R) < 0):
                raise ValueError(
                    f"rotation array {frame.name!r} has determinant -1 (a reflection, not a rotation)"
                )
        else:
            raise ValueError(f"unknown frame {frame!r}")
        S = rotate_to_global(S, R)
        E = None if E is None else rotate_to_global(E, R)
        S = 0.5 * (S + np.swapaxes(S, -1, -2))
        E = None if E is None else 0.5 * (E + np.swapaxes(E, -1, -2))
    return StepData(ids, xyz, S, E)


def _frame_label(frame) -> str:
    if frame == "cartesian":
        return "cartesian"
    if isinstance(frame, CylindricalFrame):
        return f"cylindrical(axis={tuple(frame.axis)}, origin={tuple(frame.origin)})"
    return f"rotation_array:{frame.name}"


def import_pvd(
    pvd_path: str | Path,
    out_path: str | Path,
    *,
    stress: TensorSpec,
    strain: TensorSpec | None,
    id_array: str | None,
    material_id: str,
    solver: str,
    solver_version: str,
    analysis_step: str,
    length_unit: str,
    frame: Literal["cartesian"] | CylindricalFrame | RotationArray = "cartesian",
    averaging: str = "none",
    location_kind: str | None = None,
    periodic: bool = False,
) -> HistoryMeta:
    """Convert a ``.pvd`` series to the universal file, one load step at a time.

    Peak memory is one load step. Location IDs and coordinates must agree across steps; the
    partial output file is removed if anything fails. Returns the metadata written.
    """
    from .store import HistoryWriter

    pvd_path, out_path = Path(pvd_path), Path(out_path)
    steps = read_pvd(pvd_path)
    kind = location_kind or ("node" if stress.location == "point" else "element_centroid")
    conventions = {"stress_order": ",".join(stress.order), "location": stress.location}
    if strain is not None:
        conventions.update(strain_order=",".join(strain.order), strain_shear=strain.shear)
    meta = HistoryMeta(
        material_id=material_id,
        solver=solver,
        solver_version=solver_version,
        analysis_step=analysis_step,
        adapter="multifat.fe.vtk",
        adapter_version="1.0",
        source_files=tuple(
            SourceFile(p.name, sha256_of(p)) for p in [pvd_path] + [f for _, f in steps]
        ),
        source_units={"stress": stress.unit, "length": length_unit},
        source_frame=_frame_label(frame),
        averaging=averaging,
        location_kind=kind,
        periodic=periodic,
        created=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        notes="" if id_array else "location IDs are VTK point or cell indices",
        strain_kind="total_mechanical" if strain is not None else None,
        source_conventions=conventions,
    )
    first = read_step(
        read_vtu(steps[0][1]),
        stress,
        strain,
        frame=frame,
        length_unit=length_unit,
        id_array=id_array,
    )
    cases = LoadCases(time=np.array([t for t, _ in steps]), labels=tuple(f.name for _, f in steps))
    with atomic_destination(out_path, [pvd_path] + [p for _, p in steps]) as temporary:
        with HistoryWriter(
            temporary, meta, first.location_id, first.coordinates, cases, len(steps),
            with_strain=strain is not None,
        ) as writer:  # fmt: skip
            for k, (_, path) in enumerate(steps):
                got = first if k == 0 else read_step(
                    read_vtu(path), stress, strain, frame=frame, length_unit=length_unit,
                    id_array=id_array,
                )  # fmt: skip
                if got.location_id.shape != first.location_id.shape:
                    raise ValueError(
                        f"{stress.location} count differs between load steps "
                        f"({got.location_id.shape[0]} in {path.name}, {first.location_id.shape[0]} first)"
                    )
                if not np.array_equal(got.location_id, first.location_id):
                    raise ValueError(f"location IDs differ between load steps ({path.name})")
                scale = max(float(np.max(np.abs(first.coordinates))), 1.0)
                if not np.allclose(got.coordinates, first.coordinates, rtol=0, atol=1e-9 * scale):
                    raise ValueError(f"coordinates differ between load steps ({path.name})")
                writer.write_step(k, got.S, got.E)
    return meta
