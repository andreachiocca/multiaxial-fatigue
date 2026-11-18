"""Synthetic FE exports from an analytical thick-walled tube, the oracle of the FE tests.

Elastic, axisymmetric tube under internal pressure ``p_i`` (outer pressure zero), an additional
axial force ``F`` and a torque ``T``. Stress follows the textbook Lame solution plus the axial
membrane stress and Saint-Venant torsion, in ``(r, theta, z)``; strain follows isotropic Hooke's
law (tensor shear). The tensors are rotated to global Cartesian axes with ``z`` along the tube
axis. Units are MPa, mm, N and N mm. The solution excludes end effects, plasticity, residual
stress and geometric change: it is a local gauge-section field, not a rig model.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import numpy as np

from .csvio import CSV_FORMAT
from .csvwriter import FeCsvWriter
from .model import HistoryMeta, LoadCases
from .store import HistoryWriter
from .tensors import CANONICAL, _label_to_ij, to_mm, to_mpa


@dataclass(frozen=True)
class TubeGeometry:
    """Inner and outer radius in mm."""

    a: float
    b: float

    def __post_init__(self):
        if not 0.0 < self.a < self.b:
            raise ValueError("a tube needs 0 < a < b")

    @property
    def area(self) -> float:
        return float(np.pi * (self.b**2 - self.a**2))

    @property
    def polar_moment(self) -> float:
        return float(np.pi * (self.b**4 - self.a**4) / 2.0)


@dataclass(frozen=True)
class TubeMaterial:
    """Young's modulus ``E`` in MPa and Poisson's ratio ``nu``."""

    E: float
    nu: float


@dataclass(frozen=True)
class TubeLoad:
    """One load case: internal pressure ``p_i`` (MPa), additional axial force ``F`` (N), torque ``T`` (N mm)."""

    p_i: float
    F: float
    T: float


def tube_stress_cyl(
    r: np.ndarray, load: TubeLoad, geom: TubeGeometry, *, closed_ends: bool = True
) -> np.ndarray:
    """Stress tensors ``(..., 3, 3)`` in ``(r, theta, z)`` at radii ``r``.

    ``F`` is the force carried by the wall in addition to the pressure thrust, which closed ends
    add as the Lame constant ``A`` to the axial stress.
    """
    r = np.asarray(r, dtype=float)
    a2, b2 = geom.a**2, geom.b**2
    A = load.p_i * a2 / (b2 - a2)
    B = a2 * b2 * load.p_i / (b2 - a2)
    s = np.zeros(r.shape + (3, 3))
    s[..., 0, 0] = A - B / r**2
    s[..., 1, 1] = A + B / r**2
    s[..., 2, 2] = load.F / geom.area + (A if closed_ends else 0.0)
    tau = load.T * r / geom.polar_moment
    s[..., 1, 2] = tau
    s[..., 2, 1] = tau
    return s


def tube_strain_cyl(stress: np.ndarray, material: TubeMaterial) -> np.ndarray:
    """Isotropic Hooke strain (tensor shear) of stress tensors ``(..., 3, 3)``."""
    stress = np.asarray(stress, dtype=float)
    trace = np.trace(stress, axis1=-2, axis2=-1)[..., None, None]
    return ((1.0 + material.nu) * stress - material.nu * trace * np.eye(3)) / material.E


def tube_points(
    geom: TubeGeometry, n_radii: int, n_theta: int, z: float = 0.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Locations on a ring grid: ``(location_id, coordinates, r, theta)``, radius-major.

    IDs start at 1, like solver node labels. A single radius is the mid-wall.
    """
    if n_radii < 1 or n_theta < 1:
        raise ValueError("need at least one radius and one angle")
    radii = (
        np.linspace(geom.a, geom.b, n_radii) if n_radii > 1 else np.array([(geom.a + geom.b) / 2])
    )
    angles = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    r = np.repeat(radii, n_theta)
    theta = np.tile(angles, n_radii)
    xyz = np.stack([r * np.cos(theta), r * np.sin(theta), np.full_like(r, z)], axis=1)
    return np.arange(1, len(r) + 1), xyz, r, theta


def _cylindrical_to_cartesian(T: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Rotate ``(r, theta, z)`` tensors to Cartesian axes by the angle ``theta`` about ``z``."""
    c, s = np.cos(theta), np.sin(theta)
    R = np.zeros(theta.shape + (3, 3))
    R[..., 0, 0], R[..., 0, 1] = c, -s
    R[..., 1, 0], R[..., 1, 1] = s, c
    R[..., 2, 2] = 1.0
    rotated = R @ T @ np.swapaxes(R, -1, -2)
    return 0.5 * (rotated + np.swapaxes(rotated, -1, -2))  # remove roundoff asymmetry


def tube_cartesian(
    r: np.ndarray,
    theta: np.ndarray,
    load: TubeLoad,
    geom: TubeGeometry,
    material: TubeMaterial,
    *,
    closed_ends: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Global Cartesian stress (MPa) and strain at points ``(r, theta)`` for one load case."""
    cyl = tube_stress_cyl(r, load, geom, closed_ends=closed_ends)
    S = _cylindrical_to_cartesian(cyl, np.asarray(theta, dtype=float))
    E = _cylindrical_to_cartesian(tube_strain_cyl(cyl, material), np.asarray(theta, dtype=float))
    return S, E


def tube_meta(**override) -> HistoryMeta:
    """Metadata of a synthetic tube history (stress in MPa, length in mm, Cartesian frame)."""
    kw = dict(
        material_id="synthetic_tube",
        solver="analytical thick tube",
        solver_version="1",
        analysis_step="synthetic",
        adapter="multifat.fe.synthetic",
        adapter_version="1.0",
        source_files=(),
        source_units={"stress": "MPa", "length": "mm"},
        source_frame="cartesian",
        averaging="none",
        location_kind="node",
        notes="Lame and Saint-Venant solution; not a validation of any fatigue method.",
    )
    kw.update(override)
    return HistoryMeta(**kw)


def write_tube_history(
    path,
    geom: TubeGeometry,
    material: TubeMaterial,
    loads: Sequence[TubeLoad],
    *,
    n_radii: int,
    n_theta: int,
    closed_ends: bool = True,
    dtype: str = "float64",
) -> None:
    """Write the universal file for a sequence of load cases, one load case at a time."""
    ids, xyz, r, theta = tube_points(geom, n_radii, n_theta)
    n = len(loads)
    cases = LoadCases(time=np.linspace(0.0, 1.0, n), labels=tuple(f"load{i}" for i in range(n)))
    with HistoryWriter(
        path, tube_meta(), ids, xyz, cases, n, with_strain=True, dtype=dtype
    ) as writer:
        for k, load in enumerate(loads):
            S, E = tube_cartesian(r, theta, load, geom, material, closed_ends=closed_ends)
            writer.write_step(k, S, E)


_VTK_TYPES = {"float64": "Float64", "int64": "Int64", "uint8": "UInt8"}


def _pack(T: np.ndarray, order: Sequence[str], kind: str, shear: str) -> np.ndarray:
    """Tensors ``(P, 3, 3)`` as component columns in ``order`` (engineering shear doubles strain)."""
    cols = []
    for label in order:
        i, j = _label_to_ij(label)
        col = T[:, i, j]
        cols.append(2.0 * col if (kind == "strain" and shear == "engineering" and i != j) else col)
    return np.stack(cols, axis=1)


def _data_array(name: str | None, values: np.ndarray, encoding: str, header_type: str) -> str:
    values = np.ascontiguousarray(values)
    n_comp = values.shape[1] if values.ndim == 2 else 1
    attrs = f'type="{_VTK_TYPES[values.dtype.name]}"'
    if name is not None:
        attrs += f' Name="{name}"'
    attrs += (
        f' NumberOfComponents="{n_comp}" format="{"ascii" if encoding == "ascii" else "binary"}"'
    )
    if encoding == "ascii":
        fmt = "{:.17g}" if values.dtype.kind == "f" else "{:d}"
        body = " ".join(fmt.format(x) for x in values.ravel())
    else:
        raw = values.astype(values.dtype.newbyteorder("<")).tobytes()
        size = len(raw).to_bytes(4 if header_type == "UInt32" else 8, "little")
        if encoding == "base64_separate":  # header as a block of its own
            body = base64.b64encode(size).decode() + base64.b64encode(raw).decode()
        else:  # the layout VTK writes: one stream
            body = base64.b64encode(size + raw).decode()
    return f"<DataArray {attrs}>{body}</DataArray>"


def _vtu_text(points, point_arrays, cell_arrays, encoding, header_type) -> str:
    n = len(points)
    wrap = lambda tag, arrays: (  # noqa: E731
        f"<{tag}>"
        + "".join(_data_array(k, v, encoding, header_type) for k, v in arrays.items())
        + f"</{tag}>"
    )
    cells = (
        "<Cells>"
        + _data_array("connectivity", np.arange(n, dtype=np.int64), encoding, header_type)
        + _data_array("offsets", np.arange(1, n + 1, dtype=np.int64), encoding, header_type)
        + _data_array("types", np.ones(n, dtype=np.uint8), encoding, header_type)
        + "</Cells>"
    )
    return (
        '<?xml version="1.0"?>\n'
        f'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="{header_type}">'
        f'<UnstructuredGrid><Piece NumberOfPoints="{n}" NumberOfCells="{n}">'
        + wrap("PointData", point_arrays)
        + wrap("CellData", cell_arrays)
        + "<Points>" + _data_array(None, points, encoding, header_type) + "</Points>"
        + cells
        + "</Piece></UnstructuredGrid></VTKFile>\n"
    )  # fmt: skip


def write_tube_pvd(
    directory,
    geom: TubeGeometry,
    material: TubeMaterial,
    loads: Sequence[TubeLoad],
    *,
    n_radii: int,
    n_theta: int,
    encoding: Literal["ascii", "base64"] = "ascii",
    header_type: Literal["UInt32", "UInt64"] = "UInt32",
    base64_layout: Literal["joint", "separate"] = "joint",
    frame: Literal["cartesian", "cylindrical"] = "cartesian",
    order: Sequence[str] = CANONICAL,
    shear: Literal["tensor", "engineering"] = "tensor",
    stress_unit: str = "MPa",
    length_unit: str = "mm",
    data_location: Literal["point", "cell"] = "point",
    closed_ends: bool = True,
) -> Path:
    """Write a ``.pvd`` series (one ``.vtu`` point cloud per load case) with chosen conventions.

    ``base64_layout`` chooses between the single base64 stream VTK writes and a header encoded as
    a block of its own. Arrays are ``Stress``, ``Strain`` and ``NodeID``. With ``frame="cylindrical"`` the tensors are
    left in ``(r, theta, z)`` and the importer must be told so. Returns the ``.pvd`` path.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    ids, xyz, r, theta = tube_points(geom, n_radii, n_theta)
    points = xyz / float(to_mm(1.0, length_unit))
    stress_scale = 1.0 / float(to_mpa(1.0, stress_unit))
    kind = (
        "ascii"
        if encoding == "ascii"
        else ("base64" if base64_layout == "joint" else "base64_separate")
    )
    entries = []
    for k, load in enumerate(loads):
        if frame == "cylindrical":
            S = tube_stress_cyl(r, load, geom, closed_ends=closed_ends)
            E = tube_strain_cyl(S, material)
        else:
            S, E = tube_cartesian(r, theta, load, geom, material, closed_ends=closed_ends)
        arrays = {
            "Stress": _pack(S, order, "stress", shear) * stress_scale,
            "Strain": _pack(E, order, "strain", shear),
            "NodeID": ids.astype(np.int64),
        }
        name = f"step_{k:04d}.vtu"
        point_arrays, cell_arrays = (arrays, {}) if data_location == "point" else ({}, arrays)
        text = _vtu_text(points, point_arrays, cell_arrays, kind, header_type)
        (directory / name).write_text(text, encoding="utf-8")
        entries.append(f'<DataSet timestep="{k}" part="0" file="{name}"/>')
    pvd = directory / "tube.pvd"
    pvd.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
        "<Collection>" + "".join(entries) + "</Collection></VTKFile>\n",
        encoding="utf-8",
    )
    return pvd


_CSV_COMPONENTS = ("xx", "yy", "zz", "xy", "yz", "xz")
_ABAQUS_NUMBERS = ("11", "22", "33", "12", "23", "13")
_ANSYS_SUFFIXES = ("X", "Y", "Z", "XY", "YZ", "XZ")


def _csv_tensor_names(flavour: str, group: str) -> list[str]:
    """Column names of a tensor group (``S``, ``Eel``, ``Epl`` or ``E``) in a solver flavour."""
    if flavour == "canonical":
        return [f"{group}_{c}" for c in _CSV_COMPONENTS]
    if flavour == "abaqus":
        prefix = {"S": "S", "Eel": "EE", "Epl": "PE", "E": "E"}[group]
        return [f"{prefix}{n}" for n in _ABAQUS_NUMBERS]
    prefix = {"S": "S", "Eel": "EPEL", "Epl": "EPPL", "E": "EPTO"}[group]
    return [f"{prefix}{x}" for x in _ANSYS_SUFFIXES]


def write_tube_csv(
    path,
    geom: TubeGeometry,
    material: TubeMaterial,
    loads: Sequence[TubeLoad],
    *,
    n_radii: int,
    n_theta: int,
    flavour: Literal["canonical", "abaqus", "ansys"] = "canonical",
    shear: Literal["tensor", "engineering"] = "tensor",
    stress_unit: str = "MPa",
    length_unit: str = "mm",
    order: Literal["step_major", "node_major"] = "step_major",
    with_instance: bool = False,
    with_material: bool = False,
    with_time: bool = False,
    strain: Literal["split", "total", "none"] = "split",
    analysis_type: Literal["linear_elastic", "elastic_plastic"] = "linear_elastic",
    plastic_fraction: float | None = None,
    header: dict | None = None,
    closed_ends: bool = True,
) -> Path:
    """Write the tube as an FE CSV file with the column names of a solver flavour.

    The stress and the Hooke strain are the oracle of the importer tests. ``strain="split"``
    writes the elastic strain and a plastic strain ``plastic_fraction`` times as large (default
    0.25 for ``elastic_plastic``, else 0), ``"total"`` the sum of the two, ``"none"`` no strain
    columns. Engineering shear doubles the off-diagonal strain columns only. Returns the path.
    """
    if flavour not in ("canonical", "abaqus", "ansys"):
        raise ValueError(f"flavour must be canonical, abaqus or ansys, got {flavour!r}")
    if order not in ("step_major", "node_major"):
        raise ValueError(f"order must be step_major or node_major, got {order!r}")
    if plastic_fraction is None:
        plastic_fraction = 0.25 if analysis_type == "elastic_plastic" else 0.0
    ids, xyz, r, theta = tube_points(geom, n_radii, n_theta)
    n_nodes, n_steps = len(ids), len(loads)
    coords = xyz / float(to_mm(1.0, length_unit))
    stress_scale = 1.0 / float(to_mpa(1.0, stress_unit))
    head = ["node_id", "instance", "material", "step", "time", "x", "y", "z"]
    if flavour != "canonical":
        head = [h.upper() for h in head]
    node_c, inst_c, mat_c, step_c, time_c, x_c, y_c, z_c = head
    columns = [node_c]
    columns += [inst_c] if with_instance else []
    columns += [mat_c] if with_material else []
    columns += [step_c] + ([time_c] if with_time else []) + [x_c, y_c, z_c]
    stress_names = _csv_tensor_names(flavour, "S")
    columns += stress_names
    groups: list[tuple[str, list[str]]] = []
    if strain == "split":
        groups = [
            ("Eel", _csv_tensor_names(flavour, "Eel")),
            ("Epl", _csv_tensor_names(flavour, "Epl")),
        ]
    elif strain == "total":
        groups = [("E", _csv_tensor_names(flavour, "E"))]
    elif strain != "none":
        raise ValueError(f"strain must be split, total or none, got {strain!r}")
    for _, names in groups:
        columns += names
    fields = {
        "format": CSV_FORMAT,
        "solver": f"analytical thick tube ({flavour} column names)",
        "solver_version": "1",
        "stress_unit": stress_unit,
        "length_unit": length_unit,
        "coordinate_system": "global_cartesian",
        "averaging": "none",
        "analysis_type": analysis_type,
        "step_role": "time_sample",
    }
    if analysis_type == "elastic_plastic":
        fields["cycle"] = "stabilised"
    if groups:
        fields["strain_convention"] = shear
        fields["strain_unit"] = "1"
    if with_time:
        fields["time_unit"] = "s"
    fields["n_nodes"] = n_nodes
    fields["n_steps"] = n_steps
    fields["notes"] = "Lame and Saint-Venant solution; not a validation of any fatigue method."
    fields.update(header or {})
    times = np.linspace(0.0, 1.0, n_steps) if n_steps > 1 else np.zeros(1)
    pairs = [(i, j) for i, j in ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))]

    def tensors_of(load):
        S, E = tube_cartesian(r, theta, load, geom, material, closed_ends=closed_ends)
        out = {"S": S * stress_scale}
        if strain == "split":
            out["Eel"], out["Epl"] = E, plastic_fraction * E
        elif strain == "total":
            out["E"] = (1.0 + plastic_fraction) * E
        return out

    def row_of(node: int, step: int, data: dict) -> dict:
        row: dict = {node_c: int(ids[node])}
        if with_instance:
            row[inst_c] = "TUBE-1"
        if with_material:
            row[mat_c] = "STEEL"
        row[step_c] = step + 1
        if with_time:
            row[time_c] = float(times[step])
        row[x_c], row[y_c], row[z_c] = (float(v) for v in coords[node])
        for name, (i, j) in zip(stress_names, pairs):
            row[name] = float(data["S"][node, i, j])
        for group, names in groups:
            for name, (i, j) in zip(names, pairs):
                value = float(data[group][node, i, j])
                row[name] = value * (2.0 if shear == "engineering" and i != j else 1.0)
        return row

    per_step = [tensors_of(load) for load in loads]
    with FeCsvWriter(path, fields, columns=columns) as writer:
        if order == "step_major":
            for step, data in enumerate(per_step):
                writer.write_rows(row_of(node, step, data) for node in range(n_nodes))
        else:
            for node in range(n_nodes):
                writer.write_rows(row_of(node, step, data) for step, data in enumerate(per_step))
    return Path(path)


def main(argv: Sequence[str] | None = None) -> int:
    """Write a synthetic tube history as a universal file, a ``.pvd`` series and a CSV example.

    ``python -m multifat.fe.synthetic --out DIR`` writes ``DIR/tube_history.h5`` and
    ``DIR/vtk/tube.pvd``; the two describe the same three load cases and import back to the same
    tensors. ``--csv DIR`` writes ``DIR/tube_two_steps.csv``, the two-step example of the CSV
    exchange format (``docs/fe-csv-format.md``).
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m multifat.fe.synthetic", description=main.__doc__
    )
    parser.add_argument("--out", help="directory for the HDF5 and VTK files (created if missing)")
    parser.add_argument("--csv", help="directory for the two-step CSV example (created if missing)")
    parser.add_argument("--n-radii", type=int, default=3)
    parser.add_argument("--n-theta", type=int, default=12)
    args = parser.parse_args(argv)
    if args.out is None and args.csv is None:
        parser.error("give --out, --csv or both")
    geom, material = TubeGeometry(a=6.0, b=7.5), TubeMaterial(E=210000.0, nu=0.3)
    if args.out is not None:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        loads = [
            TubeLoad(0.0, 0.0, 0.0),
            TubeLoad(50.0, 12000.0, 0.0),
            TubeLoad(50.0, 12000.0, 8e4),
        ]
        kw = dict(n_radii=args.n_radii, n_theta=args.n_theta)
        write_tube_history(out / "tube_history.h5", geom, material, loads, **kw)
        pvd = write_tube_pvd(out / "vtk", geom, material, loads, encoding="base64", **kw)
        print(f"wrote {out / 'tube_history.h5'} and {pvd}")
    if args.csv is not None:
        directory = Path(args.csv)
        directory.mkdir(parents=True, exist_ok=True)
        two_steps = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 12000.0, 8e4)]
        path = write_tube_csv(
            directory / "tube_two_steps.csv", geom, material, two_steps,
            n_radii=2, n_theta=2, strain="total",
        )  # fmt: skip
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
