"""The FE CSV exchange format: header, named columns with solver aliases and the streaming importer.

A file starts with ``# key: value`` header lines, then one comma-separated column line, then one
row per location and load step. Nothing is inferred: units, frame, averaging and strain
conventions are header keys, columns are matched by name (never by position) and a canonical column
together with its alias is refused. The format is specified in ``docs/fe-csv-format.md``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from ._atomic import atomic_destination
from .model import HistoryMeta, LoadCases, SourceFile, sha256_of
from .store import HistoryWriter
from .tensors import (
    CANONICAL,
    LENGTH_UNITS,
    STRAIN_UNITS,
    STRESS_UNITS,
    components_to_sym3,
    elastic_strain,
    to_mm,
    to_mpa,
    to_strain,
)

CSV_FORMAT = "multifat-fe-csv 1.0"
_FORMAT_NAME, _FORMAT_MAJOR = CSV_FORMAT.split()[0], CSV_FORMAT.split()[1].split(".")[0]

REQUIRED_HEADER_KEYS = (
    "format",
    "solver",
    "solver_version",
    "stress_unit",
    "length_unit",
    "coordinate_system",
    "averaging",
    "analysis_type",
    "step_role",
)
CSV_AVERAGING = ("none", "nodal_average_within_material", "nodal_average_global", "extrapolated")
CSV_ANALYSIS_TYPES = ("linear_elastic", "elastic_plastic")
CSV_STEP_ROLES = ("time_sample", "unit_case")
STRAIN_CONVENTIONS = ("tensor", "engineering")
_OPTIONAL_TEXT = ("load_description", "model", "source", "notes")
_OPTIONAL_KEYS = (
    "strain_convention",
    "strain_unit",
    "cycle",
    "time_unit",
    "temperature_unit",
    "n_nodes",
    "n_steps",
) + _OPTIONAL_TEXT


def _choice(key: str, value: str, allowed: Iterable[str]) -> str:
    allowed = tuple(allowed)
    if value not in allowed:
        raise ValueError(f"header key {key!r}: {value!r} is not one of {allowed}")
    return value


def _count(key: str, value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise ValueError(f"header key {key!r}: {value!r} is not an integer") from None
    if number < 1:
        raise ValueError(f"header key {key!r}: must be a positive integer, got {number}")
    return number


@dataclass(frozen=True)
class CsvHeader:
    """Validated header of an FE CSV file (see ``docs/fe-csv-format.md``)."""

    format: str
    solver: str
    solver_version: str
    stress_unit: str
    length_unit: str
    coordinate_system: str
    averaging: str
    analysis_type: str
    step_role: str
    strain_convention: str | None = None
    strain_unit: str | None = None
    cycle: str | None = None
    time_unit: str | None = None
    temperature_unit: str | None = None
    n_nodes: int | None = None
    n_steps: int | None = None
    load_description: str = ""
    model: str = ""
    source: str = ""
    notes: str = ""
    extra: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: Mapping[str, object]) -> "CsvHeader":
        """Validate raw ``key -> text`` pairs and build the header."""
        raw = {str(k): str(v).strip() for k, v in data.items()}
        missing = [k for k in REQUIRED_HEADER_KEYS if not raw.get(k)]
        if missing:
            raise ValueError(f"missing required header key(s): {', '.join(missing)}")
        name, _, version = raw["format"].partition(" ")
        if name != _FORMAT_NAME or version.split(".")[0] != _FORMAT_MAJOR:
            raise ValueError(
                f"header key 'format': {raw['format']!r} is not supported; expected {CSV_FORMAT!r}"
                " (or a later 1.x minor version)"
            )
        _choice("stress_unit", raw["stress_unit"], STRESS_UNITS)
        _choice("length_unit", raw["length_unit"], LENGTH_UNITS)
        _choice("coordinate_system", raw["coordinate_system"], ("global_cartesian",))
        _choice("averaging", raw["averaging"], CSV_AVERAGING)
        _choice("analysis_type", raw["analysis_type"], CSV_ANALYSIS_TYPES)
        _choice("step_role", raw["step_role"], CSV_STEP_ROLES)
        kw: dict[str, object] = {k: raw[k] for k in REQUIRED_HEADER_KEYS}
        if raw.get("strain_convention"):
            kw["strain_convention"] = _choice(
                "strain_convention", raw["strain_convention"], STRAIN_CONVENTIONS
            )
        if raw.get("strain_unit"):
            kw["strain_unit"] = _choice("strain_unit", raw["strain_unit"], STRAIN_UNITS)
        cycle = raw.get("cycle")
        if cycle:
            if cycle not in ("stabilised", "not_applicable") and not (
                cycle.isdigit() and int(cycle) >= 1
            ):
                raise ValueError(
                    f"header key 'cycle': {cycle!r} is not 'stabilised', 'not_applicable' or a "
                    "positive integer"
                )
            kw["cycle"] = cycle
        if raw["analysis_type"] == "elastic_plastic" and cycle != "stabilised":
            raise ValueError(
                "header key 'cycle': an elastic_plastic file must declare cycle: stabilised "
                f"(got {cycle!r}); the shakedown or stabilised cycle is what a fatigue assessment "
                "needs"
            )
        if raw["step_role"] == "unit_case" and raw["analysis_type"] != "linear_elastic":
            raise ValueError(
                "header key 'step_role': unit_case superposition needs analysis_type: linear_elastic"
            )
        for key in ("time_unit", "temperature_unit"):
            if raw.get(key):
                kw[key] = raw[key]
        for key in ("n_nodes", "n_steps"):
            if raw.get(key):
                kw[key] = _count(key, raw[key])
        for key in _OPTIONAL_TEXT:
            if raw.get(key):
                kw[key] = raw[key]
        known = set(REQUIRED_HEADER_KEYS) | set(_OPTIONAL_KEYS)
        kw["extra"] = {k: v for k, v in raw.items() if k not in known}
        return cls(**kw)


def read_header(path: str | Path) -> tuple[CsvHeader, int]:
    """Parse the ``# key: value`` lines at the top of a file.

    Returns
    -------
    tuple
        The validated header and the number of header lines (the column line follows).

    Notes
    -----
    Only the first colon of a line separates key and value, so a free-text value may contain
    colons. An empty ``#`` line is ignored. Anything else that starts with ``#`` and has no colon,
    a repeated key, or a file without a column line is refused.
    """
    raw: dict[str, str] = {}
    n_lines = 0
    column_line = None
    with open(path, encoding="utf-8-sig", newline="") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.startswith("#"):
                column_line = line
                break
            n_lines = number
            text = line[1:].strip()
            if not text:
                continue
            key, colon, value = text.partition(":")
            key = key.strip()
            if not colon or not key:
                raise ValueError(f"{path}: line {number}: header lines are '# key: value'")
            if key in raw:
                raise ValueError(f"{path}: line {number}: header key {key!r} given twice")
            raw[key] = value.strip()
    if not column_line or not column_line.strip():
        raise ValueError(f"{path}: no column line after the {n_lines} header lines")
    return CsvHeader.from_mapping(raw), n_lines


# --- columns -----------------------------------------------------------------------------------

_COMPONENTS = ("xx", "yy", "zz", "xy", "yz", "xz")
_ABAQUS_NUMBERS = dict(zip(_COMPONENTS, ("11", "22", "33", "12", "23", "13")))
_ANSYS_SUFFIX = dict(zip(_COMPONENTS, ("X", "Y", "Z", "XY", "YZ", "XZ")))
_GROUPS = {"S": "stress", "Eel": "elastic strain", "Epl": "plastic strain", "E": "total strain"}

CANONICAL_COLUMNS = (
    ("node_id", "instance", "material", "step", "time", "x", "y", "z")
    + tuple(f"S_{c}" for c in _COMPONENTS)
    + tuple(f"Eel_{c}" for c in _COMPONENTS)
    + tuple(f"Epl_{c}" for c in _COMPONENTS)
    + tuple(f"E_{c}" for c in _COMPONENTS)
    + ("T",)
)
REQUIRED_COLUMNS = ("node_id", "step", "x", "y", "z") + tuple(f"S_{c}" for c in _COMPONENTS)


def _alias_table() -> dict[str, str]:
    table = {name.lower(): name for name in CANONICAL_COLUMNS}
    table["node"] = "node_id"
    table["temp"] = "T"
    for comp in _COMPONENTS:
        n, a = _ABAQUS_NUMBERS[comp], _ANSYS_SUFFIX[comp].lower()
        table[f"s{n}"] = f"S_{comp}"  # Abaqus S11 .. S13: S12 -> xy, S23 -> yz, S13 -> xz
        table[f"ee{n}"] = f"Eel_{comp}"
        table[f"pe{n}"] = f"Epl_{comp}"
        table[f"e{n}"] = f"E_{comp}"
        table[f"s{a}"] = f"S_{comp}"  # Ansys SX .. SXZ (short Sxx .. Sxz coincide for shear)
        table[f"epel{a}"] = f"Eel_{comp}"
        table[f"eppl{a}"] = f"Epl_{comp}"
        table[f"epto{a}"] = f"E_{comp}"
        table[f"s{comp}"] = f"S_{comp}"
    return table


COLUMN_ALIASES: dict[str, str] = _alias_table()


@dataclass(frozen=True)
class ColumnMap:
    """Canonical column name -> the name used in the file, plus the columns left out."""

    names: dict[str, str]
    ignored: tuple[str, ...] = ()

    def _group(self, prefix: str) -> bool:
        return f"{prefix}_xx" in self.names

    @property
    def has_total(self) -> bool:
        return self._group("E")

    @property
    def has_elastic(self) -> bool:
        return self._group("Eel")

    @property
    def has_plastic(self) -> bool:
        return self._group("Epl")

    @property
    def has_split(self) -> bool:
        return self.has_elastic and self.has_plastic

    @property
    def has_elastic_only(self) -> bool:
        return self.has_elastic and not self.has_plastic and not self.has_total

    @property
    def has_plastic_only(self) -> bool:
        return self.has_plastic and not self.has_elastic and not self.has_total

    @property
    def has_strain(self) -> bool:
        return self.has_total or self.has_elastic or self.has_plastic

    @property
    def has_temperature(self) -> bool:
        return "T" in self.names

    @property
    def has_time(self) -> bool:
        return "time" in self.names

    @property
    def has_instance(self) -> bool:
        return "instance" in self.names

    @property
    def has_material(self) -> bool:
        return "material" in self.names

    def check_header(self, header: CsvHeader) -> None:
        """Refuse a header that lacks what the columns need (units, convention, strain)."""
        if self.has_strain:
            if header.strain_convention is None:
                raise ValueError(
                    "header key 'strain_convention' (tensor or engineering) is required when the "
                    "file has strain columns"
                )
            if header.strain_unit is None:
                raise ValueError(
                    "header key 'strain_unit' (1, mm/mm, m/m, percent, microstrain) is required "
                    "when the file has strain columns"
                )
        if self.has_time and header.time_unit is None:
            raise ValueError("header key 'time_unit' is required when the file has a time column")
        if header.analysis_type == "elastic_plastic" and not self.has_strain:
            raise ValueError(
                "an elastic_plastic file needs strain columns (Eel and Epl, or a total E): "
                "elastic strain cannot be derived from stress under plasticity"
            )


def resolve_columns(columns: Iterable[str], *, ignore: Iterable[str] = ()) -> ColumnMap:
    """Match the file's column names to the canonical ones, by name and case-insensitively.

    Parameters
    ----------
    columns : iterable of str
        The column line of the file.
    ignore : iterable of str
        Column names to leave out (case-insensitive), for exporter extras the format does not use.

    Returns
    -------
    ColumnMap
        Canonical-to-file names; unknown columns, a canonical column together with its alias,
        missing required columns, an incomplete strain group or a total strain given together with
        a split part are refused with a diagnostic.
    """
    skip = {str(x).strip().lower() for x in ignore}
    names: dict[str, str] = {}
    ignored: list[str] = []
    unknown: list[str] = []
    for column in columns:
        column = str(column).strip()
        key = column.lower()
        if key in skip:
            ignored.append(column)
            continue
        canonical = COLUMN_ALIASES.get(key)
        if canonical is None:
            unknown.append(column)
        elif canonical in names:
            raise ValueError(
                f"columns {names[canonical]!r} and {column!r} both mean {canonical!r}: give only one"
            )
        else:
            names[canonical] = column
    if unknown:
        raise ValueError(
            f"unknown column(s): {unknown}; rename them, or list them in ignore_columns to skip"
        )
    missing = [c for c in REQUIRED_COLUMNS if c not in names]
    if missing:
        raise ValueError(f"missing required column(s): {', '.join(missing)}")
    for prefix, what in _GROUPS.items():
        present = [f"{prefix}_{c}" for c in _COMPONENTS if f"{prefix}_{c}" in names]
        if present and len(present) != 6:
            lacking = [f"{prefix}_{c}" for c in _COMPONENTS if f"{prefix}_{c}" not in names]
            raise ValueError(
                f"incomplete {what} group {prefix}_*: missing {', '.join(lacking)} "
                "(all six components are needed)"
            )
    cmap = ColumnMap(names, tuple(ignored))
    if cmap.has_total and (cmap.has_elastic or cmap.has_plastic):
        raise ValueError(
            "give either the total strain E_* or the split Eel_*/Epl_*, not a total together "
            "with a split part"
        )
    return cmap


# --- importer ----------------------------------------------------------------------------------

_KEY_COLUMNS = ("node_id", "instance", "material", "step")
_COORDINATES = ("x", "y", "z")


def _dtypes(cmap: ColumnMap, canonical: Iterable[str]) -> dict[str, object]:
    out: dict[str, object] = {}
    for name in canonical:
        if name in ("node_id", "step"):
            out[cmap.names[name]] = "int64"
        elif name in ("instance", "material"):
            out[cmap.names[name]] = str
        else:
            out[cmap.names[name]] = "float64"
    return out


def _read_chunks(path, n_header_lines, cmap, canonical, chunk_rows):
    """Chunks of the data rows with the named canonical columns, blank numeric cells refused."""
    canonical = [c for c in canonical if c in cmap.names]
    reader = pd.read_csv(
        path,
        skiprows=n_header_lines,
        usecols=[cmap.names[c] for c in canonical],
        dtype=_dtypes(cmap, canonical),
        chunksize=chunk_rows,
        keep_default_na=False,
        na_values=[],
        encoding="utf-8-sig",
        float_precision="round_trip",
    )
    try:
        yield from reader
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{path}: cannot read the rows: {exc}") from exc


def _check_row_width(path, n_header_lines):
    """Reject ragged records before pandas can infer an implicit index or project columns."""
    with open(path, encoding="utf-8-sig") as handle:
        for _ in range(n_header_lines):
            next(handle)
        width = len(next(handle).rstrip("\r\n").split(","))
        for line_number, line in enumerate(handle, start=n_header_lines + 2):
            if line.strip() and len(line.rstrip("\r\n").split(",")) != width:
                raise ValueError(
                    f"{path}: line {line_number}: row column count differs from header ({width})"
                )


def _step_blocks(chunks, step_column):
    """Regroup chunks of a step-major file into one frame per run of equal step values."""
    parts: list[pd.DataFrame] = []
    current = None
    for chunk in chunks:
        steps = chunk[step_column].to_numpy()
        if not len(steps):
            continue
        cuts = np.flatnonzero(steps[1:] != steps[:-1]) + 1
        for piece in np.split(np.arange(len(chunk)), cuts):
            part = chunk.iloc[piece]
            value = int(steps[piece[0]])
            if current is not None and value != current:
                yield pd.concat(parts, ignore_index=True)
                parts = []
            current = value
            parts.append(part)
    if parts:
        yield pd.concat(parts, ignore_index=True)


def _keys(frame: pd.DataFrame, cmap: ColumnMap):
    n = len(frame)
    ids = frame[cmap.names["node_id"]].to_numpy(dtype=np.int64)
    inst = (
        frame[cmap.names["instance"]].to_numpy(dtype=str)
        if cmap.has_instance
        else np.full(n, "", dtype=str)
    )
    material = (
        frame[cmap.names["material"]].to_numpy(dtype=str)
        if cmap.has_material
        else np.full(n, "", dtype=str)
    )
    return inst, ids, material


def _positions(table: pd.MultiIndex, inst, ids, material) -> np.ndarray:
    """Position of each row's location in the location table (``-1`` where absent)."""
    return table.get_indexer(pd.MultiIndex.from_arrays([inst, ids, material]))


@dataclass
class _Scan:
    inst: np.ndarray
    ids: np.ndarray
    material: np.ndarray
    coordinates: np.ndarray  # mm
    table: pd.MultiIndex
    steps: list[int]
    times: list[float]


def _scan(blocks, cmap, header, path) -> _Scan:
    """First pass: validate the grid (steps ascending, same locations, keys, time) step by step."""
    scan = None
    previous_step = None
    for block in blocks:
        step = int(block[cmap.names["step"]].iloc[0])
        if previous_step is not None and step <= previous_step:
            raise ValueError(
                f"{path}: step {step} follows step {previous_step}: steps must ascend in a "
                "step-major file and each step must form one block of rows"
            )
        inst, ids, material = _keys(block, cmap)
        index = pd.MultiIndex.from_arrays([inst, ids, material])
        if not index.is_unique:
            dup = index[index.duplicated()][0]
            raise ValueError(
                f"{path}: step {step}: duplicate key (instance, node_id, material) {dup}"
            )
        xyz = to_mm(
            block[[cmap.names[c] for c in _COORDINATES]].to_numpy(float), header.length_unit
        )
        if not np.all(np.isfinite(xyz)):
            raise ValueError(f"{path}: step {step}: coordinates must be finite")
        if scan is None:
            scan = _Scan(inst, ids, material, xyz, index, [], [])
        else:
            pos = scan.table.get_indexer(index)
            if len(index) != len(scan.table) or (pos < 0).any():
                absent = scan.table.difference(index)
                new = index.difference(scan.table)
                raise ValueError(
                    f"{path}: step {step} has a different set of locations from the first step "
                    f"({len(absent)} missing, e.g. {absent[:1].tolist()}; {len(new)} new, e.g. "
                    f"{new[:1].tolist()}): the grid must be complete (a truncated or ragged file)"
                )
            aligned = np.empty_like(xyz)
            aligned[pos] = xyz
            extent = max(float(np.max(np.abs(scan.coordinates))), 1e-300)
            if np.max(np.abs(aligned - scan.coordinates)) > 1e-6 * extent:
                raise ValueError(
                    f"{path}: step {step}: coordinates differ from the first step; give the "
                    "undeformed coordinates of each location at every step"
                )
        if cmap.has_time:
            values = np.unique(block[cmap.names["time"]].to_numpy(float))
            if len(values) != 1 or not np.isfinite(values[0]):
                raise ValueError(
                    f"{path}: step {step}: the time column must hold one finite value per step"
                )
            if scan.times and not values[0] > scan.times[-1]:
                raise ValueError(
                    f"{path}: time {values[0]} at step {step} does not exceed {scan.times[-1]}: "
                    "time must increase strictly along ascending steps (the file is not re-sorted)"
                )
            scan.times.append(float(values[0]))
        scan.steps.append(step)
        previous_step = step
    if scan is None:
        raise ValueError(f"{path}: the file has a column line but no data rows")
    for key, got, what in (
        ("n_nodes", len(scan.table), "locations"),
        ("n_steps", len(scan.steps), "steps"),
    ):
        declared = getattr(header, key)
        if declared is not None and declared != got:
            raise ValueError(
                f"{path}: header declares {key}: {declared} but the file has {got} {what} "
                "(truncated or padded file)"
            )
    return scan


def _first_rows_are_one_location(path, n_header_lines, cmap) -> bool:
    keys = [c for c in ("instance", "node_id", "material") if c in cmap.names]
    head = next(_read_chunks(path, n_header_lines, cmap, keys, 2), None)
    if head is None or len(head) < 2:
        return False
    return bool((head.iloc[0] == head.iloc[1]).all())


def import_fe_csv(
    csv_path: str | Path,
    out_path: str | Path,
    *,
    elastic_constants: Mapping[str, tuple[float, float]] | None = None,
    allow_global_average: bool = False,
    ignore_columns: Sequence[str] = (),
    chunk_rows: int = 200_000,
    max_node_major_rows: int = 2_000_000,
) -> HistoryMeta:
    """Convert an FE CSV file to the universal history file, one load step at a time.

    Parameters
    ----------
    csv_path, out_path : str or Path
        The CSV file (see ``docs/fe-csv-format.md``) and the HDF5 file to write.
    elastic_constants : mapping, optional
        ``material -> (E, nu)`` with ``E`` in MPa; the key is the ``material`` column value, or the
        empty string when the file has none. Needed when the file gives the plastic strain
        without the elastic strain (the elastic part is then derived by Hooke's law from the
        stress and recorded as derived) and when elastic strain is to be derived for a
        ``linear_elastic`` file without strain columns.
    allow_global_average : bool
        Accept ``averaging: nodal_average_global`` with several materials or instances (the
        averaged values then mix materials at interface nodes); the choice is recorded.
    ignore_columns : sequence of str
        Exporter columns the format does not use (case-insensitive names).
    chunk_rows : int
        Rows per read buffer.
    max_node_major_rows : int
        A node-major file is sorted in memory; a larger one is refused with advice.

    Returns
    -------
    HistoryMeta
        The metadata written. The partial output is removed if anything fails.

    Notes
    -----
    A step-major file (all locations of one step, then the next) streams with the memory of one
    step. Steps must ascend, every step must hold the same locations, keys
    ``(instance, node_id, material, step)`` must be unique and ``time`` must increase strictly
    along ascending steps; nothing is re-sorted or repaired. The stress is converted to MPa, the
    lengths to mm and the strain to the tensor-shear convention.
    """
    csv_path, out_path = Path(csv_path), Path(out_path)
    header, n_header = read_header(csv_path)
    _check_row_width(csv_path, n_header)
    columns = pd.read_csv(csv_path, skiprows=n_header, nrows=0, encoding="utf-8-sig").columns
    cmap = resolve_columns(columns, ignore=ignore_columns)
    cmap.check_header(header)
    mode = _strain_mode(cmap, header, elastic_constants)
    wanted = list(CANONICAL_COLUMNS)

    if _first_rows_are_one_location(csv_path, n_header, cmap):
        frames = []
        n_rows = 0
        for chunk in _read_chunks(csv_path, n_header, cmap, wanted, chunk_rows):
            n_rows += len(chunk)
            if n_rows > max_node_major_rows:
                raise ValueError(
                    f"{csv_path}: this node-major file (rows grouped by location) has more than "
                    f"{max_node_major_rows} rows and would be sorted in memory; export it "
                    "step-major (all locations of one step, then the next step) to stream it"
                )
            frames.append(chunk)
        everything = pd.concat(frames, ignore_index=True)
        step_column = cmap.names["step"]

        def blocks(_columns):
            for _, group in everything.groupby(step_column, sort=True):
                yield group

    else:

        def blocks(columns_):
            chunks = _read_chunks(csv_path, n_header, cmap, columns_, chunk_rows)
            return _step_blocks(chunks, cmap.names["step"])

    scan = _scan(blocks(list(_KEY_COLUMNS) + ["time", *_COORDINATES]), cmap, header, csv_path)
    _check_averaging(header, scan, allow_global_average, csv_path)
    constants = _constants_for(mode, scan, elastic_constants)
    meta = _build_meta(csv_path, header, cmap, scan, mode, allow_global_average)
    steps = np.asarray(scan.steps, dtype=np.int64)
    cases = LoadCases(
        time=np.asarray(scan.times) if cmap.has_time else None,
        labels=tuple(f"step{s}" for s in scan.steps),
        steps=steps,
    )
    with atomic_destination(out_path, (csv_path,)) as temporary:
        with HistoryWriter(
            temporary,
            meta,
            scan.ids,
            scan.coordinates,
            cases,
            len(steps),
            with_strain=mode != "none",
            instance=scan.inst if cmap.has_instance else None,
            material=scan.material if cmap.has_material else None,
            with_strain_split=mode in ("split", "plastic_only"),
            with_temperature=cmap.has_temperature,
        ) as writer:
            for k, block in enumerate(blocks(wanted)):
                _write_block(writer, k, block, cmap, header, scan, mode, constants)
    return meta


def _strain_mode(cmap, header, elastic_constants) -> str:
    if cmap.has_total:
        return "total"
    if cmap.has_split:
        return "split"
    if cmap.has_elastic_only:
        if header.analysis_type == "elastic_plastic":
            raise ValueError(
                "an elastic_plastic file with only Eel columns is refused: the elastic strain is "
                "not the total strain under plasticity; give Epl as well, or the total E"
            )
        return "elastic_only"
    if cmap.has_plastic_only:
        if not elastic_constants:
            raise ValueError(
                "the file gives the plastic strain Epl without the elastic strain: pass "
                "elastic_constants={material: (E, nu)} so the elastic part can be derived, or "
                "export Eel as well"
            )
        return "plastic_only"
    # no strain columns (elastic-plastic was refused by check_header)
    return "derived" if elastic_constants else "none"


def _constants_for(mode, scan, elastic_constants) -> dict[str, tuple[float, float]]:
    if mode not in ("plastic_only", "derived"):
        return {}
    materials = sorted(set(scan.material.tolist()))
    missing = [m for m in materials if m not in elastic_constants]
    if missing:
        raise ValueError(
            f"elastic_constants lacks the material(s) {missing} of the file (use the key '' for "
            "a file without a material column)"
        )
    constants = {}
    for material in materials:
        E, nu = elastic_constants[material]
        elastic_strain(np.zeros((1, 3, 3)), E, nu)  # validates E and nu
        constants[material] = (float(E), float(nu))
    return constants


def _check_averaging(header, scan, allow, path) -> None:
    if header.averaging != "nodal_average_global" or allow:
        return
    n_materials = len(set(scan.material.tolist()))
    n_instances = len(set(scan.inst.tolist()))
    if n_materials > 1 or n_instances > 1:
        raise ValueError(
            f"{path}: averaging is nodal_average_global but the file has {n_materials} material(s) "
            f"and {n_instances} instance(s): averaged values at interface nodes mix materials. "
            "Export averaged within each material, or pass allow_global_average=True to accept it"
        )


def _build_meta(path, header, cmap, scan, mode, allow_global_average) -> HistoryMeta:
    materials = sorted({m for m in scan.material.tolist() if m})
    mapping = {c: a for c, a in cmap.names.items() if c != a}
    conventions = {
        "component_order": ",".join(CANONICAL),
        "column_aliases_used": json.dumps(mapping, sort_keys=True),
    }
    if header.strain_convention:
        conventions["strain_convention"] = header.strain_convention
    if header.averaging == "nodal_average_global":
        conventions["global_average_allowed"] = str(bool(allow_global_average)).lower()
    units = {"stress": header.stress_unit, "length": header.length_unit}
    for key, value in (
        ("strain", header.strain_unit),
        ("time", header.time_unit),
        ("temperature", header.temperature_unit),
    ):
        if value:
            units[key] = value
    kind = {
        "total": "total_mechanical",
        "split": "total_mechanical",
        "elastic_only": "total_mechanical",
        "plastic_only": "derived_elastic_part",
        "derived": "derived_linear_elastic",
        "none": None,
    }[mode]
    return HistoryMeta(
        material_id="+".join(materials),
        solver=header.solver,
        solver_version=header.solver_version,
        analysis_step=f"steps {scan.steps[0]}..{scan.steps[-1]}",
        adapter="multifat.fe.csvio",
        adapter_version="1.0",
        source_files=(SourceFile(Path(path).name, sha256_of(path)),),
        source_units=units,
        source_frame="global_cartesian",
        averaging=header.averaging,
        location_kind="node",
        created=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        notes=header.notes,
        strain_kind=kind,
        source_conventions=conventions,
        analysis_type=header.analysis_type,
        cycle=header.cycle or "unspecified",
        step_role=header.step_role,
        load_description=header.load_description,
    )


def _tensors(block, cmap, group, pos, n, *, kind, shear=None, unit=None):
    names = [cmap.names[f"{group}_{c}"] for c in _COMPONENTS]
    values = block[names].to_numpy(float)
    if kind == "stress":
        T = to_mpa(components_to_sym3(values, CANONICAL, kind="stress"), unit)
    else:
        T = to_strain(components_to_sym3(values, CANONICAL, kind="strain", shear=shear), unit)
    out = np.empty((n, 3, 3))
    out[pos] = T
    return out


def _write_block(writer, k, block, cmap, header, scan, mode, constants) -> None:
    n = len(scan.table)
    inst, ids, material = _keys(block, cmap)
    pos = _positions(scan.table, inst, ids, material)
    S = _tensors(block, cmap, "S", pos, n, kind="stress", unit=header.stress_unit)

    def strain(group):
        return _tensors(
            block, cmap, group, pos, n, kind="strain", shear=header.strain_convention,
            unit=header.strain_unit,
        )  # fmt: skip

    def hooke():
        out = np.empty((n, 3, 3))
        for material, (E, nu) in constants.items():
            mask = scan.material == material
            out[mask] = elastic_strain(S[mask], E, nu)
        return out

    E_total = Eel = Epl = None
    if mode == "total":
        E_total = strain("E")
    elif mode == "split":
        Eel, Epl = strain("Eel"), strain("Epl")
        E_total = Eel + Epl
    elif mode == "elastic_only":
        E_total = strain("Eel")
    elif mode == "plastic_only":
        Epl = strain("Epl")
        Eel = hooke()
        E_total = Eel + Epl
    elif mode == "derived":
        E_total = hooke()
    T = None
    if cmap.has_temperature:
        T = np.empty(n)
        T[pos] = block[cmap.names["T"]].to_numpy(float)
    writer.write_step(k, S, E_total, E_elastic=Eel, E_plastic=Epl, T=T)
