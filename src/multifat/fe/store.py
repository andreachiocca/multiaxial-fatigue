"""The universal FE history file: a streaming HDF5 writer and a chunked reader.

Layout (format version 1.1; a 1.0 file opens unchanged): root attributes ``format``, ``format_version`` and ``complete`` (true
only once every load case is written, so a crashed export is refused, never read as zeros); ``/meta`` (JSON
attribute); ``/locations/id`` and ``/locations/coordinates``; ``/load_cases/index`` with optional
``time``, ``step`` and ``labels``; optional ``/locations/instance`` and ``/locations/material``
(UTF-8 strings; a location is the triple instance, ID, material); ``/fields/stress`` and optional
``/fields/strain``, ``/fields/strain_elastic``, ``/fields/strain_plastic`` of shape ``(P, N, 6)``
in the canonical component order, chunked along locations so one location's history is
contiguous, and optional ``/fields/temperature`` of shape ``(P, N)``; an empty reserved group
``/mesh``. ``h5py`` is an optional dependency (extra ``fe``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator

import numpy as np

from .model import (
    HistoryBlock,
    HistoryMeta,
    LoadCases,
    has_duplicate_locations,
    location_keys,
)
from .tensors import six_to_sym3, sym3_to_six

FORMAT_NAME = "multifat-fe-history"
FORMAT_VERSION = "1.1"
_SYM_RTOL = 1e-9
CHUNK_BYTES = 4 * 2**20  # target size of one HDF5 chunk (the format limit is 4 GiB)


def require_h5py():
    """Import and return ``h5py``, or explain which extra to install."""
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "FE history files need h5py: install it with the optional extra, multifat[fe]"
        ) from exc
    return h5py


def _check_slab(name: str, T: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    if T.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {T.shape}")
    if not np.all(np.isfinite(T)):
        raise ValueError(f"{name} must be finite")
    scale = float(np.max(np.abs(T))) if T.size else 0.0
    if np.max(np.abs(T - np.swapaxes(T, -1, -2)), initial=0.0) > _SYM_RTOL * scale:
        raise ValueError(f"{name} tensors must be symmetric")
    return T


def _check_labels(name: str, values, count: int) -> np.ndarray | None:
    if values is None:
        return None
    arr = np.asarray(values).astype(str)
    if arr.shape != (count,):
        raise ValueError(f"{name} must have shape ({count},), got {arr.shape}")
    return arr


def _string_dataset(group, name: str, values: np.ndarray) -> None:
    h5py = require_h5py()
    group.create_dataset(name, data=np.asarray(values, dtype=object), dtype=h5py.string_dtype())


def _read_strings(group, name: str) -> np.ndarray | None:
    if name not in group:
        return None
    return np.asarray(group[name].asstr()[...], dtype=str)


class HistoryWriter:
    """Write a history file one load step at a time (memory is one step, never ``P x N``)."""

    def __init__(
        self,
        path: str | Path,
        meta: HistoryMeta,
        location_id: np.ndarray,
        coordinates: np.ndarray,
        load_cases: LoadCases,
        n_load_cases: int,
        *,
        with_strain: bool,
        dtype: str = "float64",
        chunk_locations: int | None = None,
        instance: np.ndarray | None = None,
        material: np.ndarray | None = None,
        with_strain_split: bool = False,
        with_temperature: bool = False,
    ):
        h5py = require_h5py()
        ids = np.asarray(location_id)
        if ids.ndim != 1 or ids.dtype.kind not in "iu":
            raise TypeError("location IDs must be a 1-D integer array")
        instance = _check_labels("instance", instance, len(ids))
        material = _check_labels("material", material, len(ids))
        if has_duplicate_locations(location_keys(ids, instance, material)):
            raise ValueError("duplicate location IDs")
        xyz = np.asarray(coordinates, dtype=float)
        if xyz.shape != (len(ids), 3):
            raise ValueError(f"coordinates must have shape ({len(ids)}, 3), got {xyz.shape}")
        if len(ids) < 1 or n_load_cases < 1:
            raise ValueError("a history needs at least one location and one load case")
        if dtype not in ("float64", "float32"):
            raise ValueError(f"dtype must be 'float64' or 'float32', got {dtype!r}")
        if load_cases.time is not None and len(load_cases.time) != n_load_cases:
            raise ValueError("load-case time does not match the number of load cases")
        if load_cases.steps is not None and len(load_cases.steps) != n_load_cases:
            raise ValueError("load-case steps do not match the number of load cases")
        self._n, self._p, self._with_strain = n_load_cases, len(ids), with_strain
        self._split, self._temperature = with_strain_split, with_temperature
        self._written: set[int] = set()
        self._f = f = h5py.File(path, "w")
        f.attrs["format"] = FORMAT_NAME
        f.attrs["format_version"] = FORMAT_VERSION
        f.attrs["complete"] = False  # set only when every load case has been written
        f.create_group("meta").attrs["json"] = json.dumps(meta.to_dict())
        loc = f.create_group("locations")
        loc.create_dataset("id", data=ids)
        loc.create_dataset("coordinates", data=xyz)
        if instance is not None:
            _string_dataset(loc, "instance", instance)
        if material is not None:
            _string_dataset(loc, "material", material)
        lc = f.create_group("load_cases")
        lc.create_dataset("index", data=np.arange(n_load_cases))
        if load_cases.time is not None:
            lc.create_dataset("time", data=load_cases.time)
        if load_cases.steps is not None:
            lc.create_dataset("step", data=load_cases.steps)
        if load_cases.labels is not None:
            lc.attrs["labels"] = json.dumps(list(load_cases.labels))
        f.create_group("mesh")
        fields = f.create_group("fields")
        itemsize = np.dtype(dtype).itemsize
        if (
            chunk_locations is None
        ):  # about CHUNK_BYTES per chunk, whatever the number of load cases
            chunk_locations = min(4096, max(1, CHUNK_BYTES // (n_load_cases * 6 * itemsize)))
        chunks = (min(chunk_locations, self._p), n_load_cases, 6)
        names = ["stress"] + (["strain"] if with_strain else [])
        if with_strain_split:
            names += ["strain_elastic", "strain_plastic"]
        for name in names:
            fields.create_dataset(
                name, shape=(self._p, n_load_cases, 6), dtype=dtype, chunks=chunks
            )
        if with_temperature:
            fields.create_dataset(
                "temperature",
                shape=(self._p, n_load_cases),
                dtype=dtype,
                chunks=(chunks[0], n_load_cases),
            )

    def write_step(
        self,
        step: int,
        S: np.ndarray,
        E: np.ndarray | None,
        *,
        E_elastic: np.ndarray | None = None,
        E_plastic: np.ndarray | None = None,
        T: np.ndarray | None = None,
    ) -> None:
        """Store the tensors of all locations at one load case."""
        if not 0 <= step < self._n:
            raise ValueError(f"step {step} outside 0..{self._n - 1}")
        if self._with_strain != (E is not None):
            raise ValueError("strain must be given at every step if and only if with_strain=True")
        if self._split != (E_elastic is not None) or self._split != (E_plastic is not None):
            raise ValueError(
                "elastic and plastic strain must both be given at every step if and only if "
                "with_strain_split=True (strain split)"
            )
        if self._temperature != (T is not None):
            raise ValueError(
                "temperature must be given at every step if and only if with_temperature=True"
            )
        shape = (self._p, 3, 3)
        self._f["fields/stress"][:, step, :] = sym3_to_six(_check_slab("stress", S, shape))
        if E is not None:
            self._f["fields/strain"][:, step, :] = sym3_to_six(_check_slab("strain", E, shape))
        if self._split:
            for name, arr in (("strain_elastic", E_elastic), ("strain_plastic", E_plastic)):
                self._f[f"fields/{name}"][:, step, :] = sym3_to_six(_check_slab(name, arr, shape))
        if T is not None:
            T = np.asarray(T, dtype=float)
            if T.shape != (self._p,) or not np.all(np.isfinite(T)):
                raise ValueError(f"temperature must be finite with shape ({self._p},)")
            self._f["fields/temperature"][:, step] = T
        self._written.add(step)

    def abort(self) -> None:
        """Close the file without the completeness check."""
        if self._f:
            self._f.close()

    def __enter__(self) -> "HistoryWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        missing = sorted(set(range(self._n)) - self._written)
        if exc_type is None and not missing:
            self._f.attrs["complete"] = True
        self.abort()
        if exc_type is None and missing:
            raise ValueError(f"load case(s) never written: {missing[:10]}")


class HistoryFile:
    """Read a history file lazily: metadata and locations at once, tensors by chunk."""

    def __init__(self, path: str | Path):
        h5py = require_h5py()
        self._f = f = h5py.File(path, "r")
        try:
            if f.attrs.get("format") != FORMAT_NAME:
                raise ValueError(f"{path} is not a {FORMAT_NAME} file")
            if not f.attrs.get("complete", False):
                raise ValueError(f"{path} is incomplete: its writer did not finish")
            version = str(f.attrs["format_version"])
            if version.split(".")[0] != FORMAT_VERSION.split(".")[0]:
                raise ValueError(f"unsupported format version {version}")
            self.meta = HistoryMeta.from_dict(json.loads(f["meta"].attrs["json"]))
            self.location_id = f["locations/id"][...]
            self.coordinates = f["locations/coordinates"][...]
            self.instance = _read_strings(f["locations"], "instance")
            self.material = _read_strings(f["locations"], "material")
            if has_duplicate_locations(
                location_keys(self.location_id, self.instance, self.material)
            ):
                raise ValueError(f"{path}: duplicate location IDs")
            lc = f["load_cases"]
            labels = lc.attrs.get("labels")
            self.load_cases = LoadCases(
                time=lc["time"][...] if "time" in lc else None,
                labels=None if labels is None else tuple(json.loads(labels)),
                steps=lc["step"][...] if "step" in lc else None,
            )
            fields = f["fields"]
            self._stress = fields["stress"]
            self._strain = fields["strain"] if "strain" in fields else None
            self._strain_el = fields["strain_elastic"] if "strain_elastic" in fields else None
            self._strain_pl = fields["strain_plastic"] if "strain_plastic" in fields else None
            self._temperature = fields["temperature"] if "temperature" in fields else None
        except Exception:
            f.close()
            raise
        self.n_locations, self.n_load_cases = self._stress.shape[:2]
        self.has_strain = self._strain is not None
        self.has_strain_split = self._strain_el is not None and self._strain_pl is not None
        self.has_temperature = self._temperature is not None
        chunks = self._stress.chunks
        self._chunk = chunks[0] if chunks else min(4096, self.n_locations)

    def close(self) -> None:
        self._f.close()

    def __enter__(self) -> "HistoryFile":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _read(self, field, index) -> np.ndarray:
        if isinstance(index, slice):
            return field[index]
        order = np.argsort(index)
        out = np.empty((len(index),) + field.shape[1:], dtype=field.dtype)
        if len(index):
            out[order] = field[np.asarray(index)[order]]
        return out

    def read(self, index) -> HistoryBlock:
        """Read the locations at ``index`` (a slice or an array of distinct positions)."""
        if isinstance(index, slice):
            positions = np.arange(self.n_locations)[index]
        else:
            positions = np.asarray(index, dtype=int)
            if positions.ndim != 1:
                raise ValueError("index must be 1-D")
            if len(positions) and (positions.min() < 0 or positions.max() >= self.n_locations):
                raise IndexError("location position out of range")
            if len(np.unique(positions)) != len(positions):
                raise ValueError("index contains duplicate positions")
            index = positions

        def tensors(field):
            return None if field is None else six_to_sym3(self._read(field, index).astype(float))

        T = None
        if self._temperature is not None:
            T = self._read(self._temperature, index).astype(float)
        return HistoryBlock(
            self.location_id[positions],
            self.coordinates[positions],
            tensors(self._stress),
            tensors(self._strain),
            None if self.instance is None else self.instance[positions],
            None if self.material is None else self.material[positions],
            tensors(self._strain_el),
            tensors(self._strain_pl),
            T,
        )

    def iter_chunks(
        self,
        *,
        max_locations: int | None = None,
        max_bytes: int | None = None,
        index: np.ndarray | None = None,
    ) -> Iterator[HistoryBlock]:
        """Yield blocks of locations, each with all load cases, within the given budget.

        ``max_bytes`` is an estimate, not a hard cap: it counts the stored six components, the
        expanded tensors and a copy for validation (a factor of three on the stored size) for
        stress and strain, and the numpy peak stays below three times it (tested). HDF5's own
        buffers are not counted.
        """
        n_tensor_fields = 1 + self.has_strain + 2 * self.has_strain_split
        per_location = self.n_load_cases * 6 * 8 * n_tensor_fields * 3
        sizes = [self._chunk]
        if max_locations is not None:
            sizes = [max_locations]
        if max_bytes is not None:
            sizes.append(max_bytes // per_location)
        size = max(1, min(sizes))
        if index is None:
            for start in range(0, self.n_locations, size):
                yield self.read(slice(start, min(start + size, self.n_locations)))
            return
        positions = np.asarray(index, dtype=int)
        for start in range(0, len(positions), size):
            yield self.read(positions[start : start + size])
