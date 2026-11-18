"""Data model of the FE interface: history metadata, load-case axis and blocks of locations.

An FE history is a prediction, never an experimental observation: its types are separate from
the workbook dictionaries and no calibration or comparison entry point accepts them.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

AVERAGING = (
    "none",
    "nodal_average",
    "extrapolated",
    "nodal_average_within_material",
    "nodal_average_global",
)
LOCATION_KINDS = ("node", "element_centroid", "integration_point", "element_nodal")
ANALYSIS_TYPES = ("unspecified", "linear_elastic", "elastic_plastic")
STEP_ROLES = ("time_sample", "unit_case")
STRAIN_KINDS = ("total_mechanical", "derived_linear_elastic", "derived_elastic_part")
_SYM_RTOL = 1e-9


@dataclass(frozen=True)
class SourceFile:
    """A file an adapter read, identified by name and SHA-256."""

    name: str
    sha256: str


def sha256_of(path: str | Path) -> str:
    """SHA-256 hex digest of a file, read in blocks."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class LoadCases:
    """Optional coordinates of the ordered load cases: time, solver step numbers and labels."""

    time: np.ndarray | None = None
    labels: tuple[str, ...] | None = None
    steps: np.ndarray | None = None

    def __post_init__(self):
        if self.time is not None:
            time = np.asarray(self.time, dtype=float)
            if time.ndim != 1:
                raise ValueError(f"load-case time must be 1-D, got shape {time.shape}")
            object.__setattr__(self, "time", time)
        if self.labels is not None:
            object.__setattr__(self, "labels", tuple(str(x) for x in self.labels))
        if self.steps is not None:
            steps = np.asarray(self.steps)
            if steps.ndim != 1 or steps.dtype.kind not in "iu":
                raise ValueError(f"load-case steps must be a 1-D integer array, got {steps.shape}")
            object.__setattr__(self, "steps", steps)
        lengths = {
            len(x) for x in (self.time, self.labels, self.steps) if x is not None
        }  # fmt: skip
        if len(lengths) > 1:
            raise ValueError("load-case time, labels and steps differ in length")


@dataclass(frozen=True)
class HistoryMeta:
    """Provenance and conventions of a history file (format version 1.1)."""

    material_id: str
    solver: str
    solver_version: str
    analysis_step: str
    adapter: str
    adapter_version: str
    source_files: tuple[SourceFile, ...]
    source_units: dict[str, str]
    source_frame: str
    averaging: str
    location_kind: str
    periodic: bool = False
    local_only: bool = True
    strain_kind: str | None = "total_mechanical"
    kind: str = "fe_prediction"
    format_version: str = "1.1"
    created: str = ""
    notes: str = ""
    source_conventions: dict[str, str] = field(default_factory=dict)
    analysis_type: str = "unspecified"
    cycle: str = "unspecified"
    step_role: str = "time_sample"
    load_description: str = ""

    def __post_init__(self):
        if self.averaging not in AVERAGING:
            raise ValueError(f"averaging must be one of {AVERAGING}, got {self.averaging!r}")
        if self.location_kind not in LOCATION_KINDS:
            raise ValueError(
                f"location_kind must be one of {LOCATION_KINDS}, got {self.location_kind!r}"
            )
        if self.analysis_type not in ANALYSIS_TYPES:
            raise ValueError(
                f"analysis_type must be one of {ANALYSIS_TYPES}, got {self.analysis_type!r}"
            )
        if self.step_role not in STEP_ROLES:
            raise ValueError(f"step_role must be one of {STEP_ROLES}, got {self.step_role!r}")
        if self.strain_kind is not None and self.strain_kind not in STRAIN_KINDS:
            raise ValueError(f"strain_kind must be one of {STRAIN_KINDS}, got {self.strain_kind!r}")
        object.__setattr__(self, "source_files", tuple(self.source_files))

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe dictionary."""
        out = {k: getattr(self, k) for k in self.__dataclass_fields__}
        out["source_files"] = [{"name": s.name, "sha256": s.sha256} for s in self.source_files]
        out["source_units"] = dict(self.source_units)
        return out

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "HistoryMeta":
        """Inverse of :meth:`to_dict`."""
        data = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}  # newer 1.x fields
        data["source_files"] = tuple(SourceFile(**s) for s in data["source_files"])
        return cls(**data)


def location_keys(
    location_id: np.ndarray, instance: np.ndarray | None, material: np.ndarray | None
) -> list[tuple]:
    """Identity of each location: the triple ``(instance, location_id, material)``.

    A missing ``instance`` or ``material`` counts as the empty string, so a file without them
    keys on the location ID alone.
    """
    ids = np.asarray(location_id).tolist()
    inst = [""] * len(ids) if instance is None else [str(x) for x in instance]
    names = [""] * len(ids) if material is None else [str(x) for x in material]
    return list(zip(inst, ids, names))


def has_duplicate_locations(keys: list[tuple]) -> bool:
    """True when two locations share the same identity."""
    return len(set(keys)) != len(keys)


@dataclass(frozen=True, eq=False)
class HistoryBlock:
    """Histories of ``M`` locations over ``N`` ordered load cases (MPa, mm, tensor strain).

    ``instance`` and ``material`` are optional ``(M,)`` string arrays; ``E_elastic_hist`` and
    ``E_plastic_hist`` the optional split of the strain and ``T_hist`` the optional ``(M, N)``
    temperature.
    """

    location_id: np.ndarray
    coordinates: np.ndarray
    S_hist: np.ndarray
    E_hist: np.ndarray | None = field(default=None)
    instance: np.ndarray | None = field(default=None)
    material: np.ndarray | None = field(default=None)
    E_elastic_hist: np.ndarray | None = field(default=None)
    E_plastic_hist: np.ndarray | None = field(default=None)
    T_hist: np.ndarray | None = field(default=None)

    def __post_init__(self):
        ids = np.asarray(self.location_id)
        xyz = np.asarray(self.coordinates, dtype=float)
        S = np.asarray(self.S_hist, dtype=float)
        E = None if self.E_hist is None else np.asarray(self.E_hist, dtype=float)
        Eel = None if self.E_elastic_hist is None else np.asarray(self.E_elastic_hist, dtype=float)
        Epl = None if self.E_plastic_hist is None else np.asarray(self.E_plastic_hist, dtype=float)
        T = None if self.T_hist is None else np.asarray(self.T_hist, dtype=float)
        if ids.ndim != 1:
            raise ValueError(f"location_id must be 1-D, got shape {ids.shape}")
        m = ids.shape[0]
        if xyz.shape != (m, 3):
            raise ValueError(f"coordinates must have shape ({m}, 3), got {xyz.shape}")
        if S.ndim != 4 or S.shape[0] != m or S.shape[2:] != (3, 3):
            raise ValueError(f"S_hist must have shape ({m}, N, 3, 3), got {S.shape}")
        if S.shape[1] < 1:
            raise ValueError("a history needs at least one load case")
        for name, arr in (("E_hist", E), ("E_elastic_hist", Eel), ("E_plastic_hist", Epl)):
            if arr is not None and arr.shape != S.shape:
                raise ValueError(f"{name} shape {arr.shape} differs from S_hist shape {S.shape}")
        if T is not None and T.shape != S.shape[:2]:
            raise ValueError(f"T_hist must have shape {S.shape[:2]}, got {T.shape}")
        labels = {}
        for name in ("instance", "material"):
            value = getattr(self, name)
            if value is not None:
                value = np.asarray(value).astype(str)
                if value.shape != (m,):
                    raise ValueError(f"{name} must have shape ({m},), got {value.shape}")
            labels[name] = value
        if has_duplicate_locations(location_keys(ids, labels["instance"], labels["material"])):
            raise ValueError("duplicate location IDs")
        for name, arr in (
            ("coordinates", xyz),
            ("S_hist", S),
            ("E_hist", E),
            ("E_elastic_hist", Eel),
            ("E_plastic_hist", Epl),
            ("T_hist", T),
        ):
            if arr is not None and not np.all(np.isfinite(arr)):
                raise ValueError(f"{name} must be finite")
        for name, arr in (
            ("S_hist", S),
            ("E_hist", E),
            ("E_elastic_hist", Eel),
            ("E_plastic_hist", Epl),
        ):
            if arr is None:
                continue
            scale = float(np.max(np.abs(arr))) if arr.size else 0.0
            if np.max(np.abs(arr - np.swapaxes(arr, -1, -2)), initial=0.0) > _SYM_RTOL * scale:
                raise ValueError(f"{name} tensors must be symmetric")
        object.__setattr__(self, "location_id", ids)
        object.__setattr__(self, "coordinates", xyz)
        object.__setattr__(self, "S_hist", S)
        object.__setattr__(self, "E_hist", E)
        object.__setattr__(self, "E_elastic_hist", Eel)
        object.__setattr__(self, "E_plastic_hist", Epl)
        object.__setattr__(self, "T_hist", T)
        object.__setattr__(self, "instance", labels["instance"])
        object.__setattr__(self, "material", labels["material"])

    @property
    def n_locations(self) -> int:
        return int(self.location_id.shape[0])

    @property
    def n_load_cases(self) -> int:
        return int(self.S_hist.shape[1])
