"""Time-history input contract: synchronised tensors over a time or cycle axis, with provenance.

A ``TimeHistory`` is one material point: ``N >= 2`` stress tensors (MPa, tensor components) and
optionally strain tensors (dimensionless, tensor shear), on a strictly increasing time or cycle axis.
It is built from arrays, a CSV file, a one-location history file of ``multifat.fe`` or a documented
Python function, and ``record`` writes the sampled history that an assessment used, so that a result
can be reproduced from it. Units are always stated by the caller, never assumed.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from multifat.fe.model import HistoryBlock, HistoryMeta, sha256_of
from multifat.fe.tensors import CANONICAL, components_to_sym3, to_mpa

AXIS_KINDS = ("time", "cycle")
_SYM_RTOL = 1e-9
DEFAULT_COLUMNS = {
    "time": "t",
    **{comp: "s" + comp for comp in CANONICAL},
    **{"e" + comp: "e" + comp for comp in CANONICAL},
}
# Further names accepted when the caller does not rename a column: those of the FE node CSV
# (docs/fe-csv-format.md), so one file can feed both the node importer and a time history.
_ALIAS_COLUMNS = {
    "time": ("time",),
    **{comp: ("S_" + comp,) for comp in CANONICAL},
    **{"e" + comp: ("E_" + comp,) for comp in CANONICAL},
}


def _json_safe(name: str, value: Any) -> None:
    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be JSON-serialisable: {exc}") from exc


def _metadata_hash(axis_kind: str, units: Mapping[str, str], provenance: Mapping[str, Any]) -> str:
    """SHA-256 of the metadata that the array hash does not cover (canonical JSON)."""
    text = json.dumps(
        {"axis_kind": axis_kind, "units": dict(units), "provenance": dict(provenance)},
        sort_keys=True,
        allow_nan=False,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _check_time_unit(time_unit: str | None) -> None:
    if time_unit is not None and not (isinstance(time_unit, str) and time_unit.strip()):
        raise ValueError("time_unit must be a non-empty string (for example 's' or 'cycles')")


def _check_tensors(name: str, T: np.ndarray, n: int) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    if T.shape != (n, 3, 3):
        raise ValueError(f"{name} must have shape ({n}, 3, 3), got {T.shape}")
    if not np.all(np.isfinite(T)):
        raise ValueError(f"{name} must be finite")
    scale = float(np.max(np.abs(T)))
    if np.max(np.abs(T - np.swapaxes(T, -1, -2))) > _SYM_RTOL * scale:
        raise ValueError(f"{name} tensors must be symmetric")
    return T


@dataclass(frozen=True, eq=False)
class TimeHistory:
    """Stress (and optional strain) tensors of one material point over an ordered axis.

    Parameters
    ----------
    time : ndarray
        Strictly increasing, finite axis of length ``N >= 2`` (time or cycle count, see ``axis_kind``).
    S_hist : ndarray
        Symmetric stress tensors in MPa, shape ``(N, 3, 3)``.
    E_hist : ndarray or None
        Total strain tensors (tensor shear), shape ``(N, 3, 3)``.
    axis_kind : {"time", "cycle"}
        What ``time`` measures.
    units : dict
        Must hold ``"stress"``; the constructors write ``"MPa"`` and the source unit.
    provenance : dict
        JSON-serialisable record of where the history came from.
    """

    time: np.ndarray
    S_hist: np.ndarray
    E_hist: np.ndarray | None = None
    axis_kind: str = "time"
    units: dict[str, str] = field(default_factory=lambda: {"stress": "MPa"})
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        time = np.asarray(self.time, dtype=float)
        if time.ndim != 1:
            raise ValueError(f"time must be 1-D, got shape {time.shape}")
        n = time.shape[0]
        if n < 2:
            raise ValueError("a time history needs at least two samples")
        if not np.all(np.isfinite(time)):
            raise ValueError("time must be finite")
        if np.any(np.diff(time) <= 0.0):
            raise ValueError("time must be strictly increasing (no duplicates, no reversals)")
        S = _check_tensors("S_hist", self.S_hist, n)
        E = None if self.E_hist is None else _check_tensors("E_hist", self.E_hist, n)
        if self.axis_kind not in AXIS_KINDS:
            raise ValueError(f"axis_kind must be one of {AXIS_KINDS}, got {self.axis_kind!r}")
        if "stress" not in self.units:
            raise ValueError("units must state the stress unit")
        _json_safe("provenance", self.provenance)
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "S_hist", S)
        object.__setattr__(self, "E_hist", E)

    @property
    def n_samples(self) -> int:
        return int(self.time.shape[0])

    def sha256(self) -> str:
        """SHA-256 of the sampled history: axis, stress and strain as little-endian float64."""
        digest = hashlib.sha256()
        for array in (self.time, self.S_hist):
            digest.update(np.ascontiguousarray(array, dtype="<f8").tobytes())
        digest.update(b"E" if self.E_hist is not None else b"-")
        if self.E_hist is not None:
            digest.update(np.ascontiguousarray(self.E_hist, dtype="<f8").tobytes())
        return digest.hexdigest()

    def to_block(self, location_id: int = 0) -> HistoryBlock:
        """The history as a one-location ``multifat.fe`` block (axis dropped)."""
        E = None if self.E_hist is None else self.E_hist[None]
        return HistoryBlock(np.array([location_id]), np.zeros((1, 3)), self.S_hist[None], E)

    def meta(self) -> HistoryMeta:
        """Metadata for the ``multifat.fe`` bridge: local-only, MPa, global axes."""
        return HistoryMeta(
            material_id=str(self.provenance.get("material_id", "unspecified")),
            solver=str(self.provenance.get("source", "time history")),
            solver_version="",
            analysis_step="",
            adapter="multifat.history",
            adapter_version="1.0",
            source_files=(),
            source_units={"stress": "MPa", "length": "mm"},
            source_frame="global",
            averaging="none",
            location_kind="node",
            strain_kind="total_mechanical" if self.E_hist is not None else None,
        )

    def record(self, path: str | Path) -> Path:
        """Write the sampled history (``<path>.npz``) and its metadata (``<path>.json``).

        Returns the ``.json`` path. ``from_record`` reads both and checks the hash.
        """
        base = Path(path)
        arrays = {"time": self.time, "S_hist": self.S_hist}
        if self.E_hist is not None:
            arrays["E_hist"] = self.E_hist
        np.savez(base.parent / (base.name + ".npz"), **arrays)
        info = {
            "format": "multifat-time-history",
            "axis_kind": self.axis_kind,
            "units": self.units,
            "provenance": self.provenance,
            "sha256": self.sha256(),
            "metadata_sha256": _metadata_hash(self.axis_kind, self.units, self.provenance),
        }
        out = base.parent / (base.name + ".json")
        out.write_text(json.dumps(info, indent=2, sort_keys=True), encoding="utf-8")
        return out

    @classmethod
    def from_record(cls, path: str | Path) -> "TimeHistory":
        """Read a history written by ``record`` and verify its SHA-256."""
        base = Path(path)
        if base.suffix in (".json", ".npz"):
            base = base.with_suffix("")
        info = json.loads((base.parent / (base.name + ".json")).read_text(encoding="utf-8"))
        with np.load(base.parent / (base.name + ".npz")) as data:
            E = data["E_hist"] if "E_hist" in data.files else None
            history = cls(
                data["time"],
                data["S_hist"],
                E,
                info["axis_kind"],
                info["units"],
                info["provenance"],
            )
        if history.sha256() != info["sha256"]:
            raise ValueError(f"{base}: the arrays do not match the recorded hash")
        recorded = info.get("metadata_sha256")
        if (
            recorded is not None
            and _metadata_hash(history.axis_kind, history.units, history.provenance) != recorded
        ):
            raise ValueError(
                f"{base}: the metadata (units, axis, provenance) do not match the recorded hash"
            )
        return history

    @classmethod
    def from_arrays(
        cls,
        time,
        S_hist,
        E_hist=None,
        *,
        stress_unit: str,
        axis_kind: str = "time",
        time_unit: str | None = None,
        provenance: Mapping[str, Any] | None = None,
    ) -> "TimeHistory":
        """Build from arrays; ``stress_unit`` (Pa, kPa, MPa, GPa) is required and converted to MPa.

        ``time_unit`` (for example ``"s"`` or ``"cycles"``) is recorded when given. Strain, if any, is
        taken as dimensionless tensor strain: convert microstrain or percent before building.
        """
        _check_time_unit(time_unit)
        S = to_mpa(np.asarray(S_hist, dtype=float), stress_unit)
        units = {"stress": "MPa", "source_stress_unit": stress_unit}
        if time_unit is not None:
            units["time"] = time_unit
        return cls(
            np.asarray(time, dtype=float), S, E_hist, axis_kind, units, dict(provenance or {})
        )

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        stress_unit: str,
        time_unit: str,
        axis_kind: str = "time",
        delimiter: str = ",",
        columns: Mapping[str, str] | None = None,
        strain_shear: str = "tensor",
    ) -> "TimeHistory":
        """Read a CSV with a time column and the six stress components (and optionally strain).

        Default columns: ``t``, ``sxx syy szz sxy syz sxz`` and ``exx eyy ezz exy eyz exz``; the
        names of the FE node CSV are accepted too (``time``, ``S_xx`` ... ``S_xz``, ``E_xx`` ...
        ``E_xz`` for the total strain, or ``Eel_*`` together with ``Epl_*``, whose sum is the
        total). ``columns`` renames them (keys ``time``, ``xx`` ... for stress, ``exx`` ... for
        strain) and then replaces the defaults for that quantity. Having both the default and the
        alias of one quantity is an error, and so is ``Eel_*`` alone (not the total strain under
        plasticity). Strain is read only if all six of its columns are present. Both units are
        required.
        """
        import pandas as pd

        path = Path(path)
        caller = dict(columns or {})
        frame = pd.read_csv(path, sep=delimiter)

        def pick(key: str) -> str | None:
            names = (
                (caller[key],) if key in caller else (DEFAULT_COLUMNS[key],) + _ALIAS_COLUMNS[key]
            )
            present = [c for c in names if c in frame.columns]
            if len(present) > 1:
                raise ValueError(
                    f"{path.name}: both {present[0]!r} and {present[1]!r} are present for one "
                    "quantity: keep one"
                )
            return present[0] if present else None

        time_col = pick("time")
        stress_cols = [pick(c) for c in CANONICAL]
        wanted = ["time", *CANONICAL]
        missing = [
            caller.get(key, DEFAULT_COLUMNS[key])
            for key, col in zip(wanted, [time_col, *stress_cols])
            if col is None
        ]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}; found {list(frame.columns)}")
        strain_cols = [pick("e" + c) for c in CANONICAL]
        has_strain = [c is not None for c in strain_cols]
        if any(has_strain) and not all(has_strain):
            gone = [DEFAULT_COLUMNS["e" + c] for c, ok in zip(CANONICAL, has_strain) if not ok]
            raise ValueError(f"{path.name}: strain columns incomplete, missing {gone}")
        S = components_to_sym3(frame[stress_cols].to_numpy(float), CANONICAL, kind="stress")
        split = {
            group: [f"{group}_{c}" for c in CANONICAL if f"{group}_{c}" in frame.columns]
            for group in ("Eel", "Epl")
        }
        E = None
        if all(has_strain):
            if split["Eel"] or split["Epl"]:
                raise ValueError(
                    f"{path.name}: the total strain E_* and the split Eel_*/Epl_* are both present"
                )
            E = components_to_sym3(
                frame[strain_cols].to_numpy(float), CANONICAL, kind="strain", shear=strain_shear
            )
        elif split["Eel"] or split["Epl"]:
            if len(split["Eel"]) != 6 or len(split["Epl"]) != 6:
                raise ValueError(
                    f"{path.name}: Eel_* columns alone are not the total strain under plasticity; "
                    "give the total E_*, or all of Eel_* and Epl_* (their sum is the total)"
                )
            E = sum(
                components_to_sym3(
                    frame[[f"{g}_{c}" for c in CANONICAL]].to_numpy(float),
                    CANONICAL,
                    kind="strain",
                    shear=strain_shear,
                )
                for g in ("Eel", "Epl")
            )
        _check_time_unit(time_unit)
        return cls.from_arrays(
            frame[time_col].to_numpy(float),
            S,
            E,
            stress_unit=stress_unit,
            axis_kind=axis_kind,
            time_unit=time_unit,
            provenance={"source": "csv", "file": path.name, "sha256": sha256_of(path)},
        )

    @classmethod
    def from_file(
        cls, path: str | Path, location_id: int | None = None, *, axis_kind: str = "time"
    ) -> "TimeHistory":
        """Read one location of a ``multifat.fe`` history file (needs ``h5py``, extra ``fe``)."""
        from multifat.fe.store import HistoryFile

        path = Path(path)
        with HistoryFile(path) as f:
            if f.meta.step_role == "unit_case":
                raise ValueError(
                    f"{path.name} holds unit_case fields, not a time history; "
                    "use node_histories with load_functions"
                )
            if location_id is None:
                if f.n_locations != 1:
                    raise ValueError(
                        f"{path.name} holds {f.n_locations} locations; a time history is one "
                        "location: pass location_id"
                    )
                position = np.array([0])
            else:
                position = np.flatnonzero(f.location_id == location_id)
                if len(position) == 0:
                    raise KeyError(f"location {location_id} is not in {path.name}")
                if len(position) != 1:
                    raise ValueError(
                        f"location {location_id} is ambiguous across instance/material labels; "
                        "use node_histories to select the complete location key"
                    )
            time = f.load_cases.time
            if time is None:
                raise ValueError(f"{path.name} has no time axis (load_cases/time)")
            block = f.read(position)
            provenance = {
                "source": "history_file",
                "file": path.name,
                "sha256": sha256_of(path),
                "location_id": int(block.location_id[0]),
                "instance": "" if f.instance is None else str(f.instance[position[0]]),
                "material": "" if f.material is None else str(f.material[position[0]]),
                "material_id": f.meta.material_id,
                "solver": f.meta.solver,
            }
        return cls(
            time, block.S_hist[0], None if block.E_hist is None else block.E_hist[0],
            axis_kind, {"stress": "MPa", "source_stress_unit": "MPa"}, provenance,
        )  # fmt: skip

    @classmethod
    def from_function(
        cls, fn: Callable[..., Any], time, *, name: str | None = None, **params: Any
    ) -> "TimeHistory":
        """Generate a history from a documented Python function ``fn(time, **params)``.

        ``fn`` returns the stress tensors in MPa, shape ``(N, 3, 3)``, or a dict with ``"S"`` and
        optionally ``"E"``. The function's qualified name (or ``name``) and ``params``, which must
        be JSON-serialisable, are recorded as provenance, with the sampled arrays.
        """
        _json_safe("params", params)
        time = np.asarray(time, dtype=float)
        out = fn(time, **params)
        if isinstance(out, dict):
            if "S" not in out:
                raise ValueError(
                    "a function returning a dict must give the stress under 'S' (shape (N, 3, 3))"
                )
            S, E = out["S"], out.get("E")
        else:
            S, E = out, None
        S = np.asarray(S, dtype=float)
        if S.shape != (len(time), 3, 3):
            raise ValueError(
                f"the function must return stress of shape ({len(time)}, 3, 3), got {S.shape}"
            )
        label = name or f"{getattr(fn, '__module__', '?')}.{getattr(fn, '__qualname__', repr(fn))}"
        return cls.from_arrays(
            time, S, E, stress_unit="MPa",
            provenance={"source": "function", "function": label, "params": params},
        )  # fmt: skip
