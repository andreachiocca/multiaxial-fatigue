"""Selection of FE locations by ID, mask, box and screening scalar, without loading everything.

Every function returns sorted, unique file positions, ready for ``HistoryFile.iter_chunks``.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from .model import HistoryBlock
from .store import HistoryFile


def select_by_ids(
    f: HistoryFile,
    ids: Sequence[int],
    *,
    instance: str | None = None,
    material: str | None = None,
) -> np.ndarray:
    """Positions of the given location IDs; an unknown ID raises ``KeyError``.

    A node label may occur in several instances or materials; when ``instance`` or ``material``
    is omitted every matching location is returned, so an interface node gives one position per
    material. An ID with no match under the given filters raises ``KeyError``.
    """
    ids = np.asarray(ids)
    keep = np.ones(f.n_locations, dtype=bool)
    if instance is not None:
        keep &= f.instance is not None and f.instance == instance
    if material is not None:
        keep &= f.material is not None and f.material == material
    candidates = np.flatnonzero(keep)
    hit = np.isin(f.location_id[candidates], ids)
    found = candidates[hit]
    missing = np.setdiff1d(ids, f.location_id[found])
    if len(missing):
        raise KeyError(f"unknown location IDs: {missing[:10].tolist()}")
    return found


def select_by_mask(f: HistoryFile, mask: np.ndarray) -> np.ndarray:
    """Positions where a boolean ``mask`` (one entry per location) is true."""
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != (f.n_locations,):
        raise ValueError(f"mask must have shape ({f.n_locations},), got {mask.shape}")
    return np.flatnonzero(mask)


def select_in_box(f: HistoryFile, lo: Sequence[float], hi: Sequence[float]) -> np.ndarray:
    """Positions with coordinates (mm) inside the closed box ``lo <= x <= hi``."""
    xyz = f.coordinates
    return np.flatnonzero(np.all((xyz >= np.asarray(lo)) & (xyz <= np.asarray(hi)), axis=1))


def von_mises_max(block: HistoryBlock) -> np.ndarray:
    """Largest von Mises stress of each location over the load cases (a screening aid)."""
    S = block.S_hist
    dev = S - np.trace(S, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    return np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1))).max(axis=1)


def screen_by(
    f: HistoryFile,
    scalar: Callable[[HistoryBlock], np.ndarray],
    *,
    top: int | None = None,
    threshold: float | None = None,
    max_locations: int = 4096,
) -> np.ndarray:
    """Shortlist locations by a scalar computed chunk by chunk.

    Exactly one of ``top`` (the ``top`` largest values) and ``threshold`` (values at or above it)
    is given. NaN values are ignored. The scalar is a screening aid, not a damage claim.
    """
    if (top is None) == (threshold is None):
        raise ValueError("give exactly one of top and threshold")
    if top is not None and top < 1:
        raise ValueError(f"top must be a positive number of locations, got {top}")
    keep_val = np.empty(0)
    keep_pos = np.empty(0, dtype=int)
    start = 0
    for block in f.iter_chunks(max_locations=max_locations):
        values = np.asarray(scalar(block), dtype=float)
        positions = np.arange(start, start + block.n_locations)
        start += block.n_locations
        valid = ~np.isnan(values)  # a location without a value is never shortlisted
        values, positions = values[valid], positions[valid]
        if threshold is not None:
            hit = values >= threshold
            keep_val, keep_pos = np.append(keep_val, values[hit]), np.append(
                keep_pos, positions[hit]
            )
            continue
        keep_val = np.concatenate([keep_val, values])
        keep_pos = np.concatenate([keep_pos, positions])
        if len(keep_val) > top:
            best = np.argpartition(keep_val, len(keep_val) - top)[len(keep_val) - top :]
            keep_val, keep_pos = keep_val[best], keep_pos[best]
    return np.sort(keep_pos)
