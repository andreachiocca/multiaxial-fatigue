"""Tensor conventions of the FE interface: component order, shear convention, units and frames.

Every adapter converts to the repository convention before data enters a history: symmetric 3x3
tensors in MPa (stress) or dimensionless tensor strain (shear strain = engineering shear / 2),
lengths in mm, global Cartesian axes. Nothing is inferred: orders, shear convention and units are
arguments.
"""

from __future__ import annotations

from typing import Literal, Sequence

import numpy as np

CANONICAL = ("xx", "yy", "zz", "xy", "yz", "xz")
VTK9 = ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")

_INDEX = {"x": 0, "y": 1, "z": 2}
_STRESS_TO_MPA = {
    "Pa": 1e-6,
    "kPa": 1e-3,
    "MPa": 1.0,
    "GPa": 1e3,
    "N/mm^2": 1.0,
    "psi": 6.894757293168361e-3,
    "ksi": 6.894757293168361,
}
_LENGTH_TO_MM = {"m": 1e3, "cm": 10.0, "mm": 1.0, "in": 25.4, "ft": 304.8}
_STRAIN_TO_UNITLESS = {"1": 1.0, "mm/mm": 1.0, "m/m": 1.0, "percent": 1e-2, "microstrain": 1e-6}
STRESS_UNITS = tuple(_STRESS_TO_MPA)
LENGTH_UNITS = tuple(_LENGTH_TO_MM)
STRAIN_UNITS = tuple(_STRAIN_TO_UNITLESS)
_SYM_RTOL = 1e-9


def sym3_to_six(T: np.ndarray) -> np.ndarray:
    """Pack symmetric tensors ``(..., 3, 3)`` into the canonical six components ``(..., 6)``."""
    T = np.asarray(T)
    return np.stack(
        [T[..., 0, 0], T[..., 1, 1], T[..., 2, 2], T[..., 0, 1], T[..., 1, 2], T[..., 0, 2]],
        axis=-1,
    )


def six_to_sym3(v: np.ndarray) -> np.ndarray:
    """Expand canonical six components ``(..., 6)`` into symmetric tensors ``(..., 3, 3)``."""
    v = np.asarray(v)
    T = np.zeros(v.shape[:-1] + (3, 3), dtype=v.dtype)
    for k, (i, j) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))):
        T[..., i, j] = v[..., k]
        T[..., j, i] = v[..., k]
    return T


def _label_to_ij(label: str) -> tuple[int, int]:
    if len(label) != 2 or label[0] not in _INDEX or label[1] not in _INDEX:
        raise ValueError(f"unknown tensor component label {label!r}; use xx, xy, ... zz")
    return _INDEX[label[0]], _INDEX[label[1]]


def components_to_sym3(
    a: np.ndarray,
    order: Sequence[str],
    *,
    kind: Literal["stress", "strain"],
    shear: Literal["tensor", "engineering"] = "tensor",
) -> np.ndarray:
    """Build symmetric tensors ``(M, 3, 3)`` from component columns ``(M, k)``.

    Parameters
    ----------
    a : ndarray
        One row per point, one column per entry of ``order``.
    order : sequence of str
        Component label of each column, e.g. ``CANONICAL`` or ``VTK9``. When both ``ij`` and
        ``ji`` are present they must agree.
    kind : {"stress", "strain"}
        Only strain is affected by ``shear``.
    shear : {"tensor", "engineering"}
        Shear convention of the source strain; engineering shear is halved.

    Returns
    -------
    ndarray
        Symmetric tensors in the tensor-shear convention.
    """
    if kind not in ("stress", "strain"):
        raise ValueError(f"kind must be 'stress' or 'strain', got {kind!r}")
    if shear not in ("tensor", "engineering"):
        raise ValueError(f"shear must be 'tensor' or 'engineering', got {shear!r}")
    a = np.asarray(a, dtype=float)
    if a.ndim != 2 or a.shape[1] != len(order):
        raise ValueError(f"expected {len(order)} columns for order {tuple(order)}, got {a.shape}")
    if not np.all(np.isfinite(a)):
        raise ValueError("tensor components must be finite")
    if len(set(order)) != len(order):
        raise ValueError(f"duplicate component labels in order {tuple(order)}")
    T = np.full((a.shape[0], 3, 3), np.nan)
    for col, label in enumerate(order):
        i, j = _label_to_ij(label)
        T[:, i, j] = a[:, col]
    for col, label in enumerate(order):
        i, j = _label_to_ij(label)
        if i == j:
            continue
        if np.isnan(T[:, j, i]).all():
            T[:, j, i] = a[:, col]
        else:
            scale = max(float(np.max(np.abs(a))), 1e-300)
            if np.max(np.abs(T[:, j, i] - a[:, col])) > _SYM_RTOL * scale:
                raise ValueError(
                    f"tensor is not symmetric: components {label!r} and its pair differ"
                )
    if np.isnan(T).any():
        raise ValueError(f"order {tuple(order)} does not define all nine tensor entries")
    if kind == "strain" and shear == "engineering":
        off = ~np.eye(3, dtype=bool)
        T[:, off] *= 0.5
    return T


def to_mpa(a: np.ndarray, unit: str) -> np.ndarray:
    """Scale stress values given in ``unit`` (Pa, kPa, MPa, GPa, N/mm^2, psi, ksi) to MPa."""
    if unit not in _STRESS_TO_MPA:
        raise ValueError(f"unknown stress unit {unit!r}; known: {sorted(_STRESS_TO_MPA)}")
    return np.asarray(a, dtype=float) * _STRESS_TO_MPA[unit]


def to_mm(a: np.ndarray, unit: str) -> np.ndarray:
    """Scale lengths given in ``unit`` (m, cm, mm, in, ft) to mm."""
    if unit not in _LENGTH_TO_MM:
        raise ValueError(f"unknown length unit {unit!r}; known: {sorted(_LENGTH_TO_MM)}")
    return np.asarray(a, dtype=float) * _LENGTH_TO_MM[unit]


def to_strain(a: np.ndarray, unit: str) -> np.ndarray:
    """Scale strain values given in ``unit`` (1, mm/mm, m/m, percent, microstrain) to unity."""
    if unit not in _STRAIN_TO_UNITLESS:
        raise ValueError(f"unknown strain unit {unit!r}; known: {sorted(_STRAIN_TO_UNITLESS)}")
    return np.asarray(a, dtype=float) * _STRAIN_TO_UNITLESS[unit]


def elastic_strain(S: np.ndarray, E: float, nu: float) -> np.ndarray:
    """Isotropic Hooke strain ``(1 + nu) S / E - nu tr(S) I / E`` of stress tensors in MPa.

    Parameters
    ----------
    S : ndarray
        Stress tensors ``(..., 3, 3)`` in MPa.
    E : float
        Young modulus in MPa.
    nu : float
        Poisson ratio, in the open interval (-1, 0.5).

    Returns
    -------
    ndarray
        Elastic strain tensors (tensor shear convention).
    """
    if not E > 0.0:
        raise ValueError(f"Young modulus must be positive, got {E}")
    if not -1.0 < nu < 0.5:
        raise ValueError(f"Poisson ratio must lie in (-1, 0.5), got {nu}")
    S = np.asarray(S, dtype=float)
    trace = np.trace(S, axis1=-2, axis2=-1)[..., None, None]
    return ((1.0 + nu) * S - nu * trace * np.eye(3)) / E


def rotate_to_global(T: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Rotate tensors from a local frame to the global one: ``R T R^T``.

    The columns of ``R`` are the local axes written in global coordinates.
    """
    T = np.asarray(T, dtype=float)
    R = np.asarray(R, dtype=float)
    return R @ T @ np.swapaxes(R, -1, -2)


def cylindrical_rotations(
    points: np.ndarray,
    axis: Sequence[float] = (0.0, 0.0, 1.0),
    origin: Sequence[float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Local ``(r, theta, z)`` frames at each point, as rotations ``(M, 3, 3)``.

    The columns are ``e_r``, ``e_theta = e_z x e_r`` and ``e_z`` (the unit ``axis``). A point on
    the axis has no radial direction and is rejected.
    """
    points = np.asarray(points, dtype=float)
    ez = np.asarray(axis, dtype=float)
    ez = ez / np.linalg.norm(ez)
    rel = points - np.asarray(origin, dtype=float)
    radial = rel - (rel @ ez)[:, None] * ez
    norm = np.linalg.norm(radial, axis=1)
    if np.any(norm <= 1e-12 * (1.0 + np.linalg.norm(rel, axis=1))):
        raise ValueError("a point lies on the cylinder axis: its radial direction is undefined")
    er = radial / norm[:, None]
    et = np.cross(ez, er)
    return np.stack([er, et, np.broadcast_to(ez, er.shape)], axis=-1)
