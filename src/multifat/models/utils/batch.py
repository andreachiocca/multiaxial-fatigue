"""Helpers for batched (all-planes-at-once) critical-plane evaluation.

Batched formulas must reproduce the scalar ``evaluate_on_plane`` results bit for bit,
so they repeat its operation order. One trap: CPython's ``x**2`` on a float calls the C
library ``pow``, which on some platforms (MSVC) is not correctly rounded and can differ
from ``x * x`` by one ulp; ``pypow2`` reproduces it elementwise.
"""

import numpy as np


def pypow2(x: np.ndarray) -> np.ndarray:
    """Elementwise ``float(x)**2`` exactly as CPython computes it.

    Evaluates ``v**2`` in Python rather than via ``np.frompyfunc`` to avoid
    object-dtype ufunc state while matching CPython's C ``pow`` bit for bit.
    """
    arr = np.asarray(x, dtype=float)
    return np.array([v**2 for v in arr.ravel().tolist()], dtype=float).reshape(arr.shape)


def check_planes(*arrays: np.ndarray) -> None:
    """Shape check matching the scalar path (which does not reject non-finite values)."""
    for a in arrays:
        if a.ndim != 3 or a.shape[1:] != (3, 3):
            raise ValueError(f"Batched plane tensors must have shape (P, 3, 3); got {a.shape}")


def tau_max_planes(S0r, S1r, params, tau_max=None) -> np.ndarray:
    """Per-plane maximum shear stress: path value if supplied, else the larger endpoint."""
    if tau_max is not None:
        return np.asarray(tau_max, dtype=float)
    if "_tau_max" in params:
        return np.full(S0r.shape[0], float(params["_tau_max"]))
    tau0 = np.hypot(S0r[:, 0, 2], S0r[:, 1, 2])
    tau1 = np.hypot(S1r[:, 0, 2], S1r[:, 1, 2])
    return np.maximum(tau0, tau1)


def elastic_strain_planes(S: np.ndarray, *, E: float, nu: float) -> np.ndarray:
    """Batched isotropic elastic strain, same operation order as the scalar helpers."""
    tr = np.trace(S, axis1=1, axis2=2)
    return ((1.0 + nu) / E) * S - ((nu / E) * tr)[:, None, None] * np.eye(3)
