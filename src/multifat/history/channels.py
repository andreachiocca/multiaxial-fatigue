"""Scalar channels of a tensor history, the signal a scalar counting algorithm runs on.

The default channel is the projection of the mean-removed stress path on its dominant direction. In
a proportional history, ``S(t) = lambda(t) A + S_static``, it is exactly affine in the load factor
``lambda(t)``, so counting on it counts the true load reversals even when a static stress is present.
The signed von Mises channel is not affine in ``lambda`` when ``S_static`` is nonzero; it is offered
explicitly, with its sign rule stated.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

from multifat.fe.tensors import sym3_to_six

_SHEAR_WEIGHT = np.sqrt(2.0)
_WEIGHTS = np.array([1.0, 1.0, 1.0, _SHEAR_WEIGHT, _SHEAR_WEIGHT, _SHEAR_WEIGHT])


def weighted_path(S_hist: np.ndarray) -> np.ndarray:
    """Mean-removed six-component stress path, shear weighted by sqrt(2) (Frobenius norm), ``(N, 6)``."""
    path = sym3_to_six(np.asarray(S_hist, dtype=float)) * _WEIGHTS
    return path - path.mean(axis=0)


def dominant(S_hist: np.ndarray) -> np.ndarray:
    """Projection of the stress path on its dominant singular direction, in MPa, shape ``(N,)``.

    The sign of the direction is fixed so that its largest component is positive. A constant
    history gives zeros.
    """
    path = weighted_path(S_hist)
    _, singular, rows = np.linalg.svd(path, full_matrices=False)
    direction = rows[0]
    if singular[0] == 0.0:
        return np.zeros(path.shape[0])
    if direction[int(np.argmax(np.abs(direction)))] < 0.0:
        direction = -direction
    return path @ direction


def component(S_hist: np.ndarray, i: int, j: int) -> np.ndarray:
    """One tensor component, shape ``(N,)``."""
    return np.asarray(S_hist, dtype=float)[:, i, j].copy()


def normal_stress(S_hist: np.ndarray, n) -> np.ndarray:
    """Normal stress on the plane with normal ``n`` (normalised here), shape ``(N,)``."""
    n = np.asarray(n, dtype=float)
    norm = np.linalg.norm(n)
    if norm == 0.0:
        raise ValueError("the plane normal must not be the zero vector")
    n = n / norm
    return np.einsum("i,tij,j->t", n, np.asarray(S_hist, dtype=float), n)


def signed_von_mises(S_hist: np.ndarray) -> np.ndarray:
    """Von Mises stress signed by the hydrostatic part, shape ``(N,)``.

    The sign is that of the trace; where the trace is zero (pure shear) it is that of the largest
    principal value, and positive if that is zero too. Not affine in the load factor when a static
    stress is present, so it is an explicit option, not the default.
    """
    S = np.asarray(S_hist, dtype=float)
    trace = np.trace(S, axis1=-2, axis2=-1)
    dev = S - trace[:, None, None] / 3.0 * np.eye(3)
    magnitude = np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1)))
    tol = 1e-9 * max(float(np.max(np.abs(S))), 1e-300)
    principal = np.linalg.eigvalsh(S)[:, -1]
    sign = np.where(np.abs(trace) > tol, np.sign(trace), np.where(principal < 0.0, -1.0, 1.0))
    return sign * magnitude


def resolve_channel(channel: Any, history) -> np.ndarray:
    """The scalar channel of a ``TimeHistory`` from a name, a tuple or a callable.

    Accepts ``"dominant"``, ``"von_mises"``, ``("component", i, j)``, ``("normal", n)`` or a callable
    ``f(S_hist) -> (N,)``; anything else raises ``ValueError``.
    """
    S = history.S_hist
    if isinstance(channel, str) and channel == "dominant":
        out = dominant(S)
    elif isinstance(channel, str) and channel == "von_mises":
        out = signed_von_mises(S)
    elif isinstance(channel, tuple) and channel and channel[0] == "component" and len(channel) == 3:
        out = component(S, channel[1], channel[2])
    elif isinstance(channel, tuple) and channel and channel[0] == "normal" and len(channel) == 2:
        out = normal_stress(S, channel[1])
    elif callable(channel):
        out = np.asarray(channel(S), dtype=float)
    else:
        raise ValueError(
            f"unknown channel {channel!r}: use 'dominant', 'von_mises', ('component', i, j), "
            "('normal', n) or a callable"
        )
    if out.shape != (history.n_samples,):
        raise ValueError(f"the channel must have shape ({history.n_samples},), got {out.shape}")
    if not np.all(np.isfinite(out)):
        raise ValueError("the channel must be finite")
    return out


ChannelSpec = str | tuple | Callable[[np.ndarray], np.ndarray]
