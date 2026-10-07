"""Plane quantities from actual harmonic tensors, before scalar extrema encoding."""
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=1)
def _cycle():
    t = np.arange(720) * (2 * np.pi / 720)
    return np.sin(t), np.cos(t)


def shear_maximum(mean, sin, cos):
    """Maximum magnitude of shifted shear ellipse (720 time samples)."""
    st, ct = _cycle()
    path = mean[:2, 2, None] + sin[:2, 2, None]*st + cos[:2, 2, None]*ct
    return float(np.sqrt(np.max(np.sum(path*path, axis=0))))


def strain_ranges(sin, cos):
    """Engineering shear range and normal range; legacy RSS ellipse convention."""
    return (4 * float(np.sqrt(np.sum(sin[:2, 2]**2 + cos[:2, 2]**2))),
            2 * float(np.hypot(sin[2, 2], cos[2, 2])))


def split_strain_ranges(S_sin, S_cos, E_sin, E_cos, E, nu):
    def elastic(s):
        return ((1+nu)*s - nu*np.trace(s)*np.eye(3))/E
    es, ec = elastic(S_sin), elastic(S_cos)
    ge, ne = strain_ranges(es, ec)
    gp, np_ = strain_ranges(E_sin-es, E_cos-ec)
    return ge, gp, ne, np_
