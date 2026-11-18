"""Non-proportionality of a stress history, the gate of the exact (proportional) domain.

A proportional history, ``S(t) = lambda(t) A + S_static``, traces a straight line in stress space:
its mean-removed six-component path has rank one. The measure, ``rank_ratio``, is the second singular
value of that path over the first; the verdict compares it with a tolerance (default 1e-6, a
numerical floor for sampled functions and noisy data, not a physical threshold: a user who wants a
looser or stricter reading passes ``tol``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .channels import weighted_path


@dataclass(frozen=True)
class Proportionality:
    """Ratio of the second to the first singular value of the stress path, and the verdict."""

    rank_ratio: float
    proportional: bool
    tol: float


def proportionality(S_hist: np.ndarray, *, tol: float = 1e-6) -> Proportionality:
    """Measure how far a stress history is from a straight line in stress space.

    A static stress (including a hydrostatic one) does not change the verdict. A constant history
    has ``rank_ratio`` zero.
    """
    singular = np.linalg.svd(weighted_path(S_hist), compute_uv=False)
    rank_ratio = 0.0 if singular[0] == 0.0 else float(singular[1] / singular[0])
    return Proportionality(rank_ratio=rank_ratio, proportional=rank_ratio <= tol, tol=float(tol))
