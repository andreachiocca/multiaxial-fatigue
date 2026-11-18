"""Damage accumulation into a life: a Basquin curve of the damage parameter and linear Miner.

The curve is ``DP = A N^b`` with ``b < 0``, the master curve the pipeline fits per criterion and
variant (``A_surv`` and ``b`` of its in-memory fit). It is supplied by the caller: the repository does
not persist it, and nothing here changes a calibration record. The damage of an event is
``count / N_f(DP)``; events with a non-positive damage parameter contribute nothing. Below an optional
damage-parameter limit events can be ignored; otherwise the curve is extrapolated at every level.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np

from multifat.eval.fat_cycles import fatigue_cycles


@dataclass(frozen=True)
class SNCurve:
    """``DP = A N^b`` of one criterion and variant, with the provenance of the fit."""

    A: float
    b: float
    method: str
    variant: str = ""
    design_loading: str = "tension"
    source: str = ""

    def __post_init__(self):
        if not (np.isfinite(self.A) and self.A > 0.0):
            raise ValueError(f"A must be positive and finite, got {self.A}")
        if not np.isfinite(self.b):
            raise ValueError(f"b must be finite, got {self.b}")
        if self.b >= 0.0:
            raise ValueError(f"b = {self.b} >= 0: the curve gives no usable life")

    @classmethod
    def from_fit(
        cls,
        fit: Mapping[str, float],
        *,
        method: str,
        variant: str = "",
        design_loading: str = "tension",
        source: str = "",
    ) -> "SNCurve":
        """From the pipeline's in-memory fit (``cal["fit"]`` or ``cal["fit_ext"]``): ``A_surv``, ``b``."""
        return cls(float(fit["A_surv"]), float(fit["b"]), method, variant, design_loading, source)

    def cycles_to_failure(self, dp):
        """Cycles to failure at damage parameter ``dp`` (infinite where ``dp <= 0``)."""
        dp = np.asarray(dp, dtype=float)
        out = np.full(dp.shape, np.inf)
        positive = dp > 0.0
        if np.any(positive):
            with np.errstate(
                over="ignore", under="ignore"
            ):  # a tiny dp is an effectively infinite life
                out[positive] = fatigue_cycles(dp[positive], self.A, self.b)
        return out

    def dp_at(self, n_cycles: float) -> float:
        """Damage parameter at ``n_cycles`` cycles to failure."""
        return float(self.A * n_cycles**self.b)


@dataclass(frozen=True)
class Endurance:
    """Treatment of events at low damage parameter: ``extrapolate`` or ``ignore_below`` ``dp_limit``."""

    mode: Literal["extrapolate", "ignore_below"] = "extrapolate"
    dp_limit: float | None = None

    def __post_init__(self):
        if self.mode not in ("extrapolate", "ignore_below"):
            raise ValueError(f"mode must be 'extrapolate' or 'ignore_below', got {self.mode!r}")
        if self.mode == "ignore_below" and not (
            self.dp_limit is not None and np.isfinite(self.dp_limit) and self.dp_limit > 0.0
        ):
            raise ValueError("ignore_below needs a positive finite dp_limit")

    @classmethod
    def ignore_below(cls, dp_limit: float) -> "Endurance":
        """Events with a damage parameter below ``dp_limit`` contribute no damage."""
        return cls("ignore_below", float(dp_limit))

    @classmethod
    def from_knee(cls, curve: SNCurve, n_knee: float) -> "Endurance":
        """The limit at the knee life: ``dp_limit = A n_knee^b``."""
        return cls.ignore_below(curve.dp_at(n_knee))


@dataclass(frozen=True)
class MinerResult:
    """Damage and life per event and the linear sum."""

    damage: np.ndarray
    n_f: np.ndarray
    total: float


def miner(dp, counts, curve: SNCurve, endurance: Endurance | None = None) -> MinerResult:
    """Linear (Palmgren-Miner) sum of ``count / N_f(dp)`` over events.

    Parameters
    ----------
    dp, counts : array_like
        Damage parameter and count of each event (equal shape; counts not negative).
    curve : SNCurve
        Basquin curve of the same criterion and variant.
    endurance : Endurance, optional
        Default: ``extrapolate``.
    """
    endurance = endurance or Endurance()
    dp = np.asarray(dp, dtype=float)
    counts = np.asarray(counts, dtype=float)
    if dp.shape != counts.shape:
        raise ValueError(
            f"dp and counts must have the same shape, got {dp.shape} and {counts.shape}"
        )
    if np.any(counts < 0.0) or not np.all(np.isfinite(counts)):
        raise ValueError("counts must be finite and not negative")
    n_f = curve.cycles_to_failure(dp)
    damage = np.zeros(dp.shape)
    live = np.isfinite(n_f)
    if endurance.mode == "ignore_below":
        live &= dp >= endurance.dp_limit
    damage[live] = counts[live] / n_f[live]
    return MinerResult(damage=damage, n_f=n_f, total=float(np.sum(damage)))
