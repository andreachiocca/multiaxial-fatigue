"""Declared comparison semantics for fatigue models.

Life prediction can be based on a fitted damage-parameter master curve or on
a model's native life law. Damage-parameter comparisons require a physically
meaningful DP-to-experimental-life reference; an auxiliary fatigue index alone
is not sufficient.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

LifePrediction = Literal["master_curve", "native"]


@dataclass(frozen=True)
class ComparisonCapabilities:
    """Method-level scientific support, independent of dataset eligibility."""

    life_prediction: LifePrediction = "master_curve"
    dp_comparison: bool = True

    def __post_init__(self):
        if self.life_prediction not in ("master_curve", "native"):
            raise ValueError(f"Unknown life prediction mechanism: {self.life_prediction!r}")


DEFAULT_CAPABILITIES = ComparisonCapabilities()


def capabilities_for_model(model) -> ComparisonCapabilities:
    """Obtain the registered model's comparison capabilities."""
    value = getattr(model, "comparison_capabilities", DEFAULT_CAPABILITIES)
    if not isinstance(value, ComparisonCapabilities):
        raise TypeError(f"{model.name} must declare ComparisonCapabilities")
    return value


def capabilities_for_name(name: str) -> ComparisonCapabilities:
    """Resolve a result method, including _ext and Findley amplitude labels."""
    from multifat.models.registry import get_models

    base = str(name).split(".", 1)[0].removesuffix("_ext")
    return capabilities_for_model(get_models([base])[0])
