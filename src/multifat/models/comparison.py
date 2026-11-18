"""Declared comparison semantics for fatigue models.

Life prediction can be based on a fitted damage-parameter master curve or on
a model's native life law. Damage-parameter comparisons require a physically
meaningful DP-to-experimental-life reference; an auxiliary fatigue index alone
is not sufficient.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

LifePrediction = Literal["master_curve", "native"]
DPReference = Literal["master_curve", "model_law", "none"]


@dataclass(frozen=True)
class ComparisonCapabilities:
    """Method-level scientific support, independent of dataset eligibility."""

    life_prediction: LifePrediction = "master_curve"
    dp_comparison: bool = True
    dp_reference: DPReference = "master_curve"
    load_error: bool = True

    def __post_init__(self):
        if self.life_prediction not in ("master_curve", "native"):
            raise ValueError(f"Unknown life prediction mechanism: {self.life_prediction!r}")
        if self.dp_reference not in ("master_curve", "model_law", "none"):
            raise ValueError(f"Unknown DP reference: {self.dp_reference!r}")
        if self.dp_comparison == (self.dp_reference == "none"):
            raise ValueError("DP comparison and reference declarations disagree")


DEFAULT_CAPABILITIES = ComparisonCapabilities()

# Method-level capability overrides for models owned externally or declared without newer fields.
_MODEL_CAPABILITY_OVERRIDES: dict[str, dict] = {
    # TVEIT2025 integrates damage cycle-by-cycle, making iterative load scaling prohibitively slow.
    "TVEIT2025": {"load_error": False},
}


def capabilities_for_model(model) -> ComparisonCapabilities:
    """Obtain the registered model's comparison capabilities."""
    value = getattr(model, "comparison_capabilities", DEFAULT_CAPABILITIES)
    if not isinstance(value, ComparisonCapabilities):
        raise TypeError(f"{model.name} must declare ComparisonCapabilities")
    name = getattr(model, "name", None)
    if name in _MODEL_CAPABILITY_OVERRIDES:
        return replace(value, **_MODEL_CAPABILITY_OVERRIDES[name])
    return value


def capabilities_for_name(name: str) -> ComparisonCapabilities:
    """Resolve a result method, including _ext and Findley amplitude labels."""
    from multifat.models.registry import get_models

    base = str(name).split(".", 1)[0].removesuffix("_ext")
    return capabilities_for_model(get_models([base])[0])
