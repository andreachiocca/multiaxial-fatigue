"""Model complexity declarations and characterisation classes (issue #57, decision d116).

Defines the two-axis model complexity:
  - p_fit: number of quantities fitted on fatigue data (master curve, calibrated plane
    constants, native-law parameters).
  - p_in: number of distinct experimental characterisation classes required beyond the
    uniaxial/torsion S-N calibration paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Experimental characterisation classes required beyond the S-N calibration paths.
CHARACTERIZATION_CLASSES: dict[str, str] = {
    "Sigm1": "axial fatigue limit",
    "Taum1": "torsional fatigue limit",
    "Sig0": "mean-stress (R=0) data",
    "Tau0": "mean-stress (R=0) data",
    "Su": "static strength",
    "Sy": "static strength",
    "Tu": "static strength",
    "sigma_f": "S-N curves",
    "b_fat": "S-N curves",
    "tau_f": "S-N curves",
    "b0": "S-N curves",
    "epsilon_f": "strain-life / cyclic curve",
    "c_fat": "strain-life / cyclic curve",
    "gamma_f": "strain-life / cyclic curve",
    "c0": "strain-life / cyclic curve",
    "K_prime": "strain-life / cyclic curve",
    "n_prime": "strain-life / cyclic curve",
    "E": "elastic constants",
    "nu": "elastic constants",
}

# Known calibrated parameters on fatigue data (scatter minimisation on calibration paths).
CALIBRATED_CONSTANTS: dict[str, tuple[str, ...]] = {
    "FS": ("k_FS",),
    "FIN": ("k_FI",),
    "LI": ("L_LI",),
    "MGSE_ZHU": ("k_MGSE",),
    "MSWT": ("a_MSWT",),
    "PDS_FS": ("k_PDS",),
}

# Native models fitting quantities on fatigue data.
NATIVE_FITTED: dict[str, tuple[str, ...]] = {
    "MOREL": ("p", "q", "r"),
    "TVEIT2025": ("A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT"),
    "MDC": ("tau_f", "b0"),
    "MMCCM": ("sigma_f", "b_fat", "epsilon_f", "c_fat"),
    "LIUMAHADEVAN": ("sigma_f", "b_fat", "tau_f", "b0"),
}


@dataclass(frozen=True)
class ModelComplexity:
    """Complexity metrics for a fatigue criterion."""

    p_in: int
    p_fit: int
    input_classes: tuple[str, ...]
    fitted_constants: tuple[str, ...]


def model_fitted_constants(model: Any) -> tuple[str, ...]:
    """Return the names of constants fitted on fatigue data for this model."""
    if hasattr(model, "fitted_constants"):
        return tuple(getattr(model, "fitted_constants"))
    name = getattr(model, "name", str(model))
    if name in NATIVE_FITTED:
        return NATIVE_FITTED[name]
    if name in CALIBRATED_CONSTANTS:
        return CALIBRATED_CONSTANTS[name]
    return ()


def model_input_classes(model: Any, *, stress_only: bool | None = None) -> tuple[str, ...]:
    """Return experimental characterisation classes beyond calibration S-N paths."""
    if hasattr(model, "input_classes"):
        return tuple(sorted(set(getattr(model, "input_classes"))))

    name = getattr(model, "name", str(model))
    req = set(model.required_params()) if hasattr(model, "required_params") else set()
    fitted = set(model_fitted_constants(model))
    param_inputs = req - fitted

    classes = {CHARACTERIZATION_CLASSES[p] for p in param_inputs if p in CHARACTERIZATION_CLASSES}

    if stress_only is None:
        from .registry import stress_only_models

        stress_only = name in stress_only_models()

    if not stress_only:
        classes.add("elastic constants")

    return tuple(sorted(classes))


def get_model_complexity(model: Any, *, stress_only: bool | None = None) -> ModelComplexity:
    """Return the two-axis complexity (p_in, p_fit) for a model instance or name."""
    from .cp_adapter import CPModelAdapter
    from .registry import available_cp_methods, available_direct_models, stress_only_models

    name = getattr(model, "name", str(model))
    if isinstance(model, str):
        cp_dict = available_cp_methods()
        direct_dict = available_direct_models()
        if name in direct_dict:
            model = direct_dict[name]
        elif name in cp_dict:
            model = CPModelAdapter(cp_dict[name])
        else:
            raise KeyError(f"Unknown model name: {name}")

    if stress_only is None:
        stress_only = name in stress_only_models()

    fitted = model_fitted_constants(model)
    in_classes = model_input_classes(model, stress_only=stress_only)

    caps = getattr(model, "comparison_capabilities", None)
    is_native = caps is not None and getattr(caps, "life_prediction", "") == "native"

    if is_native:
        p_fit = len(fitted)
    else:
        # Standard master-curve models fit A and b on calibration data
        p_fit = 2 + len(fitted)

    p_in = len(in_classes)
    return ModelComplexity(
        p_in=p_in,
        p_fit=p_fit,
        input_classes=in_classes,
        fitted_constants=fitted,
    )
