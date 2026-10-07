from __future__ import annotations

"""Empirical estimation of missing fatigue properties.

The Excel material files used by this repository may not contain full
strain-life parameters (σ'_f, ε'_f, b, c, etc.).

When they are missing, we estimate them using empirical relations from
Zhao et al. (2020) *Fatigue Properties Estimation and Life Prediction for
Steels under Axial, Torsional, and In-Phase Loading*.

User requirement: use formulas **not dependent on hardness**.
Therefore, we implement the *Uniform Material Law* variant (Baumel & Seeger)
as summarized in Zhao et al. (2020):

  σ'_f = 1.5 σ_u
  b    = -0.087
  ε'_f = 0.59 ψ
  c    = -0.58
  ψ = 1                       if (σ_u/E) <= 0.003
      1.375 - 125 (σ_u/E)     otherwise

These relations were developed for steels. We still use them as a fallback
when the user does not provide fatigue properties, but we emit a warning.
"""

from dataclasses import dataclass
from typing import Mapping
import math
import warnings


@dataclass(frozen=True)
class AxialFatigueProps:
    sigma_f: float  # σ'_f  [MPa]
    epsilon_f: float  # ε'_f [-]
    b: float
    c: float


@dataclass(frozen=True)
class ShearFatigueProps:
    tau_f: float  # τ'_f [MPa]
    gamma_f: float  # γ'_f [-]
    b0: float
    c0: float


def estimate_axial_fatigue_props_uniform_material_law(*, Su: float, E: float) -> AxialFatigueProps:
    """Estimate (σ'_f, ε'_f, b, c) from Su and E (no hardness)."""
    Su = float(Su)
    E = float(E)
    if Su <= 0 or E <= 0:
        raise ValueError(f"Su and E must be > 0. Got Su={Su}, E={E}")

    sigma_f = 1.5 * Su
    b = -0.087

    r = Su / E
    psi = 1.0 if r <= 0.003 else (1.375 - 125.0 * r)
    epsilon_f = 0.59 * psi
    c = -0.58
    return AxialFatigueProps(sigma_f=float(sigma_f), epsilon_f=float(epsilon_f), b=float(b), c=float(c))


def ensure_axial_fatigue_props(params: Mapping[str, float]) -> AxialFatigueProps:
    """Return axial fatigue props from params, estimating if missing."""
    # Prefer user-provided values
    have = all(k in params and params[k] is not None for k in ("sigma_f", "epsilon_f", "b_fat", "c_fat"))
    if have:
        return AxialFatigueProps(
            sigma_f=float(params["sigma_f"]),
            epsilon_f=float(params["epsilon_f"]),
            b=float(params["b_fat"]),
            c=float(params["c_fat"]),
        )

    # Estimate from Su and E
    Su = params.get("Su", None)
    E = params.get("E", None)
    if Su is None or E is None:
        raise KeyError("Missing axial fatigue properties and cannot estimate because Su and/or E are missing.")

    warnings.warn(
        "Axial fatigue properties (sigma_f, epsilon_f, b_fat, c_fat) were not present in the material data file. "
        "Estimating them using the Uniform Material Law (Zhao et al., 2020). "
        "Note: this empirical law is validated for steels.",
        RuntimeWarning,
        stacklevel=2,
    )
    return estimate_axial_fatigue_props_uniform_material_law(Su=float(Su), E=float(E))


def ensure_shear_fatigue_props(params: Mapping[str, float]) -> ShearFatigueProps:
    """Return shear fatigue props, deriving from axial if missing.

    Uses the common von Mises relations (also used in several of the attached
    multiaxial fatigue papers):

      τ'_f = σ'_f / √3
      γ'_f = √3 ε'_f
      b0 = b
      c0 = c
    """
    have = all(k in params and params[k] is not None for k in ("tau_f", "gamma_f", "b0", "c0"))
    if have:
        return ShearFatigueProps(
            tau_f=float(params["tau_f"]),
            gamma_f=float(params["gamma_f"]),
            b0=float(params["b0"]),
            c0=float(params["c0"]),
        )

    warnings.warn(
        "Shear fatigue properties (tau_f, gamma_f, b0, c0) were not present in the material data file. "
        "Deriving them from axial fatigue properties using von-Mises relations: "
        "tau_f = sigma_f/sqrt(3), gamma_f = sqrt(3)*epsilon_f, b0=b, c0=c.",
        RuntimeWarning,
        stacklevel=2,
    )

    axial = ensure_axial_fatigue_props(params)
    rt3 = math.sqrt(3.0)
    return ShearFatigueProps(
        tau_f=float(axial.sigma_f / rt3),
        gamma_f=float(axial.epsilon_f * rt3),
        b0=float(axial.b),
        c0=float(axial.c),
    )
