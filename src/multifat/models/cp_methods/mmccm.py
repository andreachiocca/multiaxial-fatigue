# multifat/models/cp_methods/mmccm.py
"""Modified Manson-Coffin Curve Method proposed by Wang and Susmel (2016),
International Journal of Fatigue, doi:10.1016/j.ijfatigue.2015.10.005 (roadmap G22),
in constant-amplitude form, using Eqs. (1)-(8):

    rho = sigma_n,max / tau_a on the plane with the largest shear strain amplitude (Eq. 1);
    gamma_a = tau'_f(rho)/G (2N_f)^b(rho) + gamma'_f(rho) (2N_f)^c(rho)             (Eq. 2)
    tau'_f(rho)/G = rho (1 + nu_e) sigma'_f/E + (1 - rho) tau'_f/G                  (Eq. 5)
    gamma'_f(rho) = rho (1 + nu_p) eps'_f + (1 - rho) gamma'_f                       (Eq. 6)
    b(rho) = b b0 / ((b0 - b) rho + b),  c(rho) = c c0 / ((c0 - c) rho + c)          (Eqs. 7-8)

Here gamma_a denotes the engineering shear strain amplitude (the longest-chord rule is used for
harmonic paths; the maximum-variance direction used in the paper gives the same half range for
a single harmonic). When several planes share the maximum gamma_a, the paper chooses the one
with the largest rho (p.5); since a larger rho lowers the curve at fixed gamma_a, the repository's
tie rule based on larger damage picks that same plane. Eq. (2) is solved for N_f by bisection in
log(2N_f) over [1, 2e12]; the damage parameter is 1/N_f (zero when N_f exceeds 1e12 cycles).
The repository uses the workbook (elastic) strains and requires the strain-life constants sigma_f,
epsilon_f, b_fat, c_fat, tau_f, gamma_f, b0 and c0 (no estimation fallback is provided). nu_e = nu,
while nu_p = ``nu_p`` when supplied and 0.5 otherwise (the paper does not specify it). No rho
cut-off is given; a non-positive coefficient or a non-negative exponent (beyond a pole of
Eqs. 7-8) yields infinite damage, so the method has no ``_ext`` variant
(``ext_variant = False``).
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from multifat.models.errors import MissingParamsError
from multifat.models.utils.batch import check_planes, pypow2

from ..base import CPPlaneResult
from ..comparison import ComparisonCapabilities

_REQUIRED_STRAIN_KEYS = (
    "E",
    "nu",
    "sigma_f",
    "epsilon_f",
    "b_fat",
    "c_fat",
    "tau_f",
    "gamma_f",
    "b0",
    "c0",
)
_MAX_LOG_CYCLES = np.log(2.0e12)


def _validated_constants(params: Mapping[str, float]) -> dict:
    values = {}
    for key in _REQUIRED_STRAIN_KEYS:
        raw_value = params.get(key, None)
        if raw_value is None or not np.isfinite(float(raw_value)):
            raise MissingParamsError(
                "MMCCM", [key], "Required strain-life constants along with E and nu"
            )
        values[key] = float(raw_value)
    plastic_nu = params.get("nu_p", None)
    values["nu_p"] = 0.5 if plastic_nu is None else float(plastic_nu)
    return values


def mmccm_curve(rho, material: dict):
    """Return (A, b, B, c) for gamma_a = A (2N)^b + B (2N)^c at the given rho, Eqs. (5)-(8)."""
    rho = np.asarray(rho, dtype=float)
    G = material["E"] / (2.0 * (1.0 + material["nu"]))
    A = (
        rho * (1.0 + material["nu"]) * material["sigma_f"] / material["E"]
        + (1.0 - rho) * material["tau_f"] / G
    )
    B = rho * (1.0 + material["nu_p"]) * material["epsilon_f"] + (1.0 - rho) * material["gamma_f"]
    with np.errstate(divide="ignore", invalid="ignore"):
        b = (
            material["b_fat"]
            * material["b0"]
            / ((material["b0"] - material["b_fat"]) * rho + material["b_fat"])
        )
        c = (
            material["c_fat"]
            * material["c0"]
            / ((material["c0"] - material["c_fat"]) * rho + material["c_fat"])
        )
    return A, b, B, c


def mmccm_life(gamma_a, rho, material: dict) -> np.ndarray:
    """Return N_f solving Eq. (2).

    Curves that are not admissible (non-positive A or B, or a non-negative exponent) receive
    zero life, which means infinite damage. An amplitude below the curve at 2N = 2e12 receives
    infinite life.
    """
    gamma_a = np.asarray(gamma_a, dtype=float)
    A, b, B, c = (np.broadcast_to(value, gamma_a.shape) for value in mmccm_curve(rho, material))
    usable = (A > 0) & (B > 0) & (b < 0) & (c < 0) & np.isfinite(b) & np.isfinite(c)

    def model_curve(log_cycles):
        return A * np.exp(b * log_cycles) + B * np.exp(c * log_cycles)

    lower = np.zeros(gamma_a.shape)
    upper = np.full(gamma_a.shape, _MAX_LOG_CYCLES)
    with np.errstate(
        invalid="ignore", over="ignore"
    ):  # Out-of-domain curves get filtered afterwards
        for _ in range(80):
            midpoint = 0.5 * (lower + upper)
            is_above = model_curve(midpoint) > gamma_a
            lower, upper = np.where(is_above, midpoint, lower), np.where(is_above, upper, midpoint)
        Nf = 0.5 * np.exp(0.5 * (lower + upper))
        Nf = np.where(model_curve(0.0) <= gamma_a, 0.5, Nf)
        Nf = np.where(model_curve(_MAX_LOG_CYCLES) > gamma_a, np.inf, Nf)
    return np.where(usable, Nf, 0.0)


class MMCCM:
    """MMCCM critical plane: maximum engineering shear strain amplitude, damage 1/N_f (Eq. 2)."""

    name = "MMCCM"
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_reference="model_law"
    )
    optional_params = {"nu_p": "0.5 (plastic Poisson's ratio)"}
    np_shear_rule = "longest_chord"
    ext_variant = False

    def required_params(self) -> set[str]:
        return set(_REQUIRED_STRAIN_KEYS)

    def native_life(self, aux, params: Mapping[str, float], fit=None):
        """Return N_f from Eq. (2); the damage parameter is 1/N_f."""
        DP = np.asarray(aux["dp"], dtype=float)
        return np.where(
            DP > 0.0,
            1.0 / np.where(DP > 0.0, DP, 1.0),
            np.inf,
        )

    def expected_dp_at_life(self, N_exp, params: Mapping[str, float]):
        """Definition of the MMCCM damage parameter, DP = 1/N_f."""
        n = np.asarray(N_exp, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(n > 0, 1.0 / n, np.nan)

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]):
        damage, metric = self.evaluate_on_planes(
            S0r=np.asarray(S0r, float)[None],
            E0r=np.asarray(E0r, float)[None],
            S1r=np.asarray(S1r, float)[None],
            E1r=np.asarray(E1r, float)[None],
            params=params,
        )
        return CPPlaneResult(damage=float(damage[0]), metric=float(metric[0]))

    def evaluate_on_planes(
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
        tau_max=None,
        strain_ranges=None,
    ):
        check_planes(S0r, E0r, S1r, E1r)
        material = _validated_constants(params)
        gamma_a = np.sqrt(pypow2(E0r[:, 0, 2] - E1r[:, 0, 2]) + pypow2(E0r[:, 1, 2] - E1r[:, 1, 2]))
        tau_a = 0.5 * np.sqrt(
            pypow2(S0r[:, 0, 2] - S1r[:, 0, 2]) + pypow2(S0r[:, 1, 2] - S1r[:, 1, 2])
        )
        sigma_n_max = np.maximum(S0r[:, 2, 2], S1r[:, 2, 2])
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = np.where(
                tau_a > 0.0,
                sigma_n_max / np.where(tau_a > 0.0, tau_a, 1.0),
                0.0,
            )
            Nf = mmccm_life(gamma_a, rho, material)
            damage = np.where(gamma_a > 0.0, 1.0 / Nf, 0.0)
        return damage, gamma_a
