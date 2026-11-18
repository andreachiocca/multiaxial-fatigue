"""Material constants derived from the workbook's own test data (Phase 3, option D of the plan).

Basquin relations sigma_a = sigma'_f (2N)^b for fully reversed uniaxial tension tests and
tau_a = tau'_f (2N)^b0 for fully reversed torsion tests (log10 least squares, taking life as
the independent variable). A relation is retained only if it can be identified: at least 5
tests, at least one decade of life, and a slope whose distance from zero exceeds twice its
standard error. Plastic strain-life terms are not identifiable from elastic HCF data;
``strain_life_fill`` takes them from the uniform material law (steel form) and always marks
them as estimates.
"""

from __future__ import annotations

import math

import numpy as np

from multifat.eval.fitting import MIN_DECADES, MIN_POINTS, MIN_T
from multifat.models.utils.fatigue_estimation import (
    estimate_axial_fatigue_props_uniform_material_law,
)

BASQUIN_KEYS = ("sigma_f", "b_fat", "tau_f", "b0")
PLASTIC_KEYS = ("epsilon_f", "c_fat", "gamma_f", "c0")


def fit_basquin(amplitude, N) -> dict:
    """Fit amplitude = coef (2N)^exp and report whether the fit is identifiable."""
    amplitude = np.asarray(amplitude, dtype=float)
    N = np.asarray(N, dtype=float)
    valid = np.isfinite(amplitude) & np.isfinite(N) & (amplitude > 0) & (N > 0)
    amplitude, N = amplitude[valid], N[valid]
    fit = {
        "n": int(amplitude.size),
        "coef": None,
        "exp": None,
        "se_exp": None,
        "decades": 0.0,
    }
    if amplitude.size >= 2:
        fit["decades"] = float(np.log10(N.max() / N.min()))
    if amplitude.size < 3 or fit["decades"] == 0.0:
        return {**fit, "identifiable": False, "reason": f"{amplitude.size} tests"}
    log_N, log_amp = np.log10(2.0 * N), np.log10(amplitude)
    slope, intercept = np.polyfit(log_N, log_amp, 1)
    residuals = log_amp - (slope * log_N + intercept)
    slope_se = math.sqrt(
        float(residuals @ residuals)
        / (amplitude.size - 2)
        / float(((log_N - log_N.mean()) ** 2).sum())
    )
    fit.update(coef=float(10.0**intercept), exp=float(slope), se_exp=float(slope_se))
    issues = []
    if amplitude.size < MIN_POINTS:
        issues.append(f"{amplitude.size} < {MIN_POINTS} tests")
    if fit["decades"] < MIN_DECADES:
        issues.append(f"{fit['decades']:.2f} < {MIN_DECADES} decades of life")
    if not slope < 0.0 or abs(slope) <= MIN_T * slope_se:
        issues.append(f"slope {slope:.4f} not distinct from 0 (SE {slope_se:.4f})")
    return {**fit, "identifiable": not issues, "reason": "; ".join(issues)}


def reversed_uniaxial_tests(data: dict, material: str) -> dict:
    """Case ids, amplitudes and lives of the fully reversed tension and torsion tests.

    Membership is ``imp.experimental_data.is_design_test``, the selector of the master curve
    and of the plane-constant calibration (#77); the level is the harmonic amplitude.
    """
    from multifat.imp.experimental_data import is_design_test

    collected = {"tension": ([], [], []), "torsion": ([], [], [])}
    for case_id, case in data[material]["cases"].items():
        harmonics = case["harmonics"]
        amp_arr = np.hypot(harmonics["S_sin"], harmonics["S_cos"])
        if is_design_test(case, "tension"):
            family, level = "tension", amp_arr[1, 1]
        elif is_design_test(case, "torsion"):
            family, level = "torsion", amp_arr[0, 1]
        else:
            continue
        id_list, amp_list, N_list = collected[family]
        id_list.append(int(case_id))
        amp_list.append(float(level))
        N_list.append(float(case["meta"]["Nf_cycles"]))
    return collected


def basquin_fill(data: dict, material: str) -> dict:
    """Fitted sigma_f, b_fat, tau_f, b0 (identifiable curves only) together with the evidence."""
    test_sets = reversed_uniaxial_tests(data, material)
    curve_fits, coefficients, design_ids = {}, {}, {}
    for family, keys in (
        ("tension", ("sigma_f", "b_fat")),
        ("torsion", ("tau_f", "b0")),
    ):
        id_list, amp_list, N_list = test_sets[family]
        fit = fit_basquin(amp_list, N_list)
        curve_fits[family] = {**fit, "case_ids": id_list}
        if fit["identifiable"]:
            coefficients[keys[0]], coefficients[keys[1]] = fit["coef"], fit["exp"]
            design_ids[family] = id_list
    return {
        "params": coefficients,
        "fits": curve_fits,
        "design_case_ids": design_ids,
    }


def strain_life_fill(Su: float | None, E: float | None) -> dict:
    """Plastic strain-life terms from the uniform material law (steels); estimates only."""
    if Su is None or E is None:
        return {}
    props = estimate_axial_fatigue_props_uniform_material_law(Su=float(Su), E=float(E))
    return {
        "epsilon_f": props.epsilon_f,
        "c_fat": props.c,
        "gamma_f": math.sqrt(3.0) * props.epsilon_f,
        "c0": props.c,
    }


def zhu_k_mgse(sf, b, tf, b0, E, nu, N_lo=1e4, N_hi=2e6, n=200) -> float:
    """Zhu et al. (2018) Eq. (12), elastic contributions only (repo elastic-strain convention),
    averaged over log-spaced lives N_lo..N_hi (the paper recommends the mean over a life range):
    k(N) = [tf^2/G (2N)^(2 b0) - (1 + nu)/2 sf^2/E (2N)^(2b)] / [sf^2/(2E) (2N)^(2b)].
    """
    G_shear = E / (2.0 * (1.0 + nu))
    two_n = 2.0 * np.logspace(np.log10(N_lo), np.log10(N_hi), n)
    shear = tf**2 / G_shear * two_n ** (2.0 * b0)
    normal = sf**2 / E * two_n ** (2.0 * b)
    return float(np.mean((shear - 0.5 * (1.0 + nu) * normal) / (0.5 * normal)))
