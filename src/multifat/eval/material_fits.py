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
# Out-of-phase Woehler curve sigma_a = sigma_f_OP (2N)^b_OP of one tension-torsion path
# tau_a = lambda_OP sigma_a at phase delta_OP (deg), the calibration curve of BOLCHOUN
OP_CURVE_KEYS = ("sigma_f_OP", "b_OP", "lambda_OP", "delta_OP")
OP_ZERO_RTOL = 1e-6  # components below this fraction of the largest amplitude are zero


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


def out_of_phase_tests(data: dict, material: str) -> dict:
    """Fully reversed tension-torsion failures out of phase, grouped by path.

    A test qualifies when its mean stress is zero, its only components are the axial normal
    stress and the axial-hoop shear, both are non-zero and their phase difference is not a
    multiple of 180 deg. The path key is (tau_a/sigma_a rounded to 0.01, phase rounded to 1 deg);
    the value lists (case ids, axial amplitudes, lives, tau_a/sigma_a, phases in deg).
    """
    groups: dict = {}
    for case_id, case in data[material]["cases"].items():
        h = case.get("harmonics") or {}
        try:
            S_mean, S_sin, S_cos = (
                np.asarray(h[k], dtype=float) for k in ("S_mean", "S_sin", "S_cos")
            )
        except KeyError:
            continue
        scale = max(np.max(np.abs(S_sin)), np.max(np.abs(S_cos)), 1e-300)
        other = np.abs(np.stack([S_sin, S_cos]))
        other[:, 1, 1] = other[:, 0, 1] = other[:, 1, 0] = 0.0
        if np.max(np.abs(S_mean)) > OP_ZERO_RTOL * scale or np.max(other) > OP_ZERO_RTOL * scale:
            continue
        sigma = complex(S_sin[1, 1], S_cos[1, 1])
        tau = complex(S_sin[0, 1], S_cos[0, 1])
        if abs(sigma) <= OP_ZERO_RTOL * scale or abs(tau) <= OP_ZERO_RTOL * scale:
            continue
        delta = float(np.degrees(np.angle(tau / sigma))) % 180.0
        if min(delta, 180.0 - delta) < 0.5:
            continue  # proportional (in phase or in counter-phase)
        key = (round(abs(tau) / abs(sigma), 2), round(delta))
        ids, amps, lives, ratios, phases = groups.setdefault(key, ([], [], [], [], []))
        ids.append(int(case_id))
        amps.append(abs(sigma))
        lives.append(float(case["meta"]["Nf_cycles"]))
        ratios.append(abs(tau) / abs(sigma))
        phases.append(delta)
    return groups


def out_of_phase_fill(data: dict, material: str) -> dict:
    """Woehler curve of the most populated out-of-phase path (OP_CURVE_KEYS), if identifiable.

    Ties between paths go to the one whose lives span more decades. The curve is retained under
    the rule of ``fit_basquin`` (5 tests, 1 decade, slope distinct from 0); its tests are then
    design points of the models that use it.
    """
    groups = out_of_phase_tests(data, material)
    if not groups:
        return {"params": {}, "fit": None, "design_case_ids": [], "reason": "no out-of-phase "
                "fully reversed tension-torsion failures"}  # fmt: skip

    def rank(item):
        _, (ids, _, lives, _, _) = item
        return (len(ids), float(np.log10(max(lives) / min(lives))))

    _, (ids, amps, lives, ratios, phases) = max(groups.items(), key=rank)
    # the path itself: the mean of the group (the rounded key only groups the tests)
    lambda_op, delta = float(np.mean(ratios)), float(np.mean(phases))
    fit = fit_basquin(amps, lives)
    fit = {**fit, "case_ids": ids, "lambda_OP": lambda_op, "delta_OP": delta}
    if not fit["identifiable"]:
        return {"params": {}, "fit": fit, "design_case_ids": [], "reason": fit["reason"]}
    params = {
        "sigma_f_OP": fit["coef"],
        "b_OP": fit["exp"],
        "lambda_OP": lambda_op,
        "delta_OP": delta,
    }
    return {"params": params, "fit": fit, "design_case_ids": ids, "reason": ""}


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
