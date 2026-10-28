"""Fatigue index error for fatigue-limit datasets (decision d057).

A fatigue-limit dataset holds limit states: each plain row is a load sitting at the fatigue
limit of its cycle basis (recorded as ``Run out = Yes``) and carries no failure lives, hence no
master curve can be built. Every model is then assessed through its fatigue index, following
Boehme & Papuga (2024), FFEMS 47, doi:10.1111/ffe.14281, p. 7: "FI_est = sigma_eq / s_1"
(Eq. 8) and the fatigue index error "Delta FI = FI_est / FI_ref - 1" (Eq. 9), with FI_ref = 1 for
a limit state. For any damage parameter DP the reference is the model's own DP at the axial
R = -1 fatigue limit:

    FI_axial = DP(case) / DP(axial, sigma_a = Sigm1) - 1

which is Eqs. (8)-(9) for a criterion whose DP is an equivalent stress equal to s_1 at that
load. FI_torsion refers instead to the torsional R = -1 limit Taum1. FI_load = 1/s - 1, where s
is the factor applied to every stress and strain (means included) at which
DP(s * case) = DP(axial reference), is the stress-scaled counterpart valid for any DP: it
coincides with FI_axial for a criterion linear in the loads and equals sqrt(DP/DP_ref) - 1 for a
quadratic one such as SWT (decision d063). A missing plane constant (k_FI, k_FS, ...) is chosen
so that both references return the same DP, which is the customary fatigue-limit calibration of
such criteria.
"""

from __future__ import annotations

import numpy as np

from multifat.imp.experimental_data import build_stress_tensor, stress_to_strain_iso

# Plane constants a fatigue-limit calibration may set, with the search bounds used for it
ENDURANCE_CONSTANTS = {
    "k_FS": (0.0, 3.0),
    "k_FI": (0.0, 3.0),
    "L_LI": (0.0, 5.0),
    "k_MGSE": (0.0, 5.0),
    "a_MSWT": (0.0, 1.0),
}


def is_fatigue_limit_dataset(block: dict) -> bool:
    """True when nothing failed outright, some plain rows ran out, and both R = -1 limits exist."""
    material = block["material"]
    limits = [
        material.get("Tensile_fatigue_limit_MPa"),
        material.get("Torsional_fatigue_limit_MPa"),
    ]
    flagged = block.get("property_issues", {})
    return (
        not block["cases"]
        and bool(block.get("runouts"))
        and all(v is not None and v > 0 for v in limits)
        and "Sigm1" not in flagged
        and "Taum1" not in flagged
    )


def case_from_amplitudes(*, axial=0.0, shear=0.0, radial=0.0, E=None, nu=None, mean=(0, 0, 0)):
    """In-phase R = -1 loading (means optional), in the tensor layout expected by the importer."""
    sig_mean = build_stress_tensor(*mean)
    sig_amp = build_stress_tensor(radial, axial, shear)
    zeros = np.zeros((3, 3))
    parts = {
        tag: stress_to_strain_iso(t, E, nu) for tag, t in (("mean", sig_mean), ("sin", sig_amp))
    }
    return {
        "S0": sig_mean + sig_amp,
        "S1": sig_mean - sig_amp,
        "E0": parts["mean"] + parts["sin"],
        "E1": parts["mean"] - parts["sin"],
        "H": {"S_mean": sig_mean, "S_sin": sig_amp, "S_cos": zeros, "E_mean": parts["mean"],
              "E_sin": parts["sin"], "E_cos": stress_to_strain_iso(zeros, E, nu)},
        "meta": {},
    }  # fmt: skip


def reference_cases(Sigm1: float, Taum1: float, E, nu) -> dict[str, dict]:
    """Axial and torsional R = -1 loads placed at the fatigue limits."""
    return {
        "axial": case_from_amplitudes(axial=Sigm1, E=E, nu=nu),
        "torsion": case_from_amplitudes(shear=Taum1, E=E, nu=nu),
    }


def case_from_block(case: dict) -> dict:
    """An imported case (failure or run-out) converted to the layout used here."""
    return {
        "S0": case[1][1]["stress"],
        "E0": case[1][1]["strain"],
        "S1": case[1][2]["stress"],
        "E1": case[1][2]["strain"],
        "H": case["harmonics"],
        "meta": case["meta"],
    }


def damage(model, params: dict, case: dict, R_list, *, variant: str = "") -> float:
    result = model.evaluate_case(
        S0=case["S0"], E0=case["E0"], S1=case["S1"], E1=case["E1"], params=params,
        R_list=R_list, harmonics=case["H"],
    )  # fmt: skip
    return float(result.values.get(variant, np.nan))


def _root_in_interval(f, lo: float, hi: float, f_lo: float, tol: float = 1e-13) -> float:
    """Bisection for the root of f inside [lo, hi], which must bracket a sign change (relative `tol`)."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0 or (hi - lo) <= tol * max(1.0, abs(mid)):
            return mid
        if (f_mid > 0) == (f_lo > 0):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def calibrate_endurance_constant(
    model, key: str, params: dict, refs: dict, bounds, R_list, *, variant: str = ""
):
    """Value of `key` inside `bounds` making both references agree on DP, or None."""

    def mismatch(p: float) -> float:
        trial = {**params, key: float(p)}
        return damage(model, trial, refs["axial"], R_list, variant=variant) - damage(
            model, trial, refs["torsion"], R_list, variant=variant
        )

    lo, hi = float(bounds[0]), float(bounds[1])
    samples = np.linspace(lo, hi, 31)
    values = np.array([mismatch(p) for p in samples])
    for a, b, fa, fb in zip(samples[:-1], samples[1:], values[:-1], values[1:]):
        if fa == 0.0:
            return float(a)
        if np.isfinite(fa) and np.isfinite(fb) and fa * fb < 0:
            return float(_root_in_interval(mismatch, a, b, fa))
    return float(samples[-1]) if values[-1] == 0.0 else None


def scaled_case(case: dict, factor: float) -> dict:
    """All stresses and strains (means included) times `factor`; being elastic, strains follow stresses."""
    return {
        "S0": factor * case["S0"],
        "S1": factor * case["S1"],
        "E0": factor * case["E0"],
        "E1": factor * case["E1"],
        "H": {k: factor * v for k, v in case["H"].items()},
        "meta": case.get("meta", {}),
    }


def load_factor(
    model, params: dict, case: dict, target: float, R_list, *, variant: str = ""
) -> float:
    """Factor s for which DP(s * case) = target (DP grows with the loads), or NaN."""

    def residual(s: float) -> float:
        return damage(model, params, scaled_case(case, s), R_list, variant=variant) - target

    r1 = residual(1.0)
    if not np.isfinite(r1) or not np.isfinite(target) or target <= 0:
        return float("nan")
    if r1 == 0.0:
        return 1.0
    lo, hi = (1.0, 2.0) if r1 < 0 else (0.5, 1.0)
    for _ in range(60):  # grow or shrink the bracket until it straddles the root
        r_lo, r_hi = residual(lo), residual(hi)
        if np.isfinite(r_lo) and np.isfinite(r_hi) and r_lo <= 0.0 <= r_hi:
            break
        lo, hi = (hi, 2.0 * hi) if r1 < 0 else (0.5 * lo, lo)
    else:
        return float("nan")
    for _ in range(200):  # bisection carried out in log scale
        mid = float(np.sqrt(lo * hi))
        r_mid = residual(mid)
        if r_mid == 0.0 or hi / lo - 1.0 < 1e-13:
            return mid
        lo, hi = (mid, hi) if r_mid < 0.0 else (lo, mid)
    return float(np.sqrt(lo * hi))


def fatigue_index(
    model, params: dict, cases: list[dict], refs: dict, R_list, *, variant: str = ""
) -> list[dict]:
    """DP, the two reference DPs and both fatigue index errors, for every case."""
    dp_axial_ref = damage(model, params, refs["axial"], R_list, variant=variant)
    dp_torsion_ref = damage(model, params, refs["torsion"], R_list, variant=variant)
    rows = []
    for case in cases:
        dp = damage(model, params, case, R_list, variant=variant)
        rows.append(
            {
                "DP": dp,
                "DP_ref_axial": dp_axial_ref,
                "DP_ref_torsion": dp_torsion_ref,
                "FI_axial": dp / dp_axial_ref - 1.0,
                "FI_torsion": dp / dp_torsion_ref - 1.0,
                "FI_load": 1.0
                / load_factor(model, params, case, dp_axial_ref, R_list, variant=variant)
                - 1.0,
            }
        )
    return rows


def fatigue_index_summary(df, material: str):
    """Per-method statistics of the fatigue index errors, taken over the finite values only.

    A model may be infinite on some loads by design (MWCM beyond its rho cut-off): such rows are
    counted in ``n_nonfinite`` and kept out of the statistics.
    """
    import pandas as pd

    records = []
    for method, group in df.groupby("Method"):
        fine = group[np.isfinite(group.FI_axial) & np.isfinite(group.FI_torsion)]
        load_err = (
            group.FI_load[np.isfinite(group.FI_load)]
            if "FI_load" in group
            else pd.Series(dtype=float)
        )
        records.append(
            {
                "Method": method,
                "Material": material,
                "n": len(group),
                "n_nonfinite": len(group) - len(fine),
                "FI_axial_mean": fine.FI_axial.mean(),
                "FI_axial_sd": fine.FI_axial.std(),
                "FI_axial_min": fine.FI_axial.min(),
                "FI_axial_max": fine.FI_axial.max(),
                "FI_torsion_mean": fine.FI_torsion.mean(),
                "FI_load_mean": load_err.mean(),
                "FI_load_sd": load_err.std(),
            }
        )
    return pd.DataFrame(records)
