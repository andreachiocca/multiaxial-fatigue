"""Load-equivalent error of a fatigue assessment (issue #37).

For a failed test with observed life N_exp, the load-equivalent error is::

    e_L = ln(1 / s*)

where s* scales the whole load history (means included; strains follow elastically) so that the
method predicts exactly N_exp. With a master curve DP = A N^b, s* solves
DP(s* sigma) = A N_exp^b; with a native life law, N_native(s* sigma) = N_exp. e_L > 0 means the
method rates the test as more damaging than it was (conservative), e_L < 0 non-conservative.

Unlike ln(N_exp / N_pred), e_L does not grow with the flatness k = -1/b of the S-N curve
(SD ln N ~ k SD e_L), and unlike log10(DP_fit / DP) it does not grow with the degree d of a DP
homogeneous in the loads (SWT and the energy parameters have d = 2). It is the life-data
counterpart of the stress-scaled fatigue index error FI_load of fatigue-limit datasets
(``multifat.eval.endurance``).

Repository choices:
  - homogeneous DPs (DP(2 x) = 2^d DP(x) and DP(x/2) = 2^-d DP(x) to 1e-9) use the closed form
    e_L = ln(DP / DP_fit) / d; other DPs and native life laws are solved by a secant iteration in
    ln s, with the bisection of ``endurance.load_factor`` as fall-back;
  - the secant steps are confined to |ln s| <= LOG_S_MAX = 10; a diverging iteration falls back
    to the bisection (issue #71: SSF extrapolated beyond its fit range on FABEST);
  - run-outs give bounds only and are not evaluated.
"""

from __future__ import annotations

import math

import numpy as np

from multifat.models.errors import LoadingScopeError

from multifat.eval.endurance import load_factor, scaled_case, scaled_damage

HOMOGENEITY_TOL = 1e-9


def homogeneity_degree(model, params: dict, case: dict, R_list, DP: float, *, variant: str = ""):
    """Degree d with DP(s x) = s^d DP(x), checked at s = 2 and 1/2, or None."""
    if not (np.isfinite(DP) and DP > 0):
        return None
    up = scaled_damage(model, params, case, 2.0, R_list, variant=variant)
    down = scaled_damage(model, params, case, 0.5, R_list, variant=variant)
    if not (np.isfinite(up) and np.isfinite(down) and up > 0 and down > 0):
        return None
    d_up, d_down = math.log2(up / DP), math.log2(DP / down)
    if d_up <= 0 or abs(d_up - d_down) > HOMOGENEITY_TOL * max(1.0, abs(d_up)):
        return None
    return d_up


# Secant steps stay within |ln s| <= LOG_S_MAX (load factors 4.5e-5 to 2.2e4); a step beyond it
# means the iteration diverges (e.g. a DP that is not monotonic in the load, such as an
# extrapolated SSF surface), and the bisection of endurance.load_factor takes over.
LOG_S_MAX = 10.0


def _solve_in_log(g, s0: float = 1.0, tol: float = 1e-12, max_iter: int = 60):
    """Root of the increasing function g(ln s) by secant steps from ln s0; None if it fails."""
    x0, x1 = math.log(s0), math.log(s0) + 0.05
    g0, g1 = g(x0), g(x1)
    for _ in range(max_iter):
        if not (np.isfinite(g0) and np.isfinite(g1)) or g1 == g0:
            return None
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        if not abs(x2) <= LOG_S_MAX:
            return None
        if abs(x2 - x1) < tol:
            return math.exp(x2)
        x0, g0, x1, g1 = x1, g1, x2, g(x2)
    return None


def master_curve_load_error(
    model, params: dict, case: dict, R_list, DP: float, DP_fit: float, *, variant: str = ""
) -> float:
    """e_L for a DP compared with the master-curve value DP_fit = A N_exp^b."""
    if not (np.isfinite(DP) and np.isfinite(DP_fit) and DP > 0 and DP_fit > 0):
        return float("nan")
    degree = homogeneity_degree(model, params, case, R_list, DP, variant=variant)
    if degree is not None:
        # s*^d DP = DP_fit, so e_L = -ln s* = ln(DP / DP_fit) / d
        return math.log(DP / DP_fit) / degree

    def g(x: float) -> float:
        value = scaled_damage(model, params, case, math.exp(x), R_list, variant=variant)
        return math.log(value) - math.log(DP_fit) if value > 0 else -np.inf

    s = _solve_in_log(g)
    if s is None:
        s = load_factor(model, params, case, DP_fit, R_list, variant=variant)
    return -math.log(s) if np.isfinite(s) and s > 0 else float("nan")


def native_life_load_error(
    model, params: dict, case: dict, R_list, fit, Nf: float, *, variant: str = ""
):
    """e_L for a native life law: N_native(s* x) = N_exp, by bisection in ln s.

    NaN when undefined, including when a trial load leaves the model's declared scope
    (LoadingScopeError): the native life is then not available at that load.
    """
    if not (np.isfinite(Nf) and Nf > 0):
        return float("nan")
    try:
        return _native_life_load_error(model, params, case, R_list, fit, Nf, variant)
    except LoadingScopeError:
        return float("nan")


def _native_life_load_error(model, params, case, R_list, fit, Nf, variant):
    """Bracketing and bisection of ``native_life_load_error``."""

    def log_life(x: float) -> float:
        scaled = scaled_case(case, math.exp(x))
        result = model.evaluate_case(
            S0=scaled["S0"], E0=scaled["E0"], S1=scaled["S1"], E1=scaled["E1"], params=params,
            R_list=R_list, harmonics=scaled["H"],
        )  # fmt: skip
        aux = {"dp": float(result.values.get(variant, np.nan)), **result.aux}
        Nf_scaled = float(np.asarray(model.native_life(aux, params, fit), dtype=float))
        if math.isnan(Nf_scaled):
            raise _UndefinedLife  # no life at this load: the search fails (#78)
        return math.log(Nf_scaled) if Nf_scaled > 0 else -np.inf

    target = math.log(Nf)
    try:
        # The life decreases with the load: find lo (life >= N_exp) and hi (life < N_exp).
        lo, hi = -0.1, 0.1
        for _ in range(40):
            if log_life(lo) >= target:
                break
            lo -= 0.5
        else:
            return float("nan")
        for _ in range(40):
            if log_life(hi) < target:
                break
            hi += 0.5
        else:
            return float("nan")
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if hi - lo < 1e-12:
                break
            lo, hi = (mid, hi) if log_life(mid) >= target else (lo, mid)
    except _UndefinedLife:
        return float("nan")
    return -0.5 * (lo + hi)


class _UndefinedLife(Exception):
    """A NaN native life inside the load search: no defined life at that trial load."""
