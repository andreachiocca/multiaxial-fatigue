"""Calibration of the Tveit et al. (2025) model TVEIT2025: A, then C, K and L by Eq. (25).

Sources: Tveit, Ottosen, Erduran & Reyes (2025), Mech. Mater. 209, 105452,
doi:10.1016/j.mechmat.2025.105452, Sect. 4 (calibration; A from fatigue limits with different
mean stresses, Eqs. (22)-(23), or estimated with the FKM guideline, Eq. (26), or the Goodman
relation, Table 1; C, K and L from Eq. (25)); and Ottosen, Stenstroem & Ristinmaa (2008),
Int. J. Fatigue 30, 996-1006, doi:10.1016/j.ijfatigue.2007.08.009, Sect. 5 (closed form of the
periodic uniaxial cycle, Eqs. (23), (24), (28), (29); simplex minimisation, Eqs. (30)-(32)).

The model is TVEIT2025 (multifat.models.invariant_methods.tveit2025). For a uniaxial cycle between
sigma_2 = sigma_m + sigma_a and sigma_4 = sigma_m - sigma_a, the periodic state of the endurance
surface (back-stress alpha_2 at the maximum, alpha_4 at the minimum) solves, with k = C A / f_m1:

    1.5 alpha_2 - (A + 1) sigma_2 + f_m1
        - (A + 1)/k ln[(1 - k (sigma_2 - 1.5 alpha_2)) / (1 - k/(A + 1) (f_m1 - 1.5 A alpha_4))] = 0   (23)
    -1.5 alpha_4 - (A - 1) sigma_4 + f_m1
        - (A - 1)/k ln[(1 - k (sigma_4 - 1.5 alpha_4)) / (1 - k/(A - 1) (f_m1 - 1.5 A alpha_2))] = 0   (24)
    1/N = K/L [exp(L beta_2) + exp(L beta_4) - 2]                                          (28)-(29)
    beta_2 = ((A + 1) sigma_2 - 1.5 alpha_2 - f_m1)/f_m1,  beta_4 = ((A - 1) sigma_4 + 1.5 alpha_4 - f_m1)/f_m1
    E(C, K, L) = sum_i (1 - ln N_i / ln N_exp,i)^2                       Tveit Eq. (25), Ottosen Eq. (30)

and A follows from a fatigue limit (sigma_m0, sigma_a0) with sigma_m0 > 0 on the uniaxial Haigh
line sigma_a + A sigma_m = f_m1, A = (f_m1 - sigma_a0) / sigma_m0 (Tveit Eqs. 22-23), or is
estimated as A = a_M Su / 1000 + b_M (FKM guideline, Tveit Eq. 26) or A = f_m1 / Su (Goodman).

Repository choices:
  - The procedure is the authors' one ("author" mode of docs/calibration-validation-protocol.md);
    it runs on any fatigue-life workbook. f_m1, t_m1 and Su are the workbook's own Summary values
    (the R = -1 rows), and the calibration tests are its axial R = -1 failures, selected with the
    shared selector multifat.imp.experimental_data.is_design_test (issue #77).
  - What a workbook cannot hold comes from an optional hand-written file,
    Calibrations/<material>/TVEIT2025_inputs.json: a fatigue limit at a tensile mean stress, the
    material type for the FKM estimate, a hand-entered A, the source of the calibration tests and
    the time-integration settings. The calibration never writes to it; the record keeps a copy
    and its SHA-256.
  - A follows the order of Tveit et al. (2025, Sect. 4): a fatigue limit at a tensile mean stress
    (the most reliable route); else a hand-entered A (for example fitted to the life of a
    mean-stress test); else the FKM estimate for a declared material type with tabulated a_M, b_M
    (FKM_MEAN_STRESS); else Goodman, which the authors use where the FKM guideline gives no
    values.
  - One source: the calibration tests are those of the source of the Summary fatigue limits (same
    DOI, else same authors), unless the inputs name another; without a match every source is used
    and the confidence says so.
  - Three points: the lowest (LOW) and highest (HIGH) stress level above f_m1 and the level closest
    in ln sigma_a to sqrt(LOW HIGH) (MID; a tie goes to the lower level); at each level the test
    whose log10 N is closest to the median of the level (a tie goes to the lowest Excel row). The
    papers prescribe three points of one S-N curve, not how they are chosen.
  - The confidence of the calibration (multifat.eval.calibration_confidence) lists what weakens
    it: failures below f_m1, limits and tests of different sources, an estimated A, points close
    together or near f_m1, single-specimen levels, a poor fit, a constant on its search bound and
    the unidentified C. Each finding says what data would remove it.
  - The closed form is solved by a damped Newton method. The physical root is the one connected
    by continuation to the endurance limit sigma_a = f_m1 / (1 + A sigma_m / sigma_a), where it
    is alpha_2 = ((A + 1) sigma_2 - f_m1) / 1.5, alpha_4 = (f_m1 - (A - 1) sigma_4) / 1.5 exactly,
    along the ray of constant sigma_m / sigma_a (40 geometric steps, a step is halved when the
    corrector strays from the predictor). A root must keep the relative stress of the load sign at
    both extremes. Below the endurance limit the life is infinite.
  - Uniaxial S-N data do not identify C: many (C, K, L) sets fit almost equally well, with C
    over three to four decades, and a multi-start search ends at an arbitrary member of that
    flat valley. Eq. (25) is therefore profiled over C: at each C of a fixed grid (8 values per
    decade, C in [1e-4, 1e3]) K and L are fitted by the simplex method of Ottosen et al. (2008)
    in ln K, ln L (K in [1e-12, 1e-1], L in [0.05, 50]); a non-finite or N <= 1 prediction costs
    10 in Eq. (25). Tveit et al. (2025) used MATLAB fmincon instead.
  - The grid points within E <= 1.05 E_min form the equivalent family; its C range is recorded
    as the identifiability diagnostic. A set is admissible when the closed-form root is unique
    at every calibration point (G1) and the time-integrated life of the model, with the
    numerical settings of the record, agrees with the closed form within 0.10 decade (G2),
    since the closed form ignores the transient. The chosen set is the admissible grid point of
    lowest E; a chosen C on a grid end is flagged as set by the search bound.
  - The constants hold together with the time-integration settings G2 used, so the record
    writes them (npts_TVEIT ... tol_damage_TVEIT; model defaults unless the inputs file sets
    them) next to A, C, K, L and the fatigue limits.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from multifat.eval.calibration_confidence import MAJOR, MINOR, Finding, report_lines, summary

# Version of this calibration procedure (point rules, objective, optimiser, family, gates, routes
# to A). Increment it when one of them changes.
PROCEDURE = 2
MODE = "author"
INPUTS_SCHEMA = "multifat/tveit2025-calibration-inputs/2"
INPUTS_FILE = "TVEIT2025_inputs.json"
INPUTS_KEYS = {
    "schema",
    "material_id",
    "material_type",
    "mean_stress_fatigue_limit",
    "A",
    "calibration_source",
    "numerical_settings",
    "notes",
}
CALIBRATIONS_DIR = Path(__file__).resolve().parents[3] / "Calibrations"

# a_M, b_M of the FKM guideline (Haenel et al. 2003) for A = a_M Su/1000 + b_M (Tveit Eq. 26), by
# material type. Extend only with values checked in the guideline.
FKM_MEAN_STRESS = {"steel": (0.35, -0.1)}

# Thresholds of the confidence findings
NARROW_SPAN = 1.10  # HIGH / LOW stress below this: the points carry little S-N trend
NEAR_LIMIT = 1.05  # LOW / f_m1 below this: the LOW point sits on the asymptote
POOR_FIT_DECADES = 0.3  # residual at a calibration point
C_SPAN_DECADES = 1.0  # C range of the equivalent sets

BOUNDS = {"C": (1e-4, 1e3), "K": (1e-12, 1e-1), "L": (0.05, 50.0)}
C_PER_DECADE = 8  # grid of the profile over C
MAX_EVALUATIONS = 2000  # objective evaluations per simplex start
PENALTY = 10.0  # Eq. (25) term of a non-finite or N <= 1 prediction
EQUIVALENT = 1.05  # equivalent family: E <= 1.05 E_min
EXACT_FIT = 1e-6  # below this E_min the family is the sets within 0.01 decade at every point
G2_DECADES = 0.10  # transient against closed form at each calibration point
CONTINUATION_STEPS = 40
START_OFFSET = 1e-3  # the branch starts at (1 + 1e-3) times the endurance amplitude

VALID = "valid"
BELOW = "below_endurance_limit"
AMBIGUOUS = "ambiguous_root"
FOLD = "out_of_validity_fold"
NO_START = "out_of_validity_no_start"


class CalibrationInputError(ValueError):
    """The calibration inputs are missing, inconsistent or do not allow the requested mode."""


# ---------------------------------------------------------------------------------------------
# Closed form of the periodic uniaxial cycle (Ottosen et al. 2008, Sect. 5)
# ---------------------------------------------------------------------------------------------


def _system(a2, a4, C, A, f_m1, s2, s4):
    """Residuals of Eqs. (23)-(24) and their Jacobian; None outside the domain of the logarithms."""
    k = C * A / f_m1
    num1 = 1.0 - k * (s2 - 1.5 * a2)
    den1 = 1.0 - k / (A + 1.0) * (f_m1 - 1.5 * A * a4)
    num2 = 1.0 - k * (s4 - 1.5 * a4)
    den2 = 1.0 - k / (A - 1.0) * (f_m1 - 1.5 * A * a2)
    if num1 <= 0.0 or den1 <= 0.0 or num2 <= 0.0 or den2 <= 0.0:
        return None
    f1 = 1.5 * a2 - (A + 1.0) * s2 + f_m1 - (A + 1.0) / k * math.log(num1 / den1)
    f2 = -1.5 * a4 - (A - 1.0) * s4 + f_m1 - (A - 1.0) / k * math.log(num2 / den2)
    j11 = 1.5 - 1.5 * (A + 1.0) / num1
    j12 = 1.5 * A / den1
    j21 = 1.5 * A / den2
    j22 = -1.5 - 1.5 * (A - 1.0) / num2
    return f1, f2, j11, j12, j21, j22


def _newton(a2, a4, C, A, f_m1, s2, s4, max_iter=60):
    """Damped Newton iteration from (a2, a4); returns (a2, a4, max residual) or None."""
    system = _system(a2, a4, C, A, f_m1, s2, s4)
    if system is None:
        return None
    target = 1e-11 * max(1.0, f_m1)
    for _ in range(max_iter):
        f1, f2, j11, j12, j21, j22 = system
        norm = max(abs(f1), abs(f2))
        if norm < target:
            break
        det = j11 * j22 - j12 * j21
        if det == 0.0 or not math.isfinite(det):
            break
        d2 = (-f1 * j22 + f2 * j12) / det
        d4 = (f1 * j21 - f2 * j11) / det
        step = 1.0
        for _ in range(40):
            trial = _system(a2 + step * d2, a4 + step * d4, C, A, f_m1, s2, s4)
            if trial is not None and max(abs(trial[0]), abs(trial[1])) < norm:
                a2, a4, system = a2 + step * d2, a4 + step * d4, trial
                break
            step *= 0.5
        else:
            break
    f1, f2 = system[0], system[1]
    return a2, a4, max(abs(f1), abs(f2))


def _accepted(root, f_m1, s2, s4) -> bool:
    """Converged root whose relative stress has the sign of the load at both extremes."""
    if root is None:
        return False
    a2, a4, residual = root
    return residual < 1e-6 * max(1.0, f_m1) and s2 - 1.5 * a2 > 0.0 and 1.5 * a4 - s4 > 0.0


def endurance_amplitude(A: float, f_m1: float, sigma_m: float, sigma_a: float) -> float:
    """Endurance amplitude on the ray of constant sigma_m / sigma_a (Haigh line, Eq. 22)."""
    denominator = 1.0 + A * sigma_m / sigma_a
    return f_m1 / denominator if denominator > 0.0 else math.nan


def _limit_root(A, f_m1, s2, s4):
    """Exact root of Eqs. (23)-(24) on the endurance limit: both logarithms vanish there."""
    return ((A + 1.0) * s2 - f_m1) / 1.5, (f_m1 - (A - 1.0) * s4) / 1.5


def _follow_branch(C, A, f_m1, slope, start, x_start, target, n_steps):
    """Continuation of the root from amplitude `start` to `target` along sigma_m = slope sigma_a.

    A step is accepted when the corrector converges, keeps the load sign of the relative stress
    and stays within max(0.2 |x_pred - x_prev|, 0.05 MPa) of the secant predictor; otherwise it is
    halved. Returns ((alpha_2, alpha_4), folded).
    """
    limit = f_m1 / (1.0 + A * slope)
    d0, d1 = start - limit, target - limit
    grid = [limit + d0 * (d1 / d0) ** (i / n_steps) for i in range(1, n_steps + 1)]
    grid[-1] = target
    s_prev, x_prev = start, x_start
    s_older = x_older = None
    consecutive = 0
    queue = list(grid)
    while queue:
        s = queue[0]
        if s == s_prev:
            return x_prev, True
        if s_older is None or s_prev == s_older:
            x_pred = x_prev
        else:
            factor = (s - s_prev) / (s_prev - s_older)
            x_pred = (
                x_prev[0] + (x_prev[0] - x_older[0]) * factor,
                x_prev[1] + (x_prev[1] - x_older[1]) * factor,
            )
        s2, s4 = slope * s + s, slope * s - s
        root = _newton(x_pred[0], x_pred[1], C, A, f_m1, s2, s4)
        guard = max(0.2 * max(abs(x_pred[0] - x_prev[0]), abs(x_pred[1] - x_prev[1])), 0.05)
        if (
            _accepted(root, f_m1, s2, s4)
            and max(abs(root[0] - x_pred[0]), abs(root[1] - x_pred[1])) <= guard
        ):
            s_older, x_older, s_prev, x_prev = s_prev, x_prev, s, (root[0], root[1])
            queue.pop(0)
            consecutive = 0
            continue
        consecutive += 1
        if consecutive > 40:
            return x_prev, True
        queue.insert(0, 0.5 * (s_prev + s))
    return x_prev, False


def _standard_starts(f_m1, s2, s4):
    """Starting points of the root census: five fixed and 25 amplitude-scaled pairs."""
    starts = [
        (50.0, -50.0),
        (0.2 * f_m1, -0.2 * f_m1),
        (0.02 * f_m1, -0.02 * f_m1),
        (0.5 * f_m1, -0.5 * f_m1),
        (1.0, -1.0),
    ]
    amp, mid = 0.5 * abs(s2 - s4), 0.5 * (s2 + s4)
    starts += [
        (f2 * amp + mid, f4 * amp)
        for f2 in (0.4, 0.2, 0.6, 0.1, 0.8)
        for f4 in (0.05, 0.0, 0.2, -0.1, 0.4)
    ]
    return starts


def _census(C, A, f_m1, s2, s4, n_grid=12):
    """Distinct admissible roots from a 12 x 12 grid over the admissible box plus 30 starts."""
    M = f_m1 / (C * A)
    span = s2 - s4
    a2_lo = max((s2 - M) / 1.5, s4 / 1.5 - span)
    a2_hi = s2 / 1.5
    if A < 1.0:
        a2_hi = min(a2_hi, (f_m1 + f_m1 * (1.0 - A) / (C * A)) / (1.5 * A))
    a4_lo = max(s4 / 1.5, (f_m1 - f_m1 * (A + 1.0) / (C * A)) / (1.5 * A))
    a4_hi = s2 / 1.5 + span
    e2 = 1e-6 * max(1.0, abs(a2_hi))
    e4 = 1e-6 * max(1.0, abs(a4_lo))
    starts = [
        (x2, x4)
        for x2 in np.linspace(a2_lo, a2_hi - e2, n_grid)
        for x4 in np.linspace(a4_lo + e4, a4_hi, n_grid)
    ]
    starts += _standard_starts(f_m1, s2, s4)
    roots: List[Tuple[float, float]] = []
    for x0 in starts:
        root = _newton(float(x0[0]), float(x0[1]), C, A, f_m1, s2, s4)
        if not _accepted(root, f_m1, s2, s4):
            continue
        if not any(
            abs(root[0] - r[0]) < 1e-5 * max(1.0, abs(root[0]))
            and abs(root[1] - r[1]) < 1e-5 * max(1.0, abs(root[1]))
            for r in roots
        ):
            roots.append((root[0], root[1]))
    return roots


def _same_root(x, roots) -> bool:
    return any(
        abs(x[0] - r[0]) < 1e-5 * max(1.0, abs(x[0]))
        and abs(x[1] - r[1]) < 1e-5 * max(1.0, abs(x[1]))
        for r in roots
    )


@dataclass(frozen=True)
class ClosedFormLife:
    """Closed-form life of the periodic uniaxial cycle (inf below the endurance limit, nan when
    the physical root is lost), its status and the back-stress at the two extremes."""

    Nf: float
    status: str
    alpha_2: float = math.nan
    alpha_4: float = math.nan


def _life_from_roots(C, K, L, A, f_m1, s2, s4, a2, a4) -> float:
    """Eqs. (28)-(29)."""
    t1 = (L / f_m1) * ((A + 1.0) * s2 - 1.5 * a2 - f_m1)
    t2 = (L / f_m1) * ((A - 1.0) * s4 + 1.5 * a4 - f_m1)
    if max(t1, t2) > 700.0:
        return 0.0
    brackets = math.exp(t1) + math.exp(t2) - 2.0
    return math.inf if brackets <= 0.0 else L / (K * brackets)


def uniaxial_life(
    C: float,
    K: float,
    L: float,
    sigma_m: float,
    sigma_a: float,
    A: float,
    f_m1: float,
    *,
    check_unique: bool = False,
) -> ClosedFormLife:
    """Closed-form life of the periodic uniaxial cycle sigma_m +- sigma_a (Eqs. 23, 24, 28, 29).

    Parameters
    ----------
    C, K, L : float
        Back-stress and damage constants of Eqs. (6)-(7).
    sigma_m, sigma_a : float
        Mean stress and amplitude of the axial cycle [MPa].
    A, f_m1 : float
        Slope of the endurance surface and fully reversed axial fatigue limit [MPa].
    check_unique : bool
        Also require that the physical root is the only admissible root, at the start of the
        branch and at the target (status "ambiguous_root" or "out_of_validity_no_start"
        otherwise). This is gate G1 of the calibration.

    Returns
    -------
    ClosedFormLife
        The life is a number only for the statuses "valid" and "below_endurance_limit".
    """
    if sigma_a <= 0.0:
        return ClosedFormLife(math.inf, BELOW)
    slope = sigma_m / sigma_a
    limit = endurance_amplitude(A, f_m1, sigma_m, sigma_a)
    if not math.isfinite(limit):
        return ClosedFormLife(math.nan, NO_START)
    if sigma_a <= limit:
        return ClosedFormLife(math.inf, BELOW)
    start = limit * (1.0 + START_OFFSET)
    if sigma_a <= start * (1.0 + 1e-9):
        start = sigma_a
    s2s, s4s = slope * start + start, slope * start - start
    guess = _limit_root(A, f_m1, slope * limit + limit, slope * limit - limit)
    root = _newton(guess[0], guess[1], C, A, f_m1, s2s, s4s)
    if not _accepted(root, f_m1, s2s, s4s):
        return ClosedFormLife(math.nan, NO_START)
    x = (root[0], root[1])
    if check_unique:
        others = [r for r in _census(C, A, f_m1, s2s, s4s) if not _same_root(r, [x])]
        if others:
            return ClosedFormLife(math.nan, NO_START)
    if start < sigma_a:
        x, folded = _follow_branch(C, A, f_m1, slope, start, x, sigma_a, CONTINUATION_STEPS)
        if folded:
            return ClosedFormLife(math.nan, FOLD)
    s2, s4 = sigma_m + sigma_a, sigma_m - sigma_a
    Nf = _life_from_roots(C, K, L, A, f_m1, s2, s4, x[0], x[1])
    if check_unique:
        others = [r for r in _census(C, A, f_m1, s2, s4) if not _same_root(r, [x])]
        if others:
            return ClosedFormLife(Nf, AMBIGUOUS, x[0], x[1])
    return ClosedFormLife(Nf, VALID, x[0], x[1])


def eq25_objective(
    C: float,
    K: float,
    L: float,
    points: Sequence[Tuple[float, float, float]],
    A: float,
    f_m1: float,
) -> float:
    """Tveit et al. (2025) Eq. (25) over points (sigma_m, sigma_a, N_exp); 10 per lost prediction."""
    total = 0.0
    for sigma_m, sigma_a, N_exp in points:
        Nf = uniaxial_life(C, K, L, sigma_m, sigma_a, A, f_m1).Nf
        if not math.isfinite(Nf) or Nf <= 1.0:
            total += PENALTY
        else:
            total += (1.0 - math.log(Nf) / math.log(N_exp)) ** 2
    return total


# ---------------------------------------------------------------------------------------------
# Simplex minimisation (Nelder & Mead 1965, as used by Ottosen et al. 2008)
# ---------------------------------------------------------------------------------------------


def _nelder_mead(fun, x0, lower, upper, max_evaluations=MAX_EVALUATIONS):
    """Bounded simplex minimisation; points outside the box are projected onto it.

    Returns (x, f(x), evaluations, converged). Converged: simplex within 1e-8 in every
    coordinate and objective values within 1e-12.
    """
    n = len(x0)
    span = upper - lower
    simplex = [np.clip(np.asarray(x0, float), lower, upper)]
    for i in range(n):
        vertex = simplex[0].copy()
        step = 0.05 * span[i]
        vertex[i] = vertex[i] + step if vertex[i] + step <= upper[i] else vertex[i] - step
        simplex.append(vertex)
    values = [fun(v) for v in simplex]
    evaluations = n + 1
    converged = False
    while evaluations < max_evaluations:
        order = np.argsort(values, kind="stable")
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        if (
            max(float(np.max(np.abs(v - simplex[0]))) for v in simplex[1:]) <= 1e-8
            and values[-1] - values[0] <= 1e-12
        ):
            converged = True
            break
        centroid = np.mean(simplex[:-1], axis=0)
        reflected = np.clip(2.0 * centroid - simplex[-1], lower, upper)
        f_r = fun(reflected)
        evaluations += 1
        if f_r < values[0]:
            expanded = np.clip(3.0 * centroid - 2.0 * simplex[-1], lower, upper)
            f_e = fun(expanded)
            evaluations += 1
            simplex[-1], values[-1] = (expanded, f_e) if f_e < f_r else (reflected, f_r)
            continue
        if f_r < values[-2]:
            simplex[-1], values[-1] = reflected, f_r
            continue
        if f_r < values[-1]:
            contracted = np.clip(centroid + 0.5 * (reflected - centroid), lower, upper)
            f_c = fun(contracted)
            evaluations += 1
            if f_c <= f_r:
                simplex[-1], values[-1] = contracted, f_c
                continue
        else:
            contracted = np.clip(centroid + 0.5 * (simplex[-1] - centroid), lower, upper)
            f_c = fun(contracted)
            evaluations += 1
            if f_c < values[-1]:
                simplex[-1], values[-1] = contracted, f_c
                continue
        for i in range(1, n + 1):
            simplex[i] = simplex[0] + 0.5 * (simplex[i] - simplex[0])
            values[i] = fun(simplex[i])
            evaluations += 1
    best = int(np.argmin(values))
    return simplex[best], values[best], evaluations, converged


def c_grid() -> List[float]:
    """Fixed C values of the profile: 8 per decade over the C bounds, both bounds included."""
    low, high = (math.log10(b) for b in BOUNDS["C"])
    n = int(round((high - low) * C_PER_DECADE))
    return [10.0 ** (low + (high - low) * j / n) for j in range(n + 1)]


@dataclass
class FitResult:
    """The profile of Eq. (25) over C, its minimum and the equivalent family (in increasing E)."""

    profile: List[Dict[str, Any]]
    E_min: float
    family: List[Dict[str, Any]] = field(default_factory=list)


def fit_ckl(points: Sequence[Tuple[float, float, float]], A: float, f_m1: float) -> FitResult:
    """Eq. (25) profiled over C: at each C of c_grid(), the best K and L (Repository choices).

    K and L are fitted by the simplex method in ln K, ln L, sweeping the grid upwards and
    downwards with the neighbouring optimum as start (plus a fixed start, K = 1e-5, L = 5); the
    better of the two sweeps is kept. The equivalent family is the grid points with
    E <= 1.05 E_min (when E_min < 1e-6: those within 0.01 decade at every point), ordered by E.
    """
    lower = np.log([BOUNDS["K"][0], BOUNDS["L"][0]])
    upper = np.log([BOUNDS["K"][1], BOUNDS["L"][1]])
    fixed = np.log([1e-5, 5.0])
    grid = c_grid()
    best: Dict[int, Dict[str, Any]] = {}
    for order in (range(len(grid)), range(len(grid) - 1, -1, -1)):
        previous = None
        for j in order:
            C = grid[j]

            def objective(x, C=C) -> float:
                K, L = np.exp(x)
                return eq25_objective(C, float(K), float(L), points, A, f_m1)

            tried = [fixed] if previous is None else [previous, fixed]
            for x0 in tried:
                x, value, _, _ = _nelder_mead(objective, x0, lower, upper)
                if j not in best or value < best[j]["E"]:
                    K, L = (float(v) for v in np.exp(x))
                    best[j] = {"C": C, "K": K, "L": L, "E": float(value)}
            previous = np.log([best[j]["K"], best[j]["L"]])
    profile = [best[j] for j in range(len(grid))]
    finite = [p for p in profile if math.isfinite(p["E"]) and p["E"] < PENALTY]
    if not finite:
        return FitResult(profile, math.nan)
    E_min = min(p["E"] for p in finite)
    if E_min >= EXACT_FIT:
        family = [p for p in finite if p["E"] <= EQUIVALENT * E_min]
    else:
        family = [
            p
            for p in finite
            if all(
                abs(math.log10(uniaxial_life(p["C"], p["K"], p["L"], m, a, A, f_m1).Nf / n)) <= 0.01
                for m, a, n in points
            )
        ]
    return FitResult(profile, E_min, sorted(family, key=lambda p: (p["E"], p["C"])))


# ---------------------------------------------------------------------------------------------
# Admissibility (G1, G2) and choice
# ---------------------------------------------------------------------------------------------


def numerical_settings(inputs: Dict[str, Any]) -> Dict[str, float]:
    """Time-integration settings of the model (constants-file names): inputs file, else defaults."""
    from multifat.models.invariant_methods.tveit2025 import SETTINGS

    declared = inputs.get("numerical_settings") or {}
    unknown = sorted(set(declared) - set(SETTINGS))
    if unknown:
        raise CalibrationInputError(f"unknown numerical_settings: {', '.join(unknown)}")
    return {key: declared.get(key, default) for key, (_, default) in SETTINGS.items()}


def _transient_life(C, K, L, sigma_m, sigma_a, A, f_m1, t_m1, settings) -> float:
    """Time-integrated life of the model on the axial cycle (axial slot sigma_yy)."""
    from multifat.models.invariant_methods.tveit2025 import SETTINGS, tveit2025_life

    S_mean = np.zeros((3, 3))
    S_sin = np.zeros((3, 3))
    S_mean[1, 1] = sigma_m
    S_sin[1, 1] = sigma_a
    keywords = {SETTINGS[key][0]: value for key, value in settings.items()}
    return tveit2025_life(S_mean, S_sin, np.zeros((3, 3)), f_m1, t_m1, A, C, K, L, **keywords).Nf


def check_admissible(
    candidate: Dict[str, Any],
    points: Sequence[Tuple[float, float, float]],
    A: float,
    f_m1: float,
    t_m1: float,
    settings: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """Gates G1 (unique closed-form root) and G2 (transient within 0.10 decade) at every point."""
    settings = numerical_settings({}) if settings is None else settings
    rows = []
    G1 = G2 = True
    for sigma_m, sigma_a, N_exp in points:
        closed = uniaxial_life(
            candidate["C"],
            candidate["K"],
            candidate["L"],
            sigma_m,
            sigma_a,
            A,
            f_m1,
            check_unique=True,
        )
        G1 = G1 and closed.status == VALID
        transient = _transient_life(
            candidate["C"],
            candidate["K"],
            candidate["L"],
            sigma_m,
            sigma_a,
            A,
            f_m1,
            t_m1,
            settings,
        )
        if math.isfinite(transient) and math.isfinite(closed.Nf) and closed.Nf > 0.0:
            gap = math.log10(transient) - math.log10(closed.Nf)
        else:
            gap = math.nan
        G2 = G2 and math.isfinite(gap) and abs(gap) <= G2_DECADES
        rows.append(
            {
                "sigma_m": sigma_m,
                "sigma_a": sigma_a,
                "N_exp": N_exp,
                "N_closed": closed.Nf,
                "status": closed.status,
                "N_transient": transient,
                "log10_transient_minus_closed": gap,
                "residual_closed": _log10_ratio(closed.Nf, N_exp),
                "residual_transient": _log10_ratio(transient, N_exp),
            }
        )
    return {"G1": bool(G1), "G2": bool(G2), "points": rows}


def _log10_ratio(N, N_exp) -> float:
    return math.log10(N / N_exp) if math.isfinite(N) and N > 0.0 else math.nan


# ---------------------------------------------------------------------------------------------
# Inputs, A and the calibration tests
# ---------------------------------------------------------------------------------------------


def read_inputs(path: Optional[Path], material: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """The optional hand-written inputs and their SHA-256 ({} and None when there is no file)."""
    if path is None or not Path(path).is_file():
        return {}, None
    path = Path(path)
    raw = path.read_bytes()
    inputs = json.loads(raw.decode("utf-8"))
    if inputs.get("schema") != INPUTS_SCHEMA:
        raise CalibrationInputError(f"{path.as_posix()}: schema must be {INPUTS_SCHEMA!r}")
    if inputs.get("material_id") != material:
        raise CalibrationInputError(
            f"{path.as_posix()} holds inputs for '{inputs.get('material_id')}', not '{material}'"
        )
    unknown = sorted(set(inputs) - INPUTS_KEYS)
    if unknown:
        raise CalibrationInputError(f"{path.as_posix()}: unknown keys {', '.join(unknown)}")
    # Hash with LF line ends, so that a Windows checkout (CRLF) gives the same digest
    return inputs, hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest()


def fatigue_limits(block: Dict[str, Any]) -> Dict[str, Any]:
    """f_m1 and t_m1 of the workbook Summary (R = -1 rows) with the rows they come from."""
    properties = block["material"]
    f_m1 = properties.get("Tensile_fatigue_limit_MPa")
    t_m1 = properties.get("Torsional_fatigue_limit_MPa")
    issues = block.get("property_issues") or {}
    for key in ("Sigm1", "Taum1"):
        if key in issues:
            raise CalibrationInputError(f"{key} of the Summary is not usable: {issues[key]}")
    if f_m1 is None or t_m1 is None:
        raise CalibrationInputError("TVEIT2025 needs Sigm1 and Taum1 in the workbook Summary")
    if not (float(f_m1) > 0.0 and float(t_m1) > 0.0):
        raise CalibrationInputError(f"Sigm1 and Taum1 must be positive, got {f_m1}, {t_m1}")
    rows = (block.get("fatigue_limit_sources") or {}).get("Sigm1") or []
    return {
        "Sigm1": float(f_m1),
        "Taum1": float(t_m1),
        "origin": "workbook Summary",
        "sources": [
            {k: row.get(k) for k in ("authors", "reference", "reference_cycles", "excel_row")}
            for row in rows
        ],
    }


def mean_stress_slope(
    inputs: Dict[str, Any], f_m1: float, properties: Dict[str, Any]
) -> Dict[str, Any]:
    """A in the order of Tveit et al. (2025, Sect. 4): mean-stress limit, given A, FKM, Goodman."""
    limit, given = inputs.get("mean_stress_fatigue_limit"), inputs.get("A")
    if limit and given:
        raise CalibrationInputError("declare either mean_stress_fatigue_limit or A, not both")
    material_type = inputs.get("material_type")
    note = ""
    if limit:
        sigma_a0 = float(limit["sigma_a"])
        if "sigma_m" in limit:
            sigma_m0 = float(limit["sigma_m"])
        elif "R" in limit:
            R = float(limit["R"])
            if not -1.0 < R < 1.0:
                raise CalibrationInputError(f"mean-stress fatigue limit needs -1 < R < 1, got {R}")
            sigma_m0 = sigma_a0 * (1.0 + R) / (1.0 - R)
        else:
            raise CalibrationInputError("mean_stress_fatigue_limit needs sigma_m or R")
        if sigma_m0 <= 0.0:
            raise CalibrationInputError("the mean-stress fatigue limit needs a tensile mean stress")
        A = (f_m1 - sigma_a0) / sigma_m0
        route, source = "mean_stress_fatigue_limit", limit.get("source", "")
        how = (
            f"fatigue limit sigma_a = {sigma_a0:g} MPa at sigma_m = {sigma_m0:.6g} MPa on the "
            "Haigh line: A = (Sigm1 - sigma_a) / sigma_m (Tveit et al. 2025, Eqs. 22-23)"
        )
    elif given:
        A, route, source = float(given["value"]), "supplied", given.get("source", "")
        how = "hand-entered value"
    else:
        Su = properties.get("Tensile_strength_MPa")
        if Su is None or not float(Su) > 0.0:
            raise CalibrationInputError(
                "A needs a fatigue limit at a tensile mean stress or a hand-entered A (inputs "
                "file), or the tensile strength Su (Summary) for an estimate"
            )
        Su = float(Su)
        if material_type in FKM_MEAN_STRESS:
            a_M, b_M = FKM_MEAN_STRESS[material_type]
            A, route = a_M * Su / 1000.0 + b_M, "fkm"
            how = (
                f"FKM estimate for {material_type}: A = a_M Su/1000 + b_M with a_M = {a_M:g}, "
                f"b_M = {b_M:g}, Su = {Su:g} MPa (Tveit et al. 2025, Eq. 26)"
            )
            source = "FKM guideline (Haenel et al. 2003); Su from the workbook Summary"
        else:
            A, route = f_m1 / Su, "goodman"
            how = (
                f"Goodman estimate A = Sigm1 / Su with Su = {Su:g} MPa (Tveit et al. 2025, Sect. 4)"
            )
            source = "Goodman relation; Sigm1 and Su from the workbook Summary"
            # As the authors do where the FKM guideline gives no values (their 7:3 brass)
            note = (
                "material type is not declared"
                if material_type is None
                else f"FKM_MEAN_STRESS has no values for material type {material_type!r}"
            )
    if not (0.0 < A < 1.0):
        raise CalibrationInputError(f"A = {A:.6g} is outside (0, 1), where Eqs. (23)-(24) hold")
    return {
        "A": float(A),
        "route": route,
        "equation": how,
        "source": source,
        "material_type": material_type,
        "goodman_because": note,
    }


def _normal(text) -> str:
    out = str(text or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if out.startswith(prefix):
            out = out[len(prefix) :]
    return out


def dataset_sources(workbook: Path) -> Dict[int, Tuple[str, str]]:
    """(authors, reference) of every Dataset row by Excel row, with the overlay applied."""
    import openpyxl

    from multifat.imp import experimental_data as exp
    from multifat.imp.overlays import apply_overlay_to_workbook

    workbook_file = openpyxl.load_workbook(workbook, data_only=True)
    try:
        if exp.OVERLAYS_ENABLED:
            with warnings.catch_warnings():  # the importer has already reported a mismatch
                warnings.simplefilter("ignore")
                apply_overlay_to_workbook(workbook_file, workbook)
        sheet = workbook_file["Dataset"]
        rows = sheet.iter_rows(values_only=True)
        header = [str(h) if h is not None else "" for h in next(rows)]
        authors = header.index("Authors") if "Authors" in header else None
        refs = [i for i, h in enumerate(header) if h.startswith("Article DOI")]
        out = {}
        for excel_row, row in enumerate(rows, start=2):
            name = row[authors] if authors is not None and authors < len(row) else None
            ref = row[refs[0]] if refs and refs[0] < len(row) else None
            out[excel_row] = (_normal(name), _normal(ref))
        return out
    finally:
        workbook_file.close()


def axial_tests(cases: Dict[int, Any], sources: Dict[int, Tuple[str, str]], failures: bool):
    """Axial R = -1 tests (shared selector) with sigma_m, sigma_a, life and source."""
    from multifat.imp.experimental_data import is_design_test

    out = []
    for case_id, case in cases.items():
        if not is_design_test(case, "tension"):
            continue
        meta = case.get("meta", {})
        if failures and meta.get("test_status", "failure") != "failure":
            continue
        axial = [float(np.asarray(case[1][k]["stress"], dtype=float)[1, 1]) for k in (1, 2)]
        row = int(meta.get("excel_row", -1))
        out.append(
            {
                "case_id": int(case_id),
                "specimen_id": meta.get("specimen_id"),
                "excel_row": row,
                "sigma_m": 0.5 * (max(axial) + min(axial)),
                "sigma_a": 0.5 * (max(axial) - min(axial)),
                "N_exp": float(meta["Nf_cycles"]),
                "source": list(sources.get(row, ("", ""))),
            }
        )
    return out


def calibration_source(
    tests: List[Dict[str, Any]], limits: Dict[str, Any], declared: Optional[str]
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """The tests of one source: the declared one, else that of the Summary fatigue limits.

    Returns the tests and how they were matched: "declared", "same reference", "same authors"
    (another publication of the authors of the limits), "no match" or "unknown" (the Summary names
    no source); in the last two cases every test is returned.
    """
    if declared:
        key = _normal(declared)
        chosen = [t for t in tests if key in t["source"][0] or key in t["source"][1]]
        return chosen, {"match": "declared", "source": declared}
    refs = {_normal(s.get("reference")) for s in limits["sources"]} - {""}
    names = {_normal(s.get("authors")) for s in limits["sources"]} - {""}
    if not refs and not names:
        return tests, {"match": "unknown", "source": ""}
    same_ref = [t for t in tests if t["source"][1] in refs]
    if same_ref:
        return same_ref, {"match": "same reference", "source": ", ".join(sorted(refs))}
    same_authors = [t for t in tests if t["source"][0] in names]
    if same_authors:
        return same_authors, {"match": "same authors", "source": ", ".join(sorted(names))}
    return tests, {"match": "no match", "source": ", ".join(sorted(refs | names))}


def three_point_rule(candidates: List[Dict[str, Any]], f_m1: float) -> List[Dict[str, Any]]:
    """LOW, MID and HIGH of the axial R = -1 failures above f_m1 (Repository choices)."""
    pool = [c for c in candidates if c["sigma_a"] > f_m1]
    levels = sorted({c["sigma_a"] for c in pool})
    if len(levels) < 3:
        raise CalibrationInputError(
            f"the procedure needs axial R = -1 failures at three stress levels above "
            f"Sigm1 = {f_m1:g} MPa; found {len(levels)} level(s)"
        )
    low, high = levels[0], levels[-1]
    centre = 0.5 * (math.log(low) + math.log(high))
    middle = sorted(levels[1:-1], key=lambda s: (abs(math.log(s) - centre), s))[0]
    chosen = []
    for role, level in (("LOW", low), ("MID", middle), ("HIGH", high)):
        replicates = [c for c in pool if c["sigma_a"] == level]
        median = float(np.median([math.log10(c["N_exp"]) for c in replicates]))
        pick = sorted(
            replicates, key=lambda c: (abs(math.log10(c["N_exp"]) - median), c["excel_row"])
        )[0]
        chosen.append(
            {
                **pick,
                "role": role,
                "replicates": len(replicates),
                "reason": (
                    f"{role} level {level:.6g} MPa ({len(replicates)} replicate(s)); log10 N "
                    f"closest to the level median {median:.4f}"
                ),
            }
        )
    return chosen


# ---------------------------------------------------------------------------------------------
# Confidence of the calibration
# ---------------------------------------------------------------------------------------------


def _input_findings(f_m1, failures, runouts, matched, slope) -> List[Finding]:
    """Findings on the limits, the source of the tests and A (Steps 1-2 of the method guide)."""
    found = []
    below = [t for t in failures if t["sigma_a"] < f_m1]
    if below:
        found.append(
            Finding(
                "failures_below_limit",
                MAJOR,
                f"{len(below)} axial R = -1 failure(s) lie below Sigm1 = {f_m1:g} MPa "
                f"(lowest {min(t['sigma_a'] for t in below):.4g} MPa)",
                "the model predicts an infinite life for every test below its fatigue limit",
                "check Sigm1 against the failures and run-outs of the same source, e.g. with a "
                "censored fit of the limit",
            )
        )
    above = [t for t in runouts if t["sigma_a"] > f_m1]
    if above:
        found.append(
            Finding(
                "runouts_above_limit",
                MINOR,
                f"{len(above)} axial R = -1 run-out(s) lie above Sigm1 = {f_m1:g} MPa "
                f"(highest {max(t['sigma_a'] for t in above):.4g} MPa)",
                "the fatigue limit may be conservative",
                "check Sigm1 against the run-outs of the same source",
            )
        )
    if matched["match"] == "no match":
        found.append(
            Finding(
                "limit_source_differs",
                MAJOR,
                f"no axial test comes from the source of the fatigue limits ({matched['source']})",
                "a limit applied to another source's tests puts tests on the wrong side of the "
                "endurance surface",
                "use the fatigue limits of the source of the calibration tests",
            )
        )
    elif matched["match"] == "same authors":
        found.append(
            Finding(
                "limit_other_publication",
                MINOR,
                "the fatigue limits come from another publication of the authors of the tests "
                f"({matched['source']})",
                "the limits and the tests may belong to different specimen batches",
                "confirm that the limits and the tests are of the same batch",
            )
        )
    if slope["route"] in ("fkm", "goodman"):
        found.append(
            Finding(  # minor: the authors' own route when no mean-stress limit is available
                "A_estimated",
                MINOR,
                f"A = {slope['A']:.4g} is estimated ({slope['route']}), not measured",
                "every mean-stress prediction depends on A",
                "a fatigue limit at a tensile mean stress (e.g. an R = 0.1 staircase) gives A "
                "directly (Tveit et al. 2025, Eq. 23)",
            )
        )
    if slope["route"] == "goodman" and slope["material_type"] is None:
        found.append(
            Finding(
                "material_type_for_FKM",
                MINOR,
                f"Goodman was used because the {slope['goodman_because']}",
                "the authors estimate A with the FKM guideline where it gives values for the "
                "material type",
                "declare material_type in the inputs file (FKM values exist for: "
                + ", ".join(sorted(FKM_MEAN_STRESS))
                + ")",
            )
        )
    if slope["route"] == "supplied":
        found.append(
            Finding(
                "A_supplied",
                MINOR,
                f"A = {slope['A']:.4g} is entered by hand ({slope['source'] or 'no source given'})",
                "the calibration cannot check how A was obtained",
                "a fatigue limit at a tensile mean stress gives A directly",
            )
        )
    return found


def _point_findings(points, f_m1) -> List[Finding]:
    """Findings on the three calibration points (Step 3 of the method guide)."""
    found = []
    off_axis = [p for p in points if abs(p["sigma_m"]) > 1e-3 * p["sigma_a"]]
    if off_axis:
        found.append(
            Finding(
                "point_not_R_minus_1",
                MAJOR,
                f"{len(off_axis)} calibration point(s) have a mean stress",
                "the closed form is fitted on fully reversed points; A then mixes into the fit",
                "use axial R = -1 failures",
            )
        )
    levels = sorted(p["sigma_a"] for p in points)
    if levels[-1] / levels[0] < NARROW_SPAN:
        found.append(
            Finding(
                "narrow_stress_range",
                MAJOR,
                f"the points span only {levels[0]:.4g} to {levels[-1]:.4g} MPa "
                f"({100 * (levels[-1] / levels[0] - 1):.1f} %)",
                "points this close carry little S-N trend, and K and L are poorly determined",
                "add failures over a wider stress range of the same source",
            )
        )
    if levels[0] / f_m1 < NEAR_LIMIT:
        found.append(
            Finding(
                "low_point_near_limit",
                MINOR,
                f"the LOW point ({levels[0]:.4g} MPa) lies within "
                f"{100 * (levels[0] / f_m1 - 1):.1f} % of Sigm1 = {f_m1:g} MPa",
                "near its asymptote the model curve is steep; a residual of a few tenths of a "
                "decade is expected there",
                "a LOW level further above the fatigue limit",
            )
        )
    single = [p for p in points if p["replicates"] == 1]
    if single:
        found.append(
            Finding(
                "single_specimen_level",
                MINOR,
                "single specimen at the "
                + ", ".join(f"{p['role']} level ({p['sigma_a']:.4g} MPa)" for p in single),
                "one specimen that disagrees with other data can bend the whole curve",
                "more replicates at that level",
            )
        )
    sources = {tuple(p["source"]) for p in points}
    if len(sources) > 1:
        found.append(
            Finding(
                "points_mixed_sources",
                MAJOR,
                f"the calibration points come from {len(sources)} sources",
                "the scatter between sources enters the constants",
                "calibrate on the tests of one source",
            )
        )
    return found


def _fit_findings(chosen, family) -> List[Finding]:
    """Findings on the fit and the chosen set (Steps 4-5 of the method guide)."""
    found = []
    residuals = [
        abs(row[key])
        for row in chosen["points"]
        for key in ("residual_closed", "residual_transient")
        if math.isfinite(row[key])
    ]
    worst = max(residuals, default=0.0)
    if worst > POOR_FIT_DECADES:
        found.append(
            Finding(
                "poor_fit",
                MINOR,
                f"a calibration point is missed by {worst:.2f} decade",
                "with Sigm1 fixed the curve cannot pass through three arbitrary points; a large "
                "residual points at the limit or at one of the points",
                "check the fatigue limit and the chosen points",
            )
        )
    # The fit works in ln K and ln L, so a value on a bound comes back within rounding of it.
    # C on its bound counts only when the equivalent sets pin C: when they span decades, E is flat
    # in C and the lowest E falls on a grid end by rounding (finding C_not_identified).
    c_low, c_high = family["C_range"]
    c_free = math.log10(c_high / c_low) > C_SPAN_DECADES
    bounds = [
        name
        for name in ("C", "K", "L")
        if any(abs(math.log(chosen[name] / b)) < 1e-9 for b in BOUNDS[name])
        and not (name == "C" and c_free)
    ]
    if bounds:
        found.append(
            Finding(
                "constant_at_bound",
                MAJOR,
                f"{', '.join(bounds)} of the chosen set lies on its search bound",
                "a constant on its bound is set by the search, not by the data",
                "more informative tests (another stress level, a mean stress, a multiaxial "
                "test); the bound is not a calibrated value",
            )
        )
    if c_free:
        found.append(
            Finding(
                "C_not_identified",
                MINOR,
                f"C varies from {c_low:.3g} to {c_high:.3g} among the equivalent sets",
                "uniaxial data cannot fix C; it matters under mean shear stress and "
                "non-proportional loading",
                "a test with a mean shear stress or 90 deg out-of-phase tension-torsion decides C",
            )
        )
    return found


# ---------------------------------------------------------------------------------------------
# Reuse of an identical calibration
# ---------------------------------------------------------------------------------------------


def calibration_identity(
    workbook: Path, block: Dict[str, Any], inputs_sha256: Optional[str]
) -> Dict[str, Any]:
    """What determines the calibration: workbook bytes, overlay, inputs, procedure and revision."""
    from multifat.models.references import IMPLEMENTATION_REVISION

    return {
        "procedure": PROCEDURE,
        "implementation_revision": IMPLEMENTATION_REVISION,
        "workbook_sha256": hashlib.sha256(Path(workbook).read_bytes()).hexdigest(),
        "overlay_sha256": (block.get("overlay") or {}).get("sha256"),
        "inputs_sha256": inputs_sha256,
    }


def reusable_calibration(
    folders: Sequence[Path], identity: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """The TVEIT2025 calibration of the earliest archived record with this identity, if any.

    The calibration is deterministic, so a record with the same identity holds the result a new
    calibration would give (docs/calibration-validation-protocol.md: identical deterministic
    reruns may reuse the same record).
    """
    for folder in folders:
        for path in sorted(Path(folder).glob("*.json")):
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if not isinstance(record, dict):
                continue
            entry = (record.get("model_calibrations") or {}).get("TVEIT2025") or {}
            if entry.get("status") == "calibrated" and entry.get("identity") == identity:
                return {
                    **entry,
                    "reused_from": {
                        "calibration_id": record.get("calibration_id", path.stem),
                        "file": path.as_posix(),
                    },
                }
    return None


# ---------------------------------------------------------------------------------------------
# Calibration of one workbook
# ---------------------------------------------------------------------------------------------


def calibrate_tveit2025(
    workbook: Path,
    data: Dict[str, Any],
    material: str,
    *,
    inputs: Optional[Path] = None,
    progress: Optional[Callable[[str], None]] = None,
    archives: Sequence[Path] = (),
) -> Dict[str, Any]:
    """Calibrate A, C, K and L of TVEIT2025 on one workbook by the authors' procedure; JSON-ready.

    Parameters
    ----------
    workbook : Path
        The workbook (its Dataset sheet names the source of each test).
    data, material : dict, str
        The imported workbook and its material ID.
    inputs : Path, optional
        The inputs file; default Calibrations/<material>/TVEIT2025_inputs.json when it exists.
    progress : callable, optional
        Receives one line of text at each step.
    archives : sequence of Path
        Further record folders of this material searched, after Calibrations/<material>/, for a
        calibration with the same identity, which is then reused instead of repeated.

    Returns
    -------
    dict
        "status" ("calibrated", "insufficient_calibration_data" or "infeasible"), "reason",
        "params" (Sigm1_TVEIT, Taum1_TVEIT, A_TVEIT, C_TVEIT, K_TVEIT, L_TVEIT and the
        time-integration settings), "source", "design_case_ids", "confidence", "identity", the
        full "record", and "reused_from" when an archived calibration was reused.
    """
    say = progress or (lambda text: None)
    path = Path(inputs) if inputs is not None else CALIBRATIONS_DIR / material / INPUTS_FILE
    record: Dict[str, Any] = {
        "procedure": PROCEDURE,
        "mode": MODE,
        "references": [
            "Tveit et al. (2025) Mech. Mater. 209, 105452, Sect. 4, Eqs. (22)-(26)",
            "Ottosen et al. (2008) Int. J. Fatigue 30, 996-1006, Sect. 5, Eqs. (23)-(32)",
        ],
    }
    out: Dict[str, Any] = {
        "params": {},
        "source": {},
        "design_case_ids": {},
        "confidence": None,
        "record": record,
    }
    try:
        given, digest = read_inputs(path, material)
        record.update(
            inputs_file=path.as_posix() if digest else None, inputs_sha256=digest, inputs=given
        )
        block = data[material]
        out["identity"] = calibration_identity(Path(workbook), block, digest)
        limits = fatigue_limits(block)
        f_m1, t_m1 = limits["Sigm1"], limits["Taum1"]
        slope = mean_stress_slope(given, f_m1, block["material"])
        settings = numerical_settings(given)
        sources = dataset_sources(Path(workbook))
        failures = axial_tests(block["cases"], sources, failures=True)
        runouts = axial_tests(block.get("runouts", {}), sources, failures=False)
        pool, matched = calibration_source(failures, limits, given.get("calibration_source"))
        record.update(
            fatigue_limits=limits,
            mean_stress_slope=slope,
            numerical_settings=settings,
            calibration_source=matched,
            candidates=len(pool),
        )
        findings = _input_findings(f_m1, failures, runouts, matched, slope)
        points = three_point_rule(pool, f_m1)
    except CalibrationInputError as exc:
        out.update(status="insufficient_calibration_data", reason=str(exc))
        return out
    findings += _point_findings(points, f_m1)
    earlier = reusable_calibration([CALIBRATIONS_DIR / material, *archives], out["identity"])
    if earlier is not None:
        # Same constants; the confidence is judged again with the current rules
        chosen, family = earlier["record"]["chosen"], earlier["record"]["family"]
        confidence = summary(findings + _fit_findings(chosen, family))
        earlier["confidence"] = earlier["record"]["confidence"] = confidence
        say(
            f"TVEIT2025 {material}: reusing calibration {earlier['reused_from']['calibration_id']}"
            " (same workbook, overlay, inputs, procedure and implementation revision)"
        )
        for line in report_lines(confidence):
            say(f"TVEIT2025 {material}: {line}")
        return earlier
    A = slope["A"]
    record["calibration_points"] = points
    ids = [p["case_id"] for p in points]
    record["validation"] = (
        "every other test of the workbook is an independent assessment of these constants; the "
        "calibration tests are design points (Flag = 0)"
    )
    triples = [(p["sigma_m"], p["sigma_a"], p["N_exp"]) for p in points]
    record["objective"] = "Eq. (25): sum (1 - ln N / ln N_exp)^2, closed-form N, 10 per lost point"
    record["optimiser"] = {
        "method": "profile over C: simplex in ln K, ln L at each C of a fixed grid",
        "C_grid": f"{C_PER_DECADE} per decade over the C bounds",
        "bounds": BOUNDS,
        "max_evaluations_per_simplex": MAX_EVALUATIONS,
        "family": f"E <= {EQUIVALENT:g} E_min",
        "G2_decades": G2_DECADES,
    }
    say(f"TVEIT2025 {material}: Eq. (25) profile over C on {len(points)} points")
    fit = fit_ckl(triples, A, f_m1)
    record["E_min"] = fit.E_min
    record["profile"] = fit.profile
    if not fit.family:
        out.update(status="infeasible", reason="no C of the grid gave a finite Eq. (25) objective")
        return out
    family = fit.family
    record["family"] = {
        "n_equivalent": len(family),
        "C_range": [min(f["C"] for f in family), max(f["C"] for f in family)],
        "K_range": [min(f["K"] for f in family), max(f["K"] for f in family)],
        "L_range": [min(f["L"] for f in family), max(f["L"] for f in family)],
    }
    checked = []
    chosen = None
    for member in family:  # in increasing E: the first admissible set is the chosen one
        gates = check_admissible(member, triples, A, f_m1, t_m1, settings)
        checked.append({"C": member["C"], "E": member["E"], **gates})
        say(f"TVEIT2025 {material}: C {member['C']:.4g} G1 {gates['G1']} G2 {gates['G2']}")
        if gates["G1"] and gates["G2"]:
            chosen = {**member, **gates}
            break
    record["admissibility_checked"] = checked
    if chosen is None:
        out.update(status="infeasible", reason="no equivalent (C, K, L) set passed gates G1 and G2")
        return out
    grid = c_grid()
    chosen["C_at_search_bound"] = chosen["C"] in (grid[0], grid[-1])
    record["chosen"] = chosen
    c_low, c_high = record["family"]["C_range"]
    record["identifiability"] = (
        f"C is not identified by uniaxial data: equivalent sets (E <= {EQUIVALENT:g} E_min) span "
        f"C {c_low:.3g} to {c_high:.3g}; the chosen C is the admissible one of lowest E"
        + (" and lies on the search bound" if chosen["C_at_search_bound"] else "")
    )
    findings += _fit_findings(chosen, record["family"])
    confidence = summary(findings)
    record["confidence"] = out["confidence"] = confidence
    for line in report_lines(confidence):
        say(f"TVEIT2025 {material}: {line}")
    out["params"] = {
        "Sigm1_TVEIT": f_m1,
        "Taum1_TVEIT": t_m1,
        "A_TVEIT": A,
        "C_TVEIT": chosen["C"],
        "K_TVEIT": chosen["K"],
        "L_TVEIT": chosen["L"],
        **settings,
    }
    fitted = f"TVEIT2025 Eq. (25), {MODE} procedure, on {len(points)} axial R = -1 failures"
    out["source"] = {
        "Sigm1_TVEIT": "fatigue limit (workbook Summary)",
        "Taum1_TVEIT": "fatigue limit (workbook Summary)",
        "A_TVEIT": slope["equation"],
        "C_TVEIT": fitted + "; " + record["identifiability"],
        "K_TVEIT": fitted,
        "L_TVEIT": fitted,
        **{key: "time-integration setting checked by gate G2" for key in settings},
    }
    out["design_case_ids"] = {name: ids for name in ("C_TVEIT", "K_TVEIT", "L_TVEIT")}
    out.update(status="calibrated", reason="")
    return out
