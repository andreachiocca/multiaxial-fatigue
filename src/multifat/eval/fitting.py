"""Power-law calibration using a damage-space shift that always stays well defined.

Repository choices:
  - the scatter of the master curve is estimated with n - 2 degrees of freedom, as for any
    two-parameter regression (ASTM E739; ``material_fits.fit_basquin`` does the same), so that
    STD_ln_DP**2 is unbiased; with n = 2 it is undefined (issue #74);
  - a master curve is identifiable only with at least ``MIN_POINTS`` tests over ``MIN_DECADES``
    decades of life and a negative slope distinct from zero (|b| > ``MIN_T`` standard errors),
    the thresholds of the Basquin fits of ``material_fits``. A model whose curve is not
    identifiable is skipped (maintainer, issue #74).
"""

from __future__ import annotations

import numpy as np

# Identifiability thresholds of a fitted S-N or master curve (shared with material_fits)
MIN_POINTS = 5
MIN_DECADES = 1.0
MIN_T = 2.0


def fit_power_law_with_survival_std(Nf, DP, stdnum=0.0):
    """Fit DP = A*Nf**b while avoiding inversion of flat or shallow curves.

    STD is the scatter expressed in ln(N); it is undefined for fits that are
    flat or increasing. STD_ln_DP is always evaluated in ln(DP). When the
    slope is negative, shifting DP downwards is algebraically the same as the
    conventional life shift. Flat and positive slopes follow that same
    downward DP-space rule, and must not be read as a life survival
    probability.
    """
    Nf = np.asarray(Nf, dtype=float).ravel()
    DP = np.asarray(DP, dtype=float).ravel()
    if Nf.shape != DP.shape:
        raise ValueError("the life and damage arrays must have the same length")
    if not np.isfinite(stdnum) or stdnum < 0:
        raise ValueError("the number of standard deviations must be finite and nonnegative")
    keep = np.isfinite(Nf) & np.isfinite(DP) & (Nf > 0) & (DP > 0)
    Nf, DP = Nf[keep], DP[keep]
    degenerate = dict(
        A=np.nan,
        b=np.nan,
        STD=np.nan,
        A_surv=np.nan,
        STD_ln_DP=np.nan,
        n=int(Nf.size),
        R2_log=np.nan,
        decades=float(np.log10(Nf.max() / Nf.min())) if Nf.size >= 2 else 0.0,
        b_se=np.nan,
        identifiable=False,
        reason=f"{Nf.size} tests with distinct positive lives",
    )
    if Nf.size < 2 or np.ptp(np.log(Nf)) < 1e-12:
        return degenerate
    log_N, log_DP = np.log(Nf), np.log(DP)
    slope, intercept = np.polyfit(log_N, log_DP, 1)
    if abs(slope) < 1e-12:
        slope, intercept = 0.0, float(np.mean(log_DP))
    leftover = log_DP - (intercept + slope * log_N)
    n = int(Nf.size)
    sse = float(np.sum(leftover**2))
    # two fitted parameters: n - 2 degrees of freedom (unbiased variance)
    sigma_dp = float(np.sqrt(sse / (n - 2))) if n > 2 else np.nan
    sigma_life = sigma_dp / abs(slope) if slope < 0 else np.nan
    ss_total = float(np.sum((log_DP - np.mean(log_DP)) ** 2))
    r2 = 1 - sse / ss_total if ss_total > 1e-24 else np.nan
    sxx = float(np.sum((log_N - log_N.mean()) ** 2))
    b_se = sigma_dp / np.sqrt(sxx) if n > 2 else np.nan
    decades = float(np.log10(Nf.max() / Nf.min()))
    issues = []
    if n < MIN_POINTS:
        issues.append(f"{n} < {MIN_POINTS} design tests")
    if decades < MIN_DECADES:
        issues.append(f"{decades:.2f} < {MIN_DECADES} decades of life")
    if not slope < 0.0 or not abs(slope) > MIN_T * b_se:
        issues.append(f"slope {slope:.4f} not distinct from 0 (SE {b_se:.4f})")
    shift = stdnum * sigma_dp if stdnum > 0 else 0.0
    return dict(
        A=float(np.exp(intercept)),
        b=float(slope),
        STD=sigma_life,
        STD_ln_DP=sigma_dp,
        A_surv=float(np.exp(intercept - shift)),
        n=n,
        R2_log=r2,
        decades=decades,
        b_se=float(b_se),
        identifiable=not issues,
        reason="; ".join(issues),
    )


def stdnum_for_survival(P_s: float) -> float:
    """Number of standard deviations z = Phi^-1(P_s) for a survival probability P_s.

    Statistical model of the master curve: ln DP = ln A + b ln N + e, with e normal and the
    standard deviation estimated from the design residuals (STD_ln_DP). For a decreasing curve
    (b < 0) the curve shifted down by z STD_ln_DP gives, at a given DP, the life
    N_Ps = N_50 exp(-z STD), i.e. the life that a fraction P_s of specimens exceeds under the same
    assumption with STD = STD_ln_DP / |b|. P_s is a survival probability of the population, not a
    confidence level: the scatter is estimated from few tests and no tolerance factor is applied.
    P_s = 0.5 gives the median curve (z = 0). Only 0.5 <= P_s < 1 is accepted, since the shift is
    defined downwards.
    """
    from statistics import NormalDist

    P_s = float(P_s)
    if not (np.isfinite(P_s) and 0.5 <= P_s < 1.0):
        raise ValueError(f"the survival probability P_s must satisfy 0.5 <= P_s < 1; got {P_s}")
    return float(NormalDist().inv_cdf(P_s))
