"""Power-law calibration with a stable damage-space survival shift."""
from __future__ import annotations
import numpy as np


def fit_power_law_with_survival_std(N, CP, stdnum=0.0):
    """Fit DP=A*N**b without inverting flat or shallow curves.

    STD is scatter in ln(N), undefined for non-decreasing/flat fits.
    STD_ln_DP is always computed in ln(DP). For negative slopes the
    downward DP shift is algebraically identical to the legacy life shift.
    Flat/positive slopes use the same downward DP-space convention; they
    must not be interpreted as a life survival probability.
    """
    N, CP = np.asarray(N, dtype=float).ravel(), np.asarray(CP, dtype=float).ravel()
    if N.shape != CP.shape:
        raise ValueError("N and CP must have equal lengths")
    if not np.isfinite(stdnum) or stdnum < 0:
        raise ValueError("stdnum must be finite and nonnegative")
    mask = np.isfinite(N) & np.isfinite(CP) & (N > 0) & (CP > 0)
    N, CP = N[mask], CP[mask]
    invalid = dict(A=np.nan, b=np.nan, STD=np.nan, A_surv=np.nan,
                   STD_ln_DP=np.nan, n=int(N.size), R2_log=np.nan)
    if N.size < 2 or np.ptp(np.log(N)) < 1e-12:
        return invalid
    x, y = np.log(N), np.log(CP)
    b, c = np.polyfit(x, y, 1)
    if abs(b) < 1e-12:
        b, c = 0.0, float(np.mean(y))
    residual_dp = y-(c+b*x)
    std_dp = float(np.std(residual_dp, ddof=1))
    std_life = std_dp/abs(b) if b < 0 else np.nan
    ss_tot = float(np.sum((y-np.mean(y))**2))
    r2 = 1-float(np.sum(residual_dp**2))/ss_tot if ss_tot > 1e-24 else np.nan
    return dict(A=float(np.exp(c)), b=float(b), STD=std_life,
                STD_ln_DP=std_dp, A_surv=float(np.exp(c-stdnum*std_dp)),
                n=int(N.size), R2_log=r2)
