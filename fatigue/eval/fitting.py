# fatigue/eval/fitting.py
from __future__ import annotations
import numpy as np
from typing import Dict, Any, Tuple

def fit_power_law_with_survival_std(N: np.ndarray, CP: np.ndarray, stdnum: float = 0.0) -> Dict[str, float]:
    """
    Fit CP = A * N^b (median) and return a 'survival shifted' intercept using stdnum.

    Residual definition matches your old script:
      1) Fit ln(CP) = ln(A) + b ln(N)
      2) N_fit = (CP/A)^(1/b)
      3) residual = ln(N) - ln(N_fit)
      4) STD = std(residual)
      5) ln(A_surv) = ln(A) + b * stdnum * STD  -> A_surv = A * exp(b*stdnum*STD)
    """
    N = np.asarray(N, dtype=float)
    CP = np.asarray(CP, dtype=float)

    m = np.isfinite(N) & np.isfinite(CP) & (N > 0) & (CP > 0)
    N = N[m]
    CP = CP[m]

    if N.size < 2:
        return {"A": np.nan, "b": np.nan, "STD": np.nan, "A_surv": np.nan, "n": int(N.size), "R2_log": np.nan}

    x = np.log(N)
    y = np.log(CP)

    # y = c + b x
    b, c = np.polyfit(x, y, 1)
    A = float(np.exp(c))

    # residuals in ln(N) space (your old approach)
    N_fit = (CP / A) ** (1.0 / b)
    residual = np.log(N) - np.log(N_fit)
    STD = float(np.std(residual, ddof=1)) if residual.size > 1 else 0.0

    # survival shift
    A_surv = float(np.exp(c + b * stdnum * STD))

    # R2 in log-space (median fit)
    yhat = c + b * x
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {"A": A, "b": float(b), "STD": STD, "A_surv": A_surv, "n": int(N.size), "R2_log": float(R2)}

