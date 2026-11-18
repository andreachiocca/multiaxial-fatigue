"""Master-curve (S-N) models relating a damage parameter DP to the life N (issue #38).

Every master-curve method converts its DP to a life through a curve fitted on the design tests
(uniaxial tension, R = -1). The models, all fitted by least squares on ln DP::

    basquin            DP = A N^b                                  Basquin (1910)
    kohout_vechet      DP = DP_inf (N / (N + C))^b                 Kohout & Vechet (2001), Eq. (H2)
    stromeyer          DP = a N^b + DP_inf                         Stromeyer (1914)

Kohout & Vechet (2001), Int J Fatigue 23:175-183, doi:10.1016/S0142-1123(00)00082-7: (H2) is the
finite-life and high-cycle form (Hx) of their function, written with the fatigue limit
DP_inf = a C^b (Eq. 10). It has the Basquin line DP = a N^b (a = DP_inf C^-b) as
asymptote for N << C, the horizontal asymptote DP_inf for N >> C, and slope b/2 in log-log
coordinates at N = C. Its inverse is N = u C / (1 - u), u = (DP / DP_inf)^(1/b), and the life is
infinite for DP <= DP_inf.

Repository choices:
  - basquin keeps the established linear fit of ln DP on ln N (``fit_power_law_with_survival_std``),
    bit for bit; it is the default;
  - the nonlinear models are fitted on ln DP (multiplicative scatter, as for Basquin) with
    ``scipy.optimize.least_squares`` (or bounded MLE with run-outs via ``scipy.optimize.minimize``)
    from several starting points derived from the Basquin fit;
    the paper does not state its regression space;
  - Kohout & Vechet and Stromeyer originally formulated their models for uniaxial stress sigma;
    in multifat the master curve maps an arbitrary damage parameter DP to life, so the fatigue-limit
    asymptote is DP_inf with the units of the evaluated DP; sigma_inf is kept as a backward-compatible alias;
  - run-outs in the design set are incorporated via a Tobit-type censored log-likelihood,
    preventing infinite-life artifacts on low-stress failures;
  - the whole-range form (W1/W2) of Kohout & Vechet needs low-cycle data or the strength sigma_1,
    which has no counterpart for an arbitrary DP: it is not offered;
  - a fit is rejected, with its reason, when there are fewer than 2 p design tests, the solver does
    not converge, a parameter stops on a search bound, or the parameter covariance is singular;
  - the survival shift multiplies the whole curve by exp(-z s), s = sqrt(SSE / (n - p)) in ln DP:
    for a curved model the corresponding life factor varies along the curve.
"""

from __future__ import annotations

import math

import numpy as np

from multifat.eval import fat_cycles as fatc
from multifat.eval.fitting import fit_power_law_with_survival_std

SN_MODELS = ("basquin", "kohout_vechet", "stromeyer")
N_PARAMS = {"basquin": 2, "kohout_vechet": 3, "stromeyer": 3}
# Search intervals of the nonlinear fits: an optimum on one of these is reported as unidentified
B_BOUNDS = (-3.0, -1e-6)
LN_C_BOUNDS = (math.log(1e2), math.log(1e12))


def fit_sn_curve(
    Nf,
    DP,
    stdnum: float = 0.0,
    model: str = "basquin",
    N_runout=None,
    DP_runout=None,
) -> dict:
    """Fit the master curve `model`; the dict carries the parameters, scatter and validity."""
    if model not in SN_MODELS:
        raise ValueError(f"unknown S-N model {model!r}; choose one of {SN_MODELS}")
    base = fit_power_law_with_survival_std(Nf, DP, stdnum=stdnum)
    if model == "basquin":
        return {**base, "sn_model": "basquin", "valid": True, "reason": ""}
    Nf = np.asarray(Nf, dtype=float).ravel()
    DP = np.asarray(DP, dtype=float).ravel()
    keep = np.isfinite(Nf) & np.isfinite(DP) & (Nf > 0) & (DP > 0)
    Nf, DP = Nf[keep], DP[keep]

    if N_runout is not None and DP_runout is not None:
        N_r = np.asarray(N_runout, dtype=float).ravel()
        DP_r = np.asarray(DP_runout, dtype=float).ravel()
        keep_r = np.isfinite(N_r) & np.isfinite(DP_r) & (N_r > 0) & (DP_r > 0)
        N_r, DP_r = N_r[keep_r], DP_r[keep_r]
    else:
        N_r = np.empty(0, dtype=float)
        DP_r = np.empty(0, dtype=float)

    p = N_PARAMS[model]
    out = {
        "sn_model": model,
        "n": int(Nf.size),
        "n_runout": int(N_r.size),
        "A": np.nan,
        "b": np.nan,
        "A_surv": np.nan,
        "STD": np.nan,
        "STD_ln_DP": np.nan,
        "R2_log": np.nan,
        "params": {},
        "param_se": {},
        "shift": np.nan,
        "valid": False,
        "reason": "",
    }
    if Nf.size < 2 * p:
        out["reason"] = f"{model} needs at least {2 * p} design tests, got {Nf.size}"
        return out
    if not (np.isfinite(base["b"]) and base["b"] < 0):
        out["reason"] = f"{model}: the Basquin starting fit has no negative slope"
        return out
    solution = _fit_nonlinear(model, Nf, DP, base, N_r, DP_r)
    if isinstance(solution, str):
        out["reason"] = solution
        return out
    theta, residual, se, s = solution
    if not np.all(np.isfinite(se)):
        out["reason"] = f"{model}: singular parameter covariance (parameters not identifiable)"
        return out
    log_DP = np.log(DP)
    ss_total = float(np.sum((log_DP - log_DP.mean()) ** 2))
    shift = math.exp(-stdnum * s)
    if model == "kohout_vechet":
        ln_dp_inf, ln_C, b = theta
        dp_inf = math.exp(ln_dp_inf)
        params = {
            "DP_inf": dp_inf,
            "sigma_inf": dp_inf,
            "C": math.exp(ln_C),
            "b": float(b),
        }
        names = ("ln_DP_inf", "ln_C", "b")
        A = dp_inf * params["C"] ** (-params["b"])  # Basquin asymptote a
    else:
        ln_a, b, ln_dp_inf = theta
        dp_inf = math.exp(ln_dp_inf)
        params = {
            "a": math.exp(ln_a),
            "b": float(b),
            "DP_inf": dp_inf,
            "sigma_inf": dp_inf,
        }
        names = ("ln_a", "b", "ln_DP_inf")
        A = params["a"]
    param_se = dict(zip(names, (float(x) for x in se)))
    if "ln_DP_inf" in param_se:
        param_se["ln_sigma_inf"] = param_se["ln_DP_inf"]
    out.update(
        A=float(A),
        b=float(params["b"]),
        A_surv=float(A * shift),
        STD_ln_DP=s,
        R2_log=1 - float(residual @ residual) / ss_total if ss_total > 1e-24 else np.nan,
        params=params,
        param_se=param_se,
        shift=float(shift),
        valid=True,
    )
    return out


def _model_ln(model: str, theta, Nf):
    if model == "kohout_vechet":
        ln_dp_inf, ln_C, b = theta
        return ln_dp_inf + b * (np.log(Nf) - np.logaddexp(np.log(Nf), ln_C))
    ln_a, b, ln_dp_inf = theta
    return np.logaddexp(ln_a + b * np.log(Nf), ln_dp_inf)


def _fit_nonlinear(model: str, Nf, DP, base, N_r=None, DP_r=None):
    """Best solution over several starts, or the reason why there is none."""
    from scipy.optimize import least_squares, minimize

    log_DP = np.log(DP)
    ln_A0, b0 = math.log(base["A"]), float(base["b"])
    ln_floor = math.log(float(DP.min())) - 1e-6
    p = N_PARAMS[model]
    dof = Nf.size - p

    if model == "kohout_vechet":
        bounds = ([-np.inf, LN_C_BOUNDS[0], B_BOUNDS[0]], [ln_floor, LN_C_BOUNDS[1], B_BOUNDS[1]])
        starts = [
            (min(ln_A0 + b0 * ln_C, ln_floor - 0.1), ln_C, b0)
            for ln_C in np.linspace(LN_C_BOUNDS[0] + 1.0, LN_C_BOUNDS[1] - 1.0, 9)
        ]
    else:
        bounds = ([-np.inf, B_BOUNDS[0], -np.inf], [np.inf, B_BOUNDS[1], ln_floor])
        starts = [(ln_A0, b0, ln_floor + shift) for shift in np.linspace(-4.0, -0.05, 8)]

    best = None
    for start in starts:
        start = np.clip(start, np.asarray(bounds[0]) + 1e-9, np.asarray(bounds[1]) - 1e-9)
        result = least_squares(
            lambda theta: _model_ln(model, theta, Nf) - log_DP,
            start,
            bounds=bounds,
            method="trf",
            x_scale="jac",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=2000,
        )
        if result.status > 0 and (best is None or result.cost < best.cost):
            best = result
    if best is None:
        return f"{model}: the least-squares solver did not converge"

    if model == "kohout_vechet":
        b_val, c_val = best.x[2], best.x[1]
        at_bound = (
            np.isclose(b_val, B_BOUNDS[0], atol=1e-6)
            or np.isclose(b_val, B_BOUNDS[1], atol=1e-6)
            or np.isclose(c_val, LN_C_BOUNDS[0], atol=1e-6)
            or np.isclose(c_val, LN_C_BOUNDS[1], atol=1e-6)
        )
    else:
        b_val = best.x[1]
        at_bound = np.isclose(b_val, B_BOUNDS[0], atol=1e-6) or np.isclose(
            b_val, B_BOUNDS[1], atol=1e-6
        )
    if at_bound:
        return (
            f"{model}: a parameter stopped on its search bound (not identifiable from these tests)"
        )

    if N_r is None or N_r.size == 0:
        residual = best.fun
        s = float(math.sqrt(float(residual @ residual) / dof))
        try:
            cov = s**2 * np.linalg.inv(best.jac.T @ best.jac)
            se = np.sqrt(np.diag(cov))
        except np.linalg.LinAlgError:
            se = np.full(p, np.nan)
        return best.x, residual, se, s

    # Censored maximum likelihood with run-outs (Tobit likelihood)
    from scipy.special import log_ndtr

    s_init = max(1e-4, math.sqrt(float(best.fun @ best.fun) / dof))
    init_params = np.append(best.x, math.log(s_init))

    if model == "kohout_vechet":
        opt_bounds = [
            (-np.inf, ln_floor),
            (LN_C_BOUNDS[0], LN_C_BOUNDS[1]),
            (B_BOUNDS[0], B_BOUNDS[1]),
            (math.log(1e-4), math.log(10.0)),
        ]
    else:
        opt_bounds = [
            (-np.inf, np.inf),
            (B_BOUNDS[0], B_BOUNDS[1]),
            (-np.inf, ln_floor),
            (math.log(1e-4), math.log(10.0)),
        ]

    def nll(params_all):
        theta = params_all[:3]
        s = math.exp(params_all[3])
        pred_f = _model_ln(model, theta, Nf)
        res_f = (log_DP - pred_f) / s
        ll_f = -0.5 * math.log(2 * math.pi) - math.log(s) - 0.5 * res_f**2
        pred_r = _model_ln(model, theta, N_r)
        z_r = (pred_r - np.log(DP_r)) / s
        ll_r = log_ndtr(z_r)
        return -(float(np.sum(ll_f)) + float(np.sum(ll_r)))

    res_opt = minimize(nll, init_params, method="L-BFGS-B", bounds=opt_bounds)
    if not res_opt.success:
        return f"{model}: censored likelihood optimization did not converge"

    theta_opt = res_opt.x[:3]
    s_opt = float(math.exp(res_opt.x[3]))
    residual = _model_ln(model, theta_opt, Nf) - log_DP

    # Numerical Hessian for parameter standard errors
    eps = 1e-4
    n_p = 3
    hess = np.zeros((n_p, n_p))
    f0 = nll(res_opt.x)
    for i in range(n_p):
        xp = res_opt.x.copy()
        xp[i] += eps
        xm = res_opt.x.copy()
        xm[i] -= eps
        hess[i, i] = (nll(xp) - 2 * f0 + nll(xm)) / (eps**2)
        for j in range(i + 1, n_p):
            xpp = res_opt.x.copy()
            xpp[i] += eps
            xpp[j] += eps
            xpm = res_opt.x.copy()
            xpm[i] += eps
            xpm[j] -= eps
            xmp = res_opt.x.copy()
            xmp[i] -= eps
            xmp[j] += eps
            xmm = res_opt.x.copy()
            xmm[i] -= eps
            xmm[j] -= eps
            hess[i, j] = (nll(xpp) - nll(xpm) - nll(xmp) + nll(xmm)) / (4 * eps**2)
            hess[j, i] = hess[i, j]

    try:
        cov = np.linalg.inv(hess)
        se = np.sqrt(np.diag(cov))
    except Exception:
        se = np.full(n_p, np.nan)

    return theta_opt, residual, se, s_opt


def curve_dp(fit: dict, Nf):
    """DP of the (shifted) master curve at the lives Nf; NaN where Nf is not positive and finite."""
    Nf = np.asarray(Nf, dtype=float)
    out = np.full(Nf.shape, np.nan)
    ok = np.isfinite(Nf) & (Nf > 0)
    model = fit.get("sn_model", "basquin")
    if model == "basquin":
        A_surv, b = float(fit["A_surv"]), float(fit["b"])
        if np.isfinite(A_surv) and np.isfinite(b):
            out[ok] = A_surv * (Nf[ok] ** b)
        return out
    if not fit.get("valid"):
        return out
    p, shift = fit["params"], fit["shift"]
    dp_inf = p.get("DP_inf", p.get("sigma_inf"))
    if model == "kohout_vechet":
        out[ok] = shift * dp_inf * (Nf[ok] / (Nf[ok] + p["C"])) ** p["b"]
    else:
        out[ok] = shift * (p["a"] * Nf[ok] ** p["b"] + dp_inf)
    return out


def curve_life(fit: dict, DP):
    """Life on the (shifted) master curve for DP; +inf at or below a fatigue-limit asymptote."""
    model = fit.get("sn_model", "basquin")
    if model == "basquin":
        return fatc.fatigue_cycles(DP, fit["A_surv"], fit["b"])
    DP = np.asarray(DP, dtype=float)
    out = np.full(DP.shape, np.nan)
    if not fit.get("valid"):
        return out
    p, shift = fit["params"], fit["shift"]
    dp_inf = p.get("DP_inf", p.get("sigma_inf"))
    limit = shift * dp_inf
    ok = np.isfinite(DP) & (DP > 0)
    below = ok & (DP <= limit)
    above = ok & (DP > limit)
    out[below] = np.inf
    if model == "kohout_vechet":
        u = (DP[above] / limit) ** (1.0 / p["b"])
        out[above] = u * p["C"] / (1.0 - u)
    else:
        out[above] = ((DP[above] / shift - dp_inf) / p["a"]) ** (1.0 / p["b"])
    return out
