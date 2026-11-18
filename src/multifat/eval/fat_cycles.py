# multifat/eval/fat_cycles.py
import numpy as np
from typing import Optional, Tuple, Union
import matplotlib.pyplot as plt


def fatigue_cycles(
    F_par: Union[float, np.ndarray],
    A: float,
    B: float,
    Nmin: Optional[float] = None,
    Nmax: Optional[float] = None,
    *,
    return_warnings: bool = False,
) -> Union[float, np.ndarray, Tuple[Union[float, np.ndarray], str]]:
    """
    Invert a Basquin-type relation to obtain the expected number of cycles:

        F_par = A * Nexp^B  ->  Nexp = (F_par / A)^(1/B)

    Parameters
    ----------
    F_par : float or np.ndarray
        One or more fatigue parameter values (for instance FS_val, SWT_val, FIN_val).
    A : float
        Basquin coefficient A (must not be zero; normally positive).
    B : float
        Basquin exponent B (must not be zero).
    Nmin, Nmax : float or None, optional
        Optional limits on Nexp. When given and violated, a warning text is emitted.
        When None, the corresponding check is not performed.
    return_warnings : bool, default False
        When True, returns (Nexp, message). Otherwise only Nexp is returned.

    Returns
    -------
    Nexp : float or np.ndarray
        Expected number of cycles to failure.
    message : str (optional)
        Provided only when return_warnings=True.
    """
    # Wrap in an array so the computation works elementwise; scalars stay scalars on output
    F_arr = np.asarray(F_par, dtype=float)

    if A == 0:
        raise ValueError("The coefficient A must not be zero.")
    if B == 0:
        raise ValueError("The exponent B must not be zero.")

    # Real-valued results require F/A > 0 whenever 1/B is not an integer
    quotient = F_arr / A
    if np.any(quotient <= 0):
        raise ValueError("F_par/A has to be > 0 so that Nexp stays real-valued.")

    Nexp = quotient ** (1.0 / B)

    # Collect any bound violations
    notes = []
    if Nmin is not None:
        if np.any(Nexp < Nmin):
            notes.append(f"Nexp below Nmin ({Nmin}).")
    if Nmax is not None:
        if np.any(Nexp > Nmax):
            notes.append(f"Nexp above Nmax ({Nmax}).")
    message = " ".join(notes)

    # Give back a Python float when the caller passed a scalar
    if np.isscalar(F_par):
        Nexp = float(Nexp)

    if return_warnings:
        return Nexp, message
    return Nexp


def plot_runouts(ax, x, y, *, color="C0", label="Run-out (recorded cycles: a lower bound)"):
    """Hollow circles with rightward arrows on a cycles axis (upstream 7a87256)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if not keep.any():
        return
    ax.scatter(x[keep], y[keep], facecolors="none", edgecolors=color, marker="o", label=label)
    for x_i, y_i in zip(x[keep], y[keep]):
        ax.annotate(
            "",
            xy=(x_i * 1.35, y_i),
            xytext=(x_i, y_i),
            arrowprops={"arrowstyle": "->", "color": color, "lw": 0.8},
        )


def fatigue_plots(
    F_par,
    N_exp,
    N_pred,
    A: float,
    B: float,
    *,
    # NEW: optional “design points” (plotted in different color)
    F_design=None,
    N_design=None,
    loglog: bool = True,
    n_curve: int = 400,
    curve_range=None,  # (N_low, N_high) optional
    error_factors=(2, 5),  # NEW: factor bands for parity plot
    show: bool = True,
    return_figax: bool = True,
    # Run-outs: hollow markers with arrows (the recorded cycle count is a lower bound), never fitted
    F_runout=None,
    N_runout=None,
    N_pred_runout=None,
    # Master curve DP(N) to draw instead of A N^B (non-Basquin S-N models, issue #38)
    curve=None,
):
    """
    Build two fatigue diagrams:
      1) Wöhler plane: the measured points (N_exp, F_par) together with the design Basquin
         curve F = A*N^B, plus the optional design points (N_design, F_design) in a distinct
         color, with the curve drawn in that same color.
      2) Parity diagram: N_exp against N_pred as a scatter, a black y = x diagonal,
         thin dashed bands at factors of (2, 5, ...), equal axis scaling, and the grid
         kept underneath the data.

    Parameters
    ----------
    F_par : array-like
        Measured fatigue parameter values (same length as N_exp and N_pred).
    N_exp : array-like
        Measured cycles to failure.
    N_pred : array-like
        Predicted cycles to failure.
    A, B : float
        Coefficients of the design curve F = A * N^B.
    F_design, N_design : array-like or None
        Optional design points to draw (equal lengths). Skipped when None.
    loglog : bool
        When True, logarithmic axes are used (recommended).
    n_curve : int
        Number of samples used to draw the design curve.
    curve_range : tuple(N_low, N_high) or None
        Cycle range spanned by the curve. When None it is inferred from the data
        (N_exp, N_pred and, if supplied, N_design).
    error_factors : tuple
        Multiplicative factors drawn as dashed bands in the parity diagram (for example 2, 5).
    show : bool
        When True, plt.show() is called.
    return_figax : bool
        When True, (fig1, ax1, fig2, ax2) is returned.

    Returns
    -------
    (fig1, ax1, fig2, ax2) when return_figax=True, otherwise None.
    """
    F_par = np.asarray(F_par, dtype=float).ravel()
    N_exp = np.asarray(N_exp, dtype=float).ravel()
    N_pred = np.asarray(N_pred, dtype=float).ravel()

    if not (len(F_par) == len(N_exp) == len(N_pred)):
        raise ValueError("F_par, N_exp, and N_pred must all have equal length.")
    if A == 0 or B == 0:
        raise ValueError("Neither A nor B may be zero.")
    if np.any(N_exp <= 0) or np.any(N_pred <= 0):
        raise ValueError("N_exp and N_pred must be > 0 to be plotted on log axes.")
    F_ro = np.asarray(F_runout if F_runout is not None else [], dtype=float).ravel()
    N_ro = np.asarray(N_runout if N_runout is not None else [], dtype=float).ravel()
    Np_ro = np.asarray(N_pred_runout if N_pred_runout is not None else [], dtype=float).ravel()
    good = np.isfinite(N_pred) & (N_pred > 0)
    if loglog and np.any(F_par <= 0):
        raise ValueError("F_par must be > 0 to plot the Wöhler plane on log-log axes.")

    has_design = (F_design is not None) and (N_design is not None)
    if has_design:
        F_design = np.asarray(F_design, dtype=float).ravel()
        N_design = np.asarray(N_design, dtype=float).ravel()
        if len(F_design) != len(N_design):
            raise ValueError("F_design and N_design must have equal length.")
        if np.any(N_design <= 0):
            raise ValueError("N_design must be > 0 to be plotted on log axes.")
        if loglog and np.any(F_design <= 0):
            raise ValueError("F_design must be > 0 to plot the Wöhler plane on log-log axes.")

    # --- Work out the extent of the curve ---
    if curve_range is None:
        pooled = [N_exp, N_pred[good], N_ro[np.isfinite(N_ro) & (N_ro > 0)]]
        if has_design:
            pooled.append(N_design)
        pooled = np.concatenate(pooled)
        c_lo = np.min(pooled) / 10.0
        c_hi = np.max(pooled) * 10.0
    else:
        c_lo, c_hi = curve_range

    c_lo = float(max(c_lo, 1e-300))
    c_hi = float(max(c_hi, c_lo * 1.0001))

    if loglog:
        N_grid = np.logspace(np.log10(c_lo), np.log10(c_hi), n_curve)
    else:
        N_grid = np.linspace(c_lo, c_hi, n_curve)

    F_curve = curve(N_grid) if curve is not None else A * (N_grid**B)

    # =========================
    # 1) Wöhler plane
    # =========================
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    # Keep the grid behind every artist
    ax1.set_axisbelow(True)

    # Measured points come first
    ax1.scatter(N_exp, F_par, label="Experimental")
    plot_runouts(ax1, N_ro, F_ro)

    # Color the curve like the design points when those exist;
    # otherwise let matplotlib pick the default color.
    curve_style = {}
    if has_design:
        # Draw the design points first so the curve can reuse their default color
        design_scatter = ax1.scatter(N_design, F_design, label="Design point(s)")
        curve_style["color"] = design_scatter.get_facecolor()[0]  # matches the design points

    # Draw the design curve (same color as the design points, when present)
    if loglog:
        ax1.loglog(N_grid, F_curve, label="Design curve", **curve_style)
    else:
        ax1.plot(N_grid, F_curve, label="Design curve", **curve_style)

    ax1.set_title("Wöhler plane: experimental + design")
    ax1.set_xlabel("N (cycles)")
    ax1.set_ylabel("Fatigue parameter F")

    ax1.grid(True, which="both")
    ax1.legend()

    # =========================
    # 2) Nexp vs Npred (parity)
    # =========================
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    ax2.set_axisbelow(True)  # grid beneath points and lines
    ax2.scatter(N_exp[good], N_pred[good], label="Predictions")
    if len(Np_ro) == len(N_ro):
        plot_runouts(ax2, N_ro, Np_ro)
    if not good.any():
        ax2.text(
            0.5,
            0.5,
            "Life predictions are unavailable for this calibration",
            transform=ax2.transAxes,
            ha="center",
            va="center",
        )

    # Extent of the parity diagram
    both = np.concatenate([N_exp, N_pred[good]])
    p_lo = np.min(both) / 10.0
    p_hi = np.max(both) * 10.0
    p_lo = float(max(p_lo, 1e-300))
    p_hi = float(max(p_hi, p_lo * 1.0001))

    if loglog:
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        diag = np.logspace(np.log10(p_lo), np.log10(p_hi), 300)
    else:
        diag = np.linspace(p_lo, p_hi, 300)

    # Black y = x diagonal
    ax2.plot(diag, diag, color="black", linestyle="--", linewidth=1.5, label="y = x")

    # Thin dashed factor bands (2 and 5 unless overridden)
    for fac in error_factors:
        fac = float(fac)
        ax2.plot(diag, fac * diag, linestyle="--", linewidth=0.8, label=f"×{fac:g}")
        ax2.plot(diag, diag / fac, linestyle="--", linewidth=0.8, label=f"÷{fac:g}")

    ax2.set_title("Parity plot: N experimental vs N expected")
    ax2.set_xlabel("N experimental (cycles)")
    ax2.set_ylabel("N expected (cycles)")

    ax2.grid(True, which="both")

    # Match the two spans and then force an equal aspect ratio
    ax2.set_xlim(p_lo, p_hi)
    ax2.set_ylim(p_lo, p_hi)
    ax2.set_aspect("equal", adjustable="box")

    ax2.legend()

    if show:
        plt.show()

    if return_figax:
        return fig1, ax1, fig2, ax2
    return None
