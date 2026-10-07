# fatigue/eval/fat_cycles.py
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
    Compute expected cycles to failure by inverting a Basquin-like law:

        F_par = A * Nexp^B  ->  Nexp = (F_par / A)^(1/B)

    Parameters
    ----------
    F_par : float or np.ndarray
        Fatigue parameter value(s) (e.g., FS_val, SWT_val, FIN_val).
    A : float
        Basquin coefficient A (must be non-zero; typically > 0).
    B : float
        Basquin exponent B (must be non-zero).
    Nmin, Nmax : float or None, optional
        Optional bounds for Nexp. If provided and exceeded, a warning message is produced.
        If None, that check is skipped.
    return_warnings : bool, default False
        If True, returns (Nexp, message). Otherwise returns only Nexp.

    Returns
    -------
    Nexp : float or np.ndarray
        Expected number of cycles to failure.
    message : str (optional)
        Only returned if return_warnings=True.
    """
    # Convert to numpy for vectorized behavior; keep scalar if scalar input
    F = np.asarray(F_par, dtype=float)

    if A == 0:
        raise ValueError("A must be non-zero.")
    if B == 0:
        raise ValueError("B must be non-zero.")

    # Domain check: to stay in reals, F/A should be > 0 for arbitrary (non-integer) 1/B
    ratio = F / A
    if np.any(ratio <= 0):
        raise ValueError("F_par/A must be > 0 to compute real-valued Nexp.")

    Nexp = ratio ** (1.0 / B)

    # Prepare optional warnings
    msg_parts = []
    if Nmin is not None:
        if np.any(Nexp < Nmin):
            msg_parts.append(f"Nexp below Nmin ({Nmin}).")
    if Nmax is not None:
        if np.any(Nexp > Nmax):
            msg_parts.append(f"Nexp above Nmax ({Nmax}).")
    msg = " ".join(msg_parts)

    # Return scalar if input was scalar
    if np.isscalar(F_par):
        Nexp = float(Nexp)

    if return_warnings:
        return Nexp, msg
    return Nexp



#def fatigue_cycles_plt(
#    F_par,
#    A: float,
#    B: float,
#    Nmin=None,
#    Nmax=None,
#    *,
#    return_warnings: bool = False,
#    plot: bool = False,
#    ax=None,
#    plot_range=None,      # (N_low, N_high) optional
#    n_plot: int = 300,
#    loglog: bool = True,
#    show: bool = True,
#    return_figax: bool = False,
#):
#    """
#    Invert Basquin-like law:
#        F_par = A * Nexp^B  ->  Nexp = (F_par / A)^(1/B)
#
#    Optional:
#      - bounds check via Nmin/Nmax
#      - plot curve and point(s)
#
#    Parameters
#    ----------
#    F_par : float or array-like
#        Fatigue parameter value(s).
#    A, B : float
#        Basquin coefficients (A != 0, B != 0).
#    Nmin, Nmax : float or None
#        Optional bounds for Nexp. If None, that check is skipped.
#    return_warnings : bool
#        If True returns (Nexp, message).
#    plot : bool
#        If True, plots the fatigue curve and the computed point(s).
#    ax : matplotlib axis or None
#        If provided, plot on this axis; otherwise create a new figure/axis.
#    plot_range : tuple(N_low, N_high) or None
#        If provided, uses this N range for the curve.
#        Else uses (Nmin,Nmax) if both provided, otherwise a range around Nexp.
#    n_plot : int
#        Number of samples for the curve.
#    loglog : bool
#        If True uses log-log scale.
#    show : bool
#        If True and a new figure is created, calls plt.show().
#    return_figax : bool
#        If True (and plot=True), also returns (fig, ax).
#
#    Returns
#    -------
#    Nexp : float or np.ndarray
#    message : str (optional)
#    fig, ax : (optional, if plot=True and return_figax=True)
#    """
#    F_in = F_par
#    F = np.asarray(F_par, dtype=float)
#
#    if A == 0:
#        raise ValueError("A must be non-zero.")
#    if B == 0:
#        raise ValueError("B must be non-zero.")
#
#    ratio = F / A
#    if np.any(ratio <= 0):
#        raise ValueError("F_par/A must be > 0 to compute real-valued Nexp.")
#
#    Nexp = ratio ** (1.0 / B)
#
#    # Optional warnings
#    msg_parts = []
#    if Nmin is not None and np.any(Nexp < Nmin):
#        msg_parts.append(f"Nexp below Nmin ({Nmin}).")
#    if Nmax is not None and np.any(Nexp > Nmax):
#        msg_parts.append(f"Nexp above Nmax ({Nmax}).")
#    msg = " ".join(msg_parts)
#
#    # Convert scalar back to Python float if scalar input
#    if np.isscalar(F_in):
#        Nexp_out = float(Nexp)
#    else:
#        Nexp_out = Nexp
#
#    fig = None
#    if plot:
#        created = False
#        if ax is None:
#            fig = plt.figure()
#            ax = fig.add_subplot(111)
#            created = True
#        else:
#            fig = ax.figure
#
#        # Determine plotting N range
#        if plot_range is not None:
#            N_low, N_high = plot_range
#        elif (Nmin is not None) and (Nmax is not None):
#            N_low, N_high = Nmin, Nmax
#        else:
#            nmin = float(np.min(Nexp))
#            nmax = float(np.max(Nexp))
#            if nmin <= 0 or nmax <= 0:
#                # should not happen given checks, but keep safe
#                nmin, nmax = 1.0, 10.0
#            if np.isclose(nmin, nmax):
#                # single point -> pick a couple decades around it
#                N_low, N_high = nmin / 100.0, nmax * 100.0
#            else:
#                N_low, N_high = nmin / 10.0, nmax * 10.0
#
#        # Guard against bad ranges
#        N_low = max(float(N_low), 1e-300)
#        N_high = max(float(N_high), N_low * 1.0001)
#
#        if loglog:
#            N_line = np.logspace(np.log10(N_low), np.log10(N_high), n_plot)
#        else:
#            N_line = np.linspace(N_low, N_high, n_plot)
#
#        F_line = A * (N_line ** B)
#
#        if loglog:
#            ax.loglog(N_line, F_line)
#        else:
#            ax.plot(N_line, F_line)
#
#        # Plot computed points
#        ax.scatter(np.ravel(Nexp), np.ravel(F), s=30)
#
#        ax.set_xlabel("N (cycles)")
#        ax.set_ylabel("F_par")
#        ax.set_title("Basquin curve and computed point(s)")
#        ax.grid(True, which="both")
#
#        if created and show:
#            plt.show()
#
#    # Return logic
#    if return_warnings and plot and return_figax:
#        return Nexp_out, msg, fig, ax
#    if return_warnings:
#        return Nexp_out, msg
#    if plot and return_figax:
#        return Nexp_out, fig, ax
#    return Nexp_out



import numpy as np
import matplotlib.pyplot as plt


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
    curve_range=None,          # (N_low, N_high) optional
    error_factors=(2, 5),      # NEW: factor bands for parity plot
    show: bool = True,
    return_figax: bool = True,
):
    """
    Create two fatigue plots:
      1) Wöhler plane: experimental points (N_exp, F_par) + design Basquin curve F=A*N^B
         plus optional design points (N_design, F_design) in a different color,
         and the curve drawn with the same color as the design points.
      2) Parity plot: N_exp vs N_pred scatter + black diagonal y=x
         plus thin dashed factor-of-(2,5,...) bands, equal axis aspect, grid under data.

    Parameters
    ----------
    F_par : array-like
        Experimental fatigue parameter values (same length as N_exp and N_pred).
    N_exp : array-like
        Experimental cycles to failure.
    N_pred : array-like
        Predicted cycles to failure.
    A, B : float
        Design curve coefficients in F = A * N^B.
    F_design, N_design : array-like or None
        Optional design points to plot (same lengths). If None, skipped.
    loglog : bool
        If True uses log scales (recommended).
    n_curve : int
        Samples for the design curve.
    curve_range : tuple(N_low, N_high) or None
        Curve N-range. If None inferred from data (N_exp, N_pred, N_design if provided).
    error_factors : tuple
        Multiplicative error factors to draw as dashed bands in parity plot (e.g. 2,5).
    show : bool
        If True, calls plt.show().
    return_figax : bool
        If True, returns (fig1, ax1, fig2, ax2).

    Returns
    -------
    (fig1, ax1, fig2, ax2) if return_figax=True, else None.
    """
    F_par = np.asarray(F_par, dtype=float).ravel()
    N_exp = np.asarray(N_exp, dtype=float).ravel()
    N_pred = np.asarray(N_pred, dtype=float).ravel()

    if not (len(F_par) == len(N_exp) == len(N_pred)):
        raise ValueError("F_par, N_exp, and N_pred must have the same length.")
    if A == 0 or B == 0:
        raise ValueError("A and B must be non-zero.")
    if np.any(N_exp <= 0) or np.any(N_pred <= 0):
        raise ValueError("N_exp and N_pred must be > 0 for log plotting.")
    if loglog and np.any(F_par <= 0):
        raise ValueError("F_par must be > 0 for log-log Wöhler plotting.")

    have_design = (F_design is not None) and (N_design is not None)
    if have_design:
        F_design = np.asarray(F_design, dtype=float).ravel()
        N_design = np.asarray(N_design, dtype=float).ravel()
        if len(F_design) != len(N_design):
            raise ValueError("F_design and N_design must have the same length.")
        if np.any(N_design <= 0):
            raise ValueError("N_design must be > 0 for log plotting.")
        if loglog and np.any(F_design <= 0):
            raise ValueError("F_design must be > 0 for log-log Wöhler plotting.")

    # --- Determine curve range ---
    if curve_range is None:
        N_all = [N_exp, N_pred]
        if have_design:
            N_all.append(N_design)
        N_all = np.concatenate(N_all)
        N_low = np.min(N_all) / 10.0
        N_high = np.max(N_all) * 10.0
    else:
        N_low, N_high = curve_range

    N_low = float(max(N_low, 1e-300))
    N_high = float(max(N_high, N_low * 1.0001))

    if loglog:
        N_line = np.logspace(np.log10(N_low), np.log10(N_high), n_curve)
    else:
        N_line = np.linspace(N_low, N_high, n_curve)

    F_line = A * (N_line ** B)

    # =========================
    # 1) Wöhler plane
    # =========================
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    # Put grid below everything
    ax1.set_axisbelow(True)

    # Plot experimental points first
    if loglog:
        ax1.scatter(N_exp, F_par, label="Experimental")
    else:
        ax1.scatter(N_exp, F_par, label="Experimental")

    # Choose curve color = design points color (if present),
    # otherwise just plot the curve normally.
    curve_kwargs = {}
    if have_design:
        # Get matplotlib default next color by plotting design points first, then use that color
        des = ax1.scatter(N_design, F_design, label="Design point(s)")
        curve_kwargs["color"] = des.get_facecolor()[0]  # same as design points

    # Plot the design curve (same color as design points if present)
    if loglog:
        ax1.loglog(N_line, F_line, label="Design curve", **curve_kwargs)
    else:
        ax1.plot(N_line, F_line, label="Design curve", **curve_kwargs)

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

    ax2.set_axisbelow(True)  # grid below points/lines
    ax2.scatter(N_exp, N_pred, label="Predictions")

    # Determine parity plot range
    xy_low = min(np.min(N_exp), np.min(N_pred)) / 10.0
    xy_high = max(np.max(N_exp), np.max(N_pred)) * 10.0
    xy_low = float(max(xy_low, 1e-300))
    xy_high = float(max(xy_high, xy_low * 1.0001))

    if loglog:
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        line = np.logspace(np.log10(xy_low), np.log10(xy_high), 300)
    else:
        line = np.linspace(xy_low, xy_high, 300)

    # Black diagonal (y=x)
    ax2.plot(line, line, color="black", linestyle="--", linewidth=1.5, label="y = x")

    # Thin dashed error-factor bands (2 and 5 by default)
    for f in error_factors:
        f = float(f)
        ax2.plot(line, f * line, linestyle="--", linewidth=0.8, label=f"×{f:g}")
        ax2.plot(line, line / f, linestyle="--", linewidth=0.8, label=f"÷{f:g}")

    ax2.set_title("Parity plot: N experimental vs N expected")
    ax2.set_xlabel("N experimental (cycles)")
    ax2.set_ylabel("N expected (cycles)")

    ax2.grid(True, which="both")

    # Equal axis factor (make x/y ranges match, then use equal aspect)
    ax2.set_xlim(xy_low, xy_high)
    ax2.set_ylim(xy_low, xy_high)
    ax2.set_aspect("equal", adjustable="box")

    ax2.legend()

    if show:
        plt.show()

    if return_figax:
        return fig1, ax1, fig2, ax2
    return None

