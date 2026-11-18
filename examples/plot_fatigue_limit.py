"""Biaxial Fatigue Limit Envelopes Plot.

Failure boundaries for fully reversed, proportional (in-phase) axial + torsional loading in
normalized stress space (sigma_a / sigma_{-1} vs tau_a / tau_{-1}), calibrated on the 42CrMo4_Santus
fatigue limits of ``Experimental_Data/6_42CrMo4_Santus.xlsx`` (Benedetti et al. 2022):

- Gough-Pollard ellipse
- Maximum shear stress (Tresca), calibrated on sigma_{-1}
- Crossland (J2 amplitude + peak hydrostatic stress), calibrated on sigma_{-1} and tau_{-1}
- Findley critical plane, calibrated on sigma_{-1} and tau_{-1}

The dataset holds no combined-loading fatigue-limit tests, so only the two calibration points
are drawn; no combined-loading data are shown.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Ensure repository root is on sys.path for direct script execution
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import DOCS_IMG_DIR, INK, SERIES, apply_homepage_style, homepage_title

SIGMA_M1 = 433.0  # MPa, tensile fatigue limit at R = -1 (workbook Summary, Benedetti 2022)
TAU_M1 = 344.0  # MPa, torsional fatigue limit at R = -1 (same source)


def tresca_envelope(r: float, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """sqrt((sigma/2)^2 + tau^2) = sigma_{-1}/2  ->  x^2 + (2 r y)^2 = 1."""
    x = np.linspace(0.0, 1.0, n)
    return x, np.sqrt(1.0 - x**2) / (2.0 * r)


def crossland_envelope(r: float, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """sqrt(J2a) + a * sigma_H,max = b with sigma_H,max = sigma_a / 3 and b = tau_{-1}.

    Calibrating on sigma_{-1} and tau_{-1} gives a = 3 (r - 1/sqrt(3)).
    """
    a = 3.0 * (r - 1.0 / np.sqrt(3.0))
    x = np.linspace(0.0, 1.0, n)
    rhs = (1.0 - a * x / (3.0 * r)) ** 2 - x**2 / (3.0 * r**2)
    return x, np.sqrt(np.maximum(rhs, 0.0))


def findley_envelope(r: float, n_dir: int = 181) -> tuple[np.ndarray, np.ndarray, float]:
    """max over planes of (tau_a + k sigma_n,max) = tau_{-1} sqrt(1 + k^2).

    Stress is expressed in units of tau_{-1}; the axial amplitude is x / r.
    """
    k = (2.0 * r - 1.0) / (2.0 * np.sqrt(r * (1.0 - r)))
    th = np.linspace(0.0, np.pi, 1500)
    xs, ys = [], []
    for phi in np.linspace(0.0, np.pi / 2.0, n_dir):
        ux, uy = np.cos(phi), np.sin(phi)
        sig = 0.5 * (ux / r) * (1.0 + np.cos(2.0 * th)) + uy * np.sin(2.0 * th)
        tau = np.abs(-0.5 * (ux / r) * np.sin(2.0 * th) + uy * np.cos(2.0 * th))
        rho = np.sqrt(1.0 + k**2) / np.max(tau + k * np.maximum(0.0, sig))
        xs.append(rho * ux)
        ys.append(rho * uy)
    return np.array(xs), np.array(ys), k


def plot_fatigue_limit_envelopes(
    out_path: Path | None = None,
) -> Path:
    """Normalized biaxial fatigue limit envelopes of four criteria, homepage style."""
    apply_homepage_style()
    if out_path is None:
        out_path = DOCS_IMG_DIR / "fatigue_limit_envelopes.png"
    r = TAU_M1 / SIGMA_M1
    x = np.linspace(0.0, 1.0, 400)
    x_tr, y_tr = tresca_envelope(r)
    x_cr, y_cr = crossland_envelope(r)
    x_fi, y_fi, k_fi = findley_envelope(r)
    # (x, y, legend, direct label, colour, dash, anchor x on the curve, label position)
    curves = [
        (x, np.sqrt(1.0 - x**2), "Gough-Pollard ellipse", "Gough-Pollard", SERIES[0], "-",
         0.62, (0.80, 0.98)),
        (x_cr, y_cr, "Crossland", "Crossland", SERIES[1], "--", 0.80, (0.98, 0.80)),
        (x_fi, y_fi, f"Findley critical plane ($k = {k_fi:.2f}$)", "Findley", SERIES[3], ":",
         0.40, (0.16, 0.70)),
        (x_tr, y_tr, r"Maximum shear stress (Tresca, on $\sigma_{-1}$)", "Tresca", SERIES[2],
         "-.", 0.45, (0.20, 0.44)),
    ]  # fmt: skip

    fig, ax = plt.subplots(figsize=(8.6, 6.6), dpi=300)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.87, bottom=0.10)
    for xs, ys, legend, label, colour, dash, at, where in curves:
        ax.plot(xs, ys, color=colour, linewidth=2.4, linestyle=dash, label=legend,
                solid_capstyle="round", dash_capstyle="round")  # fmt: skip
        k = int(np.argmin(np.abs(xs - at)))
        ax.annotate(label, (xs[k], ys[k]), xytext=where, textcoords="data", ha="center",
                    va="center", fontsize=10, fontweight="semibold", color=INK["primary"],
                    arrowprops=dict(arrowstyle="-", color=colour, lw=1.0,
                                    shrinkA=4, shrinkB=2))  # fmt: skip
    ax.scatter([1.0, 0.0], [0.0, 1.0], s=70, color=INK["primary"], edgecolors="white",
               linewidths=1.6, zorder=5, clip_on=False,
               label=rf"Fatigue limits, $\sigma_{{-1}}$ = {SIGMA_M1:.0f}, $\tau_{{-1}}$ = {TAU_M1:.0f} MPa")  # fmt: skip
    ax.annotate(r"$\sigma_{-1}$", (1.0, 0.0), xytext=(8, 8), textcoords="offset points",
                fontsize=11, color=INK["primary"])  # fmt: skip
    ax.annotate(r"$\tau_{-1}$", (0.0, 1.0), xytext=(8, 6), textcoords="offset points",
                fontsize=11, color=INK["primary"])  # fmt: skip
    ax.set_xlim(0.0, 1.12)
    ax.set_ylim(0.0, 1.12)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Axial amplitude  $\sigma_a / \sigma_{-1}$", labelpad=6)
    ax.set_ylabel(r"Shear amplitude  $\tau_a / \tau_{-1}$", labelpad=6)
    ax.tick_params(length=0)
    ax.legend(loc="lower left", fontsize=8.8, handlelength=2.6)
    homepage_title(
        fig,
        "Where the criteria disagree",
        "Fatigue limit envelopes for in-phase tension-torsion, R = -1, calibrated on the "
        + rf"42CrMo4_Santus limits ($\tau_{{-1}}/\sigma_{{-1}} = {r:.2f}$)",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    p = plot_fatigue_limit_envelopes()
    print(f"Generated fatigue limit plot: {p}")
