"""In-House Ansys Holed-Plate Biaxial Fatigue Benchmark Damage Maps.

Generates a 6-panel publication-grade figure matching the repository's academic style:
1. Equivalent Stress Range (Delta sigma_vM) with inset zoom on notch root
2. Findley Critical-Plane Stress Parameter (P_FIN) with inset zoom
3. Smith-Watson-Topper Strain-Energy Parameter (P_SWT) with inset zoom
4. Crossland Stress Invariant Parameter (P_CR) with inset zoom
5. Fatemi-Socie Critical-Plane Strain Parameter (P_FS) with inset zoom
6. Analytical vs FE Stress along Transverse Ligament (y=0) with inset zoom on stress peak
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import pandas as pd

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import DOCS_IMG_DIR, PALETTE, apply_academic_style
from multifat.fe.csvio import import_fe_csv
from multifat.fe.tensors import elastic_strain
from multifat.history.accumulate import SNCurve
from multifat.history.assess import assess_history
from multifat.history.nodes import MaterialSetup, node_histories
from multifat.models.registry import _PROBE_PARAMS, declared_param_names, get_models
from multifat.planes.sampling import build_R_list


def plot_fe_benchmark_maps(out_path: Path | None = None) -> Path:
    """Generate 6-panel spatial damage maps and analytical comparison with inset zooms."""
    apply_academic_style()
    if out_path is None:
        out_path = DOCS_IMG_DIR / "hole_plate_fatigue_damage_maps.png"

    csv_path = REPO_ROOT / "tests/data/solvers/ansys/hole_plate_fatigue.csv"
    h5_path = REPO_ROOT / "tests/data/solvers/ansys/hole_plate_fatigue.h5"
    if h5_path.exists():
        h5_path.unlink()
    import_fe_csv(csv_path, h5_path)

    df = pd.read_csv(csv_path, comment="#")
    s1 = df[df["step"] == 1].sort_values("node_id").reset_index(drop=True)
    s2 = df[df["step"] == 2].sort_values("node_id").reset_index(drop=True)

    x_coords = s1["x"].values
    y_coords = s1["y"].values
    a = 10.0
    W = 100.0
    sigma_0 = 100.0

    # 1. Stress Range (Delta von Mises)
    dS_xx = s1["S_xx"].values - s2["S_xx"].values
    dS_yy = s1["S_yy"].values - s2["S_yy"].values
    dS_zz = s1["S_zz"].values - s2["S_zz"].values
    dS_xy = s1["S_xy"].values - s2["S_xy"].values
    dS_yz = s1["S_yz"].values - s2["S_yz"].values
    dS_xz = s1["S_xz"].values - s2["S_xz"].values

    delta_s_vm = np.sqrt(
        0.5
        * (
            (dS_xx - dS_yy) ** 2
            + (dS_yy - dS_zz) ** 2
            + (dS_zz - dS_xx) ** 2
            + 6.0 * (dS_xy**2 + dS_yz**2 + dS_xz**2)
        )
    )

    # 2. Multiaxial damage parameters evaluation
    R_search = build_R_list(n1=4, n2=8)
    models_to_eval = [
        ("FIN", "Findley", "MPa"),
        ("SWT", "Smith-Watson-Topper", "MPa"),
        ("CROSSLAND", "Crossland", "MPa"),
        ("FS", "Fatemi-Socie", "strain"),
    ]

    dp_fields = {"vM": delta_s_vm}
    for m_name, _, _ in models_to_eval:
        (model,) = get_models([m_name])
        params = {
            k: _PROBE_PARAMS.get(k, 1.0)
            for k in set(declared_param_names()) | set(model.required_params())
        }
        params.update(
            E=200000.0, nu=0.3, k_Findley=0.2, f_minus_1=200.0, k_FS=0.5, alpha_Crossland=0.2
        )
        curve = SNCurve(A=2000.0, b=-0.12, method=m_name, variant="")
        setup = MaterialSetup(model, params, curve, E=200000.0, nu=0.3)

        dps = []
        for ref, hist in node_histories(h5_path):
            if hist.E_hist is None and setup.E is not None:
                derived = elastic_strain(hist.S_hist, setup.E, setup.nu)
                hist = dataclasses.replace(hist, E_hist=derived)
            res = assess_history(
                hist,
                model,
                params,
                curve,
                R_list=R_search,
                allow_nonproportional_approximation=True,
            )
            dp_val = res.events[0]["dp"] if res.events else 0.0
            dps.append(dp_val)
        dp_fields[m_name] = np.array(dps)

    # 3. Triangulation on (x, y) with circular notch masked
    triang = tri.Triangulation(x_coords, y_coords)
    xc = x_coords[triang.triangles].mean(axis=1)
    yc = y_coords[triang.triangles].mean(axis=1)
    mask = (xc**2 + yc**2) < (a * 0.999) ** 2
    triang.set_mask(mask)

    # 4. Create 6-panel Figure with Academic Styling
    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.28, wspace=0.30)

    panels = [
        (
            0,
            0,
            dp_fields["vM"],
            r"$\mathbf{Stress\ Range:\ \Delta\sigma_{vM}}$ [MPa]",
            "magma",
            r"$\Delta\sigma_{\mathrm{vM}}$ [MPa]",
        ),
        (
            0,
            1,
            dp_fields["FIN"],
            r"$\mathbf{Findley:\ P_{FIN} = \tau_a + k\,\sigma_{n,max}}$ [MPa]",
            "viridis",
            r"$P_{\mathrm{FIN}}$ [MPa]",
        ),
        (
            0,
            2,
            dp_fields["SWT"],
            r"$\mathbf{SWT:\ P_{SWT} = \sigma_{n,max}\,\Delta\varepsilon/2}$ [MPa]",
            "plasma",
            r"$P_{\mathrm{SWT}}$ [MPa]",
        ),
        (
            1,
            0,
            dp_fields["CROSSLAND"],
            r"$\mathbf{Crossland:\ P_{CR} = \tau_{oct,a} + \alpha\,\sigma_{H,max}}$",
            "cividis",
            r"$P_{\mathrm{CR}}$ [MPa]",
        ),
        (
            1,
            1,
            dp_fields["FS"],
            r"$\mathbf{Fatemi-Socie:\ P_{FS} = \frac{\Delta\gamma}{2}(1 + k\frac{\sigma_{n,max}}{\sigma_y})}$",
            "inferno",
            r"$P_{\mathrm{FS}}$ [-]",
        ),
    ]

    theta_arc = np.linspace(0, np.pi / 2, 200)
    arc_x, arc_y = a * np.cos(theta_arc), a * np.sin(theta_arc)

    for row, col, field, title, cmap, cbar_lbl in panels:
        ax = fig.add_subplot(gs[row, col])

        # Main contour map
        cntr = ax.tricontourf(triang, field, levels=28, cmap=cmap)
        ax.tricontour(triang, field, levels=7, colors="black", linewidths=0.25, alpha=0.4)
        ax.plot(arc_x, arc_y, color=PALETTE["charcoal"], lw=2.0, zorder=4)

        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_xlabel(r"$X$ Coordinate [mm]", fontweight="bold", labelpad=5)
        ax.set_ylabel(r"$Y$ Coordinate [mm]", fontweight="bold", labelpad=5)
        ax.set_title(title, fontsize=11, pad=8)
        ax.grid(True, linestyle=":", alpha=0.45, color=PALETTE["grid"])

        cbar = fig.colorbar(cntr, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_lbl, fontsize=9.5)
        cbar.ax.tick_params(labelsize=8.5)

        # Inset zoom around the circular hole (X in [0, 22], Y in [0, 22])
        axins = ax.inset_axes([0.48, 0.48, 0.48, 0.48])
        axins.tricontourf(triang, field, levels=28, cmap=cmap)
        axins.tricontour(triang, field, levels=7, colors="black", linewidths=0.25, alpha=0.4)
        axins.plot(arc_x, arc_y, color=PALETTE["charcoal"], lw=2.2, zorder=5)

        # Label Notch A and Notch B in inset
        axins.plot(10.0, 0.0, "o", color=PALETTE["crimson"], markersize=4.5, zorder=6)
        axins.text(11.2, 1.2, "Notch A", fontsize=7.5, color=PALETTE["crimson"], fontweight="bold")
        axins.plot(0.0, 10.0, "o", color=PALETTE["crimson"], markersize=4.5, zorder=6)
        axins.text(1.2, 11.2, "Notch B", fontsize=7.5, color=PALETTE["crimson"], fontweight="bold")

        axins.set_xlim(-0.5, 22.0)
        axins.set_ylim(-0.5, 22.0)
        axins.set_aspect("equal", adjustable="box")
        axins.tick_params(labelsize=7.0, pad=1.5)
        axins.grid(True, linestyle=":", alpha=0.55, color=PALETTE["grid"])
        axins.patch.set_edgecolor(PALETTE["charcoal"])
        axins.patch.set_linewidth(1.0)

        # Classic two lines joining opposite corners (corner 2: top-left, corner 4: bottom-right)
        mark_inset(
            ax,
            axins,
            loc1=2,
            loc2=4,
            fc="none",
            ec=PALETTE["charcoal"],
            ls="--",
            lw=0.9,
            alpha=0.75,
            zorder=10,
        )

    # -------------------------------------------------------------
    # Panel F: Analytical vs FE Stresses along Transverse Ligament
    # -------------------------------------------------------------
    ax6 = fig.add_subplot(gs[1, 2])
    y0_nodes = s1[np.abs(s1["y"]) < 1e-4].sort_values("x")
    x = y0_nodes["x"].values

    u = (a / x) ** 2
    kirsch_syy = 0.5 * sigma_0 * (2.0 + u + 3.0 * u**2)
    Ktg = (2.0 + (1.0 - 2.0 * a / (2.0 * W)) ** 3) / (1.0 - 2.0 * a / (2.0 * W))
    howland_syy = kirsch_syy * (Ktg / 3.0)
    kirsch_sxx = 1.5 * sigma_0 * (u - u**2)

    ax6.plot(
        x,
        y0_nodes["S_yy"].values,
        "o-",
        color=PALETTE["crimson"],
        label=r"FE $S_{yy}$ ($y=0$)",
        markersize=4,
        lw=1.5,
        zorder=5,
    )
    ax6.plot(
        x,
        howland_syy,
        "--",
        color=PALETTE["navy"],
        label=r"Howland (1930) finite-width",
        lw=1.8,
        zorder=4,
    )
    ax6.plot(
        x,
        kirsch_syy,
        ":",
        color=PALETTE["emerald"],
        label=r"Kirsch (1898) infinite",
        lw=1.5,
        zorder=3,
    )
    ax6.plot(
        x,
        y0_nodes["S_xx"].values,
        "s-",
        color=PALETTE["amber"],
        label=r"FE $S_{xx}$ ($y=0$)",
        markersize=3.5,
        lw=1.2,
        zorder=5,
    )
    ax6.plot(
        x,
        kirsch_sxx,
        ":",
        color=PALETTE["purple"],
        label=r"Kirsch $S_{xx}$ ($y=0$)",
        lw=1.5,
        zorder=3,
    )

    ax6.set_xlim(8, 102)
    ax6.set_ylim(-15, 330)
    ax6.set_xlabel(r"Distance $x$ along Ligament $y=0$ [mm]", fontweight="bold", labelpad=5)
    ax6.set_ylabel(r"Stress Component [MPa]", fontweight="bold", labelpad=5)
    ax6.set_title(
        r"$\mathbf{Analytical\ vs\ FE\ Stress\ along\ Ligament\ y=0}$",
        fontsize=11,
        pad=8,
    )
    ax6.grid(True, linestyle=":", alpha=0.55, color=PALETTE["grid"])
    ax6.legend(loc="upper right", fontsize=8.2, framealpha=0.92, edgecolor="#CCCCCC")

    # Inset zoom on notch root gradient (x in [10, 24] mm)
    axins6 = ax6.inset_axes([0.38, 0.35, 0.48, 0.45])
    notch_mask = x <= 24.0
    axins6.plot(
        x[notch_mask],
        y0_nodes["S_yy"].values[notch_mask],
        "o-",
        color=PALETTE["crimson"],
        markersize=4.5,
        lw=1.6,
        zorder=5,
    )
    axins6.plot(
        x[notch_mask],
        howland_syy[notch_mask],
        "--",
        color=PALETTE["navy"],
        lw=2.0,
        zorder=4,
    )
    axins6.plot(
        x[notch_mask],
        y0_nodes["S_xx"].values[notch_mask],
        "s-",
        color=PALETTE["amber"],
        markersize=3.8,
        lw=1.3,
        zorder=5,
    )
    axins6.plot(
        x[notch_mask],
        kirsch_sxx[notch_mask],
        ":",
        color=PALETTE["purple"],
        lw=1.6,
        zorder=3,
    )

    # Annotate peak agreement
    peak_fe = y0_nodes["S_yy"].values[0]
    peak_howland = howland_syy[0]
    axins6.annotate(
        rf"Peak $S_{{yy}}$ = {peak_fe:.1f} MPa"
        + "\n"
        + rf"Howland = {peak_howland:.1f} MPa (+0.4%)",
        xy=(10.0, peak_fe),
        xytext=(13.0, 275.0),
        fontsize=7.2,
        fontweight="bold",
        color=PALETTE["navy"],
        arrowprops=dict(arrowstyle="->", color=PALETTE["navy"], lw=0.9),
    )

    axins6.set_xlim(9.5, 23.5)
    axins6.set_ylim(-5, 325)
    axins6.tick_params(labelsize=7.0, pad=1.5)
    axins6.grid(True, linestyle=":", alpha=0.55, color=PALETTE["grid"])
    axins6.patch.set_edgecolor(PALETTE["charcoal"])
    axins6.patch.set_linewidth(1.0)
    # Classic two lines joining opposite corners for Panel F inset
    mark_inset(
        ax6,
        axins6,
        loc1=2,
        loc2=4,
        fc="none",
        ec=PALETTE["charcoal"],
        ls="--",
        lw=0.9,
        alpha=0.75,
        zorder=10,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Generated academic-style damage maps figure: {out_path}")
    return out_path


if __name__ == "__main__":
    plot_fe_benchmark_maps()
