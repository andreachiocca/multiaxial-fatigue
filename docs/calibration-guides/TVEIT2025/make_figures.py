"""Figures of the TVEIT2025 calibration guide, drawn from the calibrations of the repository.

Run from the repository root:

    python docs/calibration-guides/TVEIT2025/make_figures.py

- ``fig_haigh_routes.pdf``: schematic Haigh diagram of the routes to A (no data).
- ``fig_example_<name>.pdf``: one panel per workbook of a worked example. It shows the closed-form
  uniaxial R = -1 curves of the equivalent sets of the profile over C (E <= 1.05 E_min), the
  chosen set, the slowest and the fastest C among them, the axial R = -1 failures and run-outs of
  the workbook and the three calibration points.

The calibration is the one the pipeline uses: the committed record in ``Calibrations/<material>/``
is reused, and a workbook without one (Al7075_T651_Wang) is calibrated, in about 20 s.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples"))
from style import INK, SERIES  # noqa: E402

from multifat.eval import tveit2025_calibration as tc  # noqa: E402
from multifat.imp.experimental_data import load_material_workbook  # noqa: E402

DATA = ROOT / "Experimental_Data"
OUT = HERE / "figures"
CHOSEN, SLOW, FAST = SERIES[0], SERIES[1], SERIES[2]
FAMILY = "#c9c8c2"
PDF_METADATA = {"CreationDate": None, "Creator": None, "Producer": None}  # same bytes on rerun

# figure name -> (workbook, panel label) of each panel
EXAMPLES = {
    "Al7075T6": [("Al7075_T6.xlsx", "Al7075-T6 (Pedranz et al.)")],
    "Ti6Al4V": [
        ("15_Ti6Al4V_Berto.xlsx", "Ti6Al4V_Berto"),
        ("1_Ti6Al4V_Meneghetti.xlsx", "Ti6Al4V_Meneghetti"),
    ],
    "42CrMo4": [
        ("6_42CrMo4_Santus.xlsx", "42CrMo4_Santus (limits of Benedetti et al.)"),
        ("14_42CrMo4_Reis.xlsx", "42CrMo4_Reis (Anes et al.)"),
    ],
    "Al7075T651": [("18_Al7075_T651_Wang.xlsx", "Al7075-T651 (Wang et al.)")],
}

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": INK["secondary"],
        "axes.labelcolor": INK["secondary"],
        "xtick.color": INK["secondary"],
        "ytick.color": INK["secondary"],
        "mathtext.fontset": "dejavusans",
    }
)


def calibration(workbook: Path):
    """The pipeline's calibration of a workbook, with its imported data."""
    data = load_material_workbook(workbook)
    material = next(iter(data))
    result = tc.calibrate_tveit2025(workbook, data, material)
    if result["status"] != "calibrated":
        raise RuntimeError(f"{workbook.name}: {result['status']}: {result['reason']}")
    return data[material], result


def tests(workbook: Path, block, record):
    """Failures of the calibration source, failures of other sources, run-outs: (sigma_a, N)."""
    sources = tc.dataset_sources(workbook)
    failures = tc.axial_tests(block["cases"], sources, failures=True)
    runouts = tc.axial_tests(block.get("runouts", {}), sources, failures=False)
    declared = (record.get("inputs") or {}).get("calibration_source")
    pool, _ = tc.calibration_source(failures, record["fatigue_limits"], declared)
    ids = {t["case_id"] for t in pool}
    other = [t for t in failures if t["case_id"] not in ids]
    return [[(t["sigma_a"], t["N_exp"]) for t in group] for group in (pool, other, runouts)]


def curve(s, A, f_m1, smax):
    """Closed-form uniaxial R = -1 S-N curve of the set s (dict with C, K, L)."""
    amplitudes = np.geomspace(f_m1 * 1.0005, smax, 70)
    lives = [tc.uniaxial_life(s["C"], s["K"], s["L"], 0.0, a, A, f_m1).Nf for a in amplitudes]
    lives = np.array(lives, dtype=float)
    ok = np.isfinite(lives) & (lives > 1.0)
    return lives[ok], amplitudes[ok]


def panel(ax, workbook: Path, label: str) -> str:
    """Draw one workbook; return a one-line summary for the console."""
    block, result = calibration(workbook)
    record = result["record"]
    f_m1, A = record["fatigue_limits"]["Sigm1"], record["mean_stress_slope"]["A"]
    family = [p for p in record["profile"] if p["E"] <= tc.EQUIVALENT * record["E_min"]]
    if len(family) != record["family"]["n_equivalent"]:
        raise RuntimeError(f"{workbook.name}: the equivalent family does not match the record")
    chosen = record["chosen"]
    slow, fast = min(family, key=lambda p: p["C"]), max(family, key=lambda p: p["C"])
    pool, other, runouts = tests(workbook, block, record)
    stresses = [s for s, _ in pool + other + runouts] + [f_m1]
    smin, smax = min(stresses) * 0.96, max(stresses) * 1.06

    for s in family:
        ax.plot(*curve(s, A, f_m1, smax), color=FAMILY, lw=0.8, zorder=1)
    # The curves coincide: the chosen set underneath, the slowest and fastest C on top with
    # alternating dashes, so that all three stay visible.
    ax.plot(*curve(chosen, A, f_m1, smax), color=CHOSEN, lw=3.2, zorder=3)
    for s, colour, dashes in ((slow, SLOW, (0, (4, 4))), (fast, FAST, (4, (4, 4)))):
        ax.plot(*curve(s, A, f_m1, smax), color=colour, lw=1.5, ls=dashes, zorder=3.5)
    if other:
        ax.scatter(*zip(*[(n, s) for s, n in other]), s=22, color=INK["muted"], lw=0, zorder=4)
    if runouts:
        ax.scatter(
            *zip(*[(n, s) for s, n in runouts]),
            s=34,
            marker=">",
            facecolors="white",
            edgecolors=INK["secondary"],
            lw=1.0,
            zorder=4,
        )
    ax.scatter(
        *zip(*[(n, s) for s, n in pool]),
        s=30,
        facecolors="white",
        edgecolors=INK["primary"],
        lw=1.0,
        zorder=4,
    )
    for p in record["calibration_points"]:
        ax.scatter([p["N_exp"]], [p["sigma_a"]], s=55, marker="s", color=INK["primary"], zorder=5)
        ax.annotate(
            p["role"],
            (p["N_exp"], p["sigma_a"]),
            xytext=(5, 7),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            color=INK["primary"],
            zorder=6,
        )
    ax.axhline(f_m1, color=INK["secondary"], lw=1.0, ls=":", zorder=2)
    ax.text(1.15e3, f_m1, rf" $\sigma_{{-1}}$ = {f_m1:g} MPa", va="bottom", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(1e3, 1e8)
    ax.set_ylim(smin, smax)
    ax.set_xlabel("cycles to failure $N$")
    ax.set_ylabel(r"axial $\sigma_a$ [MPa], $R=-1$")
    level = result["confidence"]["level"]
    ax.set_title(
        f"{label}: confidence {level}\n"
        f"{len(family)} equivalent sets, C from {slow['C']:.2g} to {fast['C']:.2g}",
        fontsize=9,
        color=INK["primary"],
        loc="left",
    )
    lines = [
        f"{name}: C = {s['C']:.3g}, K = {s['K']:.3g}, L = {s['L']:.3g}"
        for name, s in (("chosen", chosen), ("slowest C", slow), ("fastest C", fast))
    ]
    ax.text(
        0.98,
        0.96,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        linespacing=1.5,
    )
    ax.grid(True, color=INK["grid"], lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return f"{workbook.name}: {len(family)} sets, C {slow['C']:.3g}-{fast['C']:.3g}, {level}"


def example(name: str, panels) -> None:
    """One figure, one panel per workbook, with a shared legend."""
    height = 3.9 * len(panels) + 0.8  # inches; the bottom 0.8 in holds the legend
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.3, height), squeeze=False)
    for ax, (workbook, label) in zip(axes[:, 0], panels):
        print(panel(ax, DATA / workbook, label))
    handles = [
        Line2D([], [], color=CHOSEN, lw=3, label="chosen set"),
        Line2D([], [], color=SLOW, lw=1.5, ls="--", label="slowest C"),
        Line2D([], [], color=FAST, lw=1.5, ls="--", label="fastest C"),
        Line2D([], [], color=FAMILY, lw=1.5, label="other equivalent sets"),
        Line2D([], [], color=INK["primary"], marker="s", ls="", ms=6, label="calibration points"),
        Line2D(
            [],
            [],
            color=INK["primary"],
            marker="o",
            mfc="white",
            ls="",
            ms=5,
            label="failures, calibration source",
        ),
        Line2D([], [], color=INK["muted"], marker="o", ls="", ms=4, label="failures, other source"),
        Line2D(
            [], [], color=INK["secondary"], marker=">", mfc="white", ls="", ms=5, label="run-outs"
        ),
    ]
    fig.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.8 / height, 1, 1))
    fig.savefig(OUT / f"fig_example_{name}.pdf", metadata=PDF_METADATA)
    plt.close(fig)


def haigh() -> None:
    """Schematic Haigh diagram: the routes to A differ only in how the slope is obtained."""
    f_m1, su = 1.0, 3.0
    sm = np.linspace(0, su * 1.02, 50)
    m0, a0 = 0.9, 0.58
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    routes = (
        ((f_m1 - a0) / m0, CHOSEN, "-", r"route 1, a limit at a mean stress: $A=(\sigma_{-1}-\sigma_{a,0})/\sigma_{m,0}$"),
        (0.2, FAST, "-.", r"route 3, FKM: $A=M_\sigma=a_M\,S_u/1000+b_M$"),
        (f_m1 / su, SLOW, "--", r"route 4, Goodman: $A=\sigma_{-1}/S_u$"),
    )  # fmt: skip
    for A, colour, style, label in routes:
        amplitude = f_m1 - A * sm
        ok = amplitude >= 0
        ax.plot(sm[ok], amplitude[ok], color=colour, lw=2, ls=style, label=label)
    ax.scatter([m0], [a0], s=45, color=INK["primary"], zorder=5)
    ax.annotate(
        "fatigue limit at a\n" + r"mean stress $(\sigma_{m,0},\ \sigma_{a,0})$",
        (m0, a0),
        xytext=(0.12, 0.28),
        fontsize=8,
        va="center",
        arrowprops=dict(arrowstyle="-", color=INK["secondary"], lw=0.8, shrinkB=5),
    )
    ax.scatter([0, su], [f_m1, 0], s=45, color=INK["primary"], zorder=5, clip_on=False)
    ax.annotate(r"$\sigma_{-1}$  ($R=-1$)", (0, f_m1), xytext=(8, 6), textcoords="offset points")
    ax.annotate(r"$S_u$", (su, 0), xytext=(-4, 8), textcoords="offset points", ha="right")
    ax.set_xlim(0, su * 1.05)
    ax.set_ylim(0, f_m1 * 1.18)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"mean stress $\sigma_m$")
    ax.set_ylabel(r"stress amplitude $\sigma_a$")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=7.8, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "fig_haigh_routes.pdf", metadata=PDF_METADATA)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams["pdf.fonttype"] = 42
    haigh()
    for name, panels in EXAMPLES.items():
        example(name, panels)
    print("figures written to", OUT.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()
