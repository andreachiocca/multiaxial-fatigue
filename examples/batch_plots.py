"""Figures of the cross-batch comparison (issue #71), one message per figure.

``all_figures`` draws, into ``docs/images/batch_comparison_<family>_<name>.png``:

- ``sn``: the R = -1 S-N data of the batches per loading path, with their Basquin fits;
- ``strength``: the strength S(1e5) of every batch and series per path with 95 % intervals, and
  the range between the same-steel series of the FABEST database as a yardstick;
- ``offsets``: the batch offsets of the method errors with joint-bootstrap intervals, against
  their method-free prediction;
- ``interaction``: bias of every method on one batch against another batch; a pure batch offset
  moves all methods along a line of unit slope;
- ``transfer_bias``: median bias of predicting one batch from another (conservative or not);
- ``designer``: the load factor that covers 95 % of the tests, own batch against transfer.

The batches are distinguished by marker as well as by colour; the tables behind every figure
are in ``examples/data/batch_comparison/<family>/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.batch_reference import N_REF, reversed_tests  # noqa: E402
from examples.style import (  # noqa: E402
    BATCH_COLOURS,
    BATCH_MARKERS,
    DIVERGING,
    DOCS_IMG_DIR,
    apply_academic_style,
)

PATHS = ["A", "T", "AT-IP", "AT-NP"]
PATH_TEXT = {
    "A": "axial (A)",
    "T": "torsion (T)",
    "AT-IP": "in-phase tension-torsion (AT-IP)",
    "AT-NP": "90° out-of-phase (AT-NP)",
}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dddddd"


def _short(batch: str, family: str) -> str:
    return batch.removeprefix(family + "_") or batch


def _style(batches) -> dict:
    return {b: (BATCH_COLOURS[k % 4], BATCH_MARKERS[k % 4]) for k, b in enumerate(batches)}


def _save(fig, family: str, name: str) -> Path:
    path = DOCS_IMG_DIR / f"batch_comparison_{family}_{name}.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def figure_sn(family: str, batches, series: pd.DataFrame) -> Path | None:
    """S-N data and pooled Basquin fits of the batches, one panel per loading path."""
    tests = pd.concat([reversed_tests(b) for b in batches], ignore_index=True)
    paths = [p for p in PATHS if tests[tests["Path"] == p]["Batch"].nunique() >= 2]
    if not paths:
        return None
    style = _style(batches)
    fig, axes = plt.subplots(1, len(paths), figsize=(3.1 * len(paths) + 0.6, 3.6), sharey=True)
    axes = np.atleast_1d(axes)
    pooled = series[series["Series"] == "all"]
    for ax, path in zip(axes, paths):
        for b in batches:
            g = tests[(tests["Batch"] == b) & (tests["Path"] == path)]
            if g.empty:
                continue
            colour, marker = style[b]
            ax.scatter(g["N"], g["amp"], s=16, marker=marker, facecolor="none",
                       edgecolor=colour, linewidth=0.9, label=_short(b, family))  # fmt: skip
            fit = pooled[(pooled["Batch"] == b) & (pooled["Path"] == path)]
            if len(fit):
                r = fit.iloc[0]
                N = np.array([r["N_min"], r["N_max"]])
                ax.plot(N, r["S_ref"] * (N / N_REF) ** r["b"], color=colour, lw=1.4)
        ax.axvline(N_REF, color=MUTED, lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(PATH_TEXT[path], fontsize=9.5)
        ax.set_xlabel(r"cycles to failure $N_f$")
        ax.grid(which="both", color=GRID, lw=0.4)
    ticks = [300, 400, 500, 600, 800]
    axes[0].set_yticks(ticks, [str(t) for t in ticks])
    axes[0].set_yticks([], minor=True)
    axes[0].set_ylabel(r"amplitude $\sigma_a$ ($\tau_a$ in torsion), MPa")
    handles, labels = axes[0].get_legend_handles_labels()
    extra = {}
    for ax in axes[1:]:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            extra.setdefault(lab, h)
    for lab, h in zip(labels, handles):
        extra[lab] = h
    order = [_short(b, family) for b in batches if _short(b, family) in extra]
    fig.legend([extra[k] for k in order], order, loc="upper center", ncol=len(order),
               frameon=False, bbox_to_anchor=(0.5, 0.93))  # fmt: skip
    fig.suptitle(f"{family}: fully reversed S-N data of the batches (line: Basquin fit; "
                 r"dotted: $N = 10^5$)", y=0.995, fontsize=11)  # fmt: skip
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return _save(fig, family, "sn")


def figure_strength(family: str, batches, series: pd.DataFrame, fabest=None) -> Path | None:
    """S(1e5) of every batch (filled) and series (open) per path, with the FABEST yardstick."""
    paths = [p for p in PATHS if (series["Path"] == p).any()]
    if not paths:
        return None
    style = _style(batches)
    fig, ax = plt.subplots(figsize=(7.4, 0.95 * len(paths) + 1.9))
    step = 0.8 / max(len(batches), 1)
    for i, path in enumerate(paths):
        if fabest is not None and len(fabest):
            load = {"A": "Axial", "T": "Torsion"}.get(path)
            band = fabest[(fabest["Load"] == load) & (fabest["group"] == "all unnotched series")]
            if len(band):
                r = band.iloc[0]
                ax.fill_betweenx([i - 0.46, i + 0.46], r["S_min"], r["S_max"], color="#e6e5e0",
                                 lw=0, zorder=0)  # fmt: skip
                ax.text(r["S_min"], i - 0.5, f"FABEST database: range of {int(r['series'])} "
                        f"series of one steel, {int(r['labs'])} laboratories",
                        ha="left", va="bottom", fontsize=6.5, color=MUTED)  # fmt: skip
        for k, b in enumerate(batches):
            g = series[(series["Batch"] == b) & (series["Path"] == path)]
            y = i + (k - (len(batches) - 1) / 2) * step
            colour, marker = style[b]
            for r in g.itertuples():
                pooled = r.Series == "all"
                ax.errorbar(
                    r.S_ref, y,
                    xerr=[[r.S_ref - r.S_lo], [min(r.S_hi, 2 * r.S_ref) - r.S_ref]],
                    fmt=marker, ms=6 if pooled else 4, color=colour,
                    mfc=colour if pooled else "white", lw=1.2 if pooled else 0.7,
                    alpha=1.0 if pooled else 0.9, zorder=3 if pooled else 2,
                    label=_short(b, family) if pooled and i == 0 else None,
                )  # fmt: skip
    ax.set_yticks(range(len(paths)), [PATH_TEXT[p] for p in paths])
    ax.invert_yaxis()
    ax.set_xlabel(r"strength at $10^5$ cycles, $S(10^5)$, MPa" "\n"
                  r"($\tau_a$ in torsion, $\sigma_a$ otherwise; bars: 95 % bootstrap interval)",
                  fontsize=9.5)  # fmt: skip
    ax.grid(axis="x", color=GRID, lw=0.5)
    handles, labels = ax.get_legend_handles_labels()
    seen = dict(zip(labels, handles))
    order = [_short(b, family) for b in batches if _short(b, family) in seen]
    ax.legend([seen[k] for k in order], order, loc="lower right", frameon=True, fontsize=8,
              title="filled: all tests; open: one series", title_fontsize=7)  # fmt: skip
    ax.set_title(f"{family}: strength of the batches per loading path", fontsize=11)
    fig.tight_layout()
    return _save(fig, family, "strength")


def figure_offsets(family: str, batches, offsets: pd.DataFrame) -> Path | None:
    """Batch offsets of the method errors (joint-bootstrap interval) and their prediction."""
    if offsets.empty:
        return None
    width = offsets["batches"].str.count(r"\+")
    largest = offsets[width == width.groupby(offsets["Path"]).transform("max")]
    paths = [p for p in PATHS + ["T+mean"] if (largest["Path"] == p).any()]
    style = _style(batches)
    fig, ax = plt.subplots(figsize=(7.2, 0.9 * len(paths) + 2.3))
    step = 0.7 / max(len(batches), 1)
    for i, path in enumerate(paths):
        g = largest[largest["Path"] == path]
        for k, b in enumerate(batches):
            r = g[g["Batch"] == b]
            if r.empty:
                continue
            r = r.iloc[0]
            y = i + (k - (len(batches) - 1) / 2) * step
            colour, marker = style[b]
            ax.errorbar(r["offset"], y, xerr=[[r["offset"] - r["offset_lo"]],
                        [r["offset_hi"] - r["offset"]]], fmt=marker, ms=6.5, color=colour,
                        lw=1.3, capsize=2.5, label=_short(b, family) if i == 0 else None)  # fmt: skip
            if np.isfinite(r.get("predicted", np.nan)):
                ax.plot(r["predicted"], y, marker="x", ms=7, mew=1.6, color=INK, ls="none",
                        zorder=4, label="from the S-N curves alone" if (i, k) == (0, 0) else None)  # fmt: skip
            ax.text(1.015, y, f"{r['offset']:+.3f}", va="center", fontsize=7.5, color=MUTED,
                    transform=ax.get_yaxis_transform())  # fmt: skip
    ax.axvline(0.0, color=MUTED, lw=0.8)
    ax.set_yticks(range(len(paths)), [PATH_TEXT.get(p, p) for p in paths])
    ax.invert_yaxis()
    ax.set_xlabel(r"batch offset $\beta_b$ of the load error $e_L=\ln(1/s^*)$" "\n"
                  "(mean over the methods; positive: more conservative on that batch;\n"
                  "bars: 95 % joint-bootstrap interval)", fontsize=9)  # fmt: skip
    ax.grid(axis="x", color=GRID, lw=0.5)
    handles, labels = ax.get_legend_handles_labels()
    seen = dict(zip(labels, handles))
    order = [_short(b, family) for b in batches if _short(b, family) in seen]
    order += [k for k in seen if k not in order]
    fig.legend([seen[k] for k in order], order, loc="upper center", ncol=len(order),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.94))  # fmt: skip
    fig.suptitle(f"{family}: batch offsets, each method calibrated on its own batch",
                 fontsize=11, y=0.995)  # fmt: skip
    fig.tight_layout(rect=(0, 0, 0.97, 0.9))
    return _save(fig, family, "offsets")


def figure_interaction(family: str, batches, cells: pd.DataFrame) -> Path | None:
    """Method biases on one batch against another: a pure offset is a line of unit slope."""
    panels = []
    for path in PATHS:
        g = cells[(cells["Path"] == path) & (cells["n"] >= 3)]
        table = g.pivot(index="Method", columns="Batch", values="bias")
        present = [b for b in batches if b in table]
        for a, b in zip(present[:-1], present[1:]):
            pair = table[[a, b]].dropna()
            if len(pair) >= 5:
                panels.append((path, a, b, pair))
    if not panels:
        return None
    fig, axes = plt.subplots(1, len(panels), figsize=(3.1 * len(panels) + 0.4, 3.7))
    axes = np.atleast_1d(axes)
    for ax, (path, a, b, pair) in zip(axes, panels):
        x, y = pair[a].to_numpy(), pair[b].to_numpy()
        shift = float(np.mean(y - x))
        lim = [min(x.min(), y.min()) - 0.05, max(x.max(), y.max()) + 0.05]
        ax.plot(lim, lim, color=MUTED, lw=0.7, ls="--", label="no batch effect")
        ax.plot(lim, [lim[0] + shift, lim[1] + shift], color=INK, lw=1.0,
                label=f"common offset {shift:+.2f}")  # fmt: skip
        ax.scatter(x, y, s=14, color=BATCH_COLOURS[0], edgecolor="white", linewidth=0.4, zorder=3)
        residual = y - x - shift
        for j in np.argsort(-np.abs(residual))[:3]:
            ax.annotate(pair.index[j], (x[j], y[j]), xytext=(4, 4), textcoords="offset points",
                        fontsize=6.5, color=INK)  # fmt: skip
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_title(f"{PATH_TEXT[path]}\n{len(pair)} methods; scatter about the line "
                     f"{residual.std(ddof=1):.3f}", fontsize=8.5)  # fmt: skip
        ax.set_xlabel(rf"bias $\bar e_L$ on {_short(a, family)}")
        ax.set_ylabel(rf"bias $\bar e_L$ on {_short(b, family)}")
        ax.grid(color=GRID, lw=0.4)
        ax.legend(fontsize=6.5, loc="upper left", frameon=True)
    fig.suptitle(f"{family}: every point is one method, calibrated on each batch in turn",
                 fontsize=11)  # fmt: skip
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return _save(fig, family, "interaction")


def figure_transfer_bias(family: str, batches, penalty: pd.DataFrame) -> Path | None:
    """Median bias over the methods, design batch (rows) -> assessed batch (columns)."""
    if penalty.empty:
        return None
    n = len(batches)
    grid = np.full((n, n), np.nan)
    text = [["" for _ in batches] for _ in batches]
    for i, design in enumerate(batches):
        for j, target in enumerate(batches):
            if i == j:
                g = penalty[penalty["Target"] == target].drop_duplicates("Method")
                values = g["bias_own"]
            else:
                g = penalty[(penalty["Design"] == design) & (penalty["Target"] == target)]
                values = g["bias_transfer"]
            if len(values) >= 3:
                q1, med, q3 = np.percentile(values, [25, 50, 75])
                grid[i, j] = med
                text[i][j] = f"{med:+.2f}\n[{q1:+.2f}, {q3:+.2f}]\n{len(values)} methods"
    if not np.isfinite(grid).any():
        return None
    span = max(0.05, float(np.nanmax(np.abs(grid))))
    cmap = LinearSegmentedColormap.from_list("bias", DIVERGING)
    fig, ax = plt.subplots(figsize=(1.75 * n + 2.4, 1.45 * n + 1.5))
    image = ax.imshow(grid, cmap=cmap, norm=TwoSlopeNorm(0.0, -span, span))
    for i in range(n):
        for j in range(n):
            if text[i][j]:
                strong = abs(grid[i, j]) > 0.6 * span
                ax.text(j, i, text[i][j], ha="center", va="center", fontsize=7.5,
                        color="white" if strong else INK,
                        fontweight="bold" if i == j else "normal")  # fmt: skip
    names = [_short(b, family) for b in batches]
    ax.set_xticks(range(n), names)
    ax.set_yticks(range(n), names)
    ax.set_xlabel("batch assessed (tests)")
    ax.set_ylabel("batch used for calibration (properties, constants, master curve)")
    for k in range(n + 1):
        ax.axhline(k - 0.5, color="white", lw=2)
        ax.axvline(k - 0.5, color="white", lw=2)
    bar = fig.colorbar(image, ax=ax, shrink=0.8)
    bar.set_label(r"median bias $\bar e_L$ over the methods" "\n"
                  r"$\leftarrow$ non-conservative  |  conservative $\rightarrow$")  # fmt: skip
    ax.set_title(f"{family}: predicting one batch from another\n(diagonal: own batch; "
                 "brackets: quartiles over the methods; same tests)", fontsize=10)  # fmt: skip
    fig.tight_layout()
    return _save(fig, family, "transfer_bias")


def figure_designer(family: str, batches, designer: pd.DataFrame) -> Path | None:
    """Load factor covering 95 % of the tests: own batch (grey) against transfer (colour)."""
    if designer.empty:
        return None
    pairs = [(d, t) for t in batches for d in batches if d != t]
    pairs = [
        p for p in pairs if ((designer["Design"] == p[0]) & (designer["Target"] == p[1])).any()
    ]
    style = _style(batches)
    fig, ax = plt.subplots(figsize=(1.25 * len(pairs) + 2.2, 4.0))
    for k, (design, target) in enumerate(pairs):
        g = designer[(designer["Design"] == design) & (designer["Target"] == target)]
        for offset, column, colour in ((-0.18, "gamma_own", "#9a9993"),
                                       (0.18, "gamma_transfer", style[design][0])):  # fmt: skip
            values = g[column].to_numpy(float)
            box = ax.boxplot(values, positions=[k + offset], widths=0.3, patch_artist=True,
                             showfliers=False, medianprops=dict(color=INK, lw=1.3))  # fmt: skip
            box["boxes"][0].set(facecolor=colour, alpha=0.55, edgecolor=colour)
            ax.text(k + offset, np.percentile(values, 75) + 0.012, f"{np.median(values):.2f}",
                    ha="center", fontsize=7, color=INK)  # fmt: skip
    ax.axhline(1.0, color=MUTED, lw=0.8)
    ax.set_xticks(range(len(pairs)),
                  [f"{_short(t, family)}\nfrom {_short(d, family)}" for d, t in pairs],
                  fontsize=8)  # fmt: skip
    ax.set_ylabel(r"load factor $\gamma_{95}=\exp[-(\bar e_L-1.645\,s)]$" "\n"
                  "needed to cover 95 % of the tests")  # fmt: skip
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.plot([], [], color="#9a9993", lw=6, alpha=0.55, label="calibrated on the batch itself")
    for b in batches:
        ax.plot([], [], color=style[b][0], lw=6, alpha=0.55,
                label=f"calibrated on {_short(b, family)}")  # fmt: skip
    fig.legend(*ax.get_legend_handles_labels(), fontsize=7.5, loc="lower center", ncol=2,
               frameon=False)  # fmt: skip
    ax.set_title(f"{family}: what a designer loses without tests on the batch\n"
                 "(box: quartiles over the methods; number: median)", fontsize=10)  # fmt: skip
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return _save(fig, family, "designer")


def all_figures(family, batches, *, reference, offsets, cells, penalty, designer) -> list[Path]:
    apply_academic_style()
    made = [
        figure_sn(family, batches, reference["series"]),
        figure_strength(family, batches, reference["series"], reference.get("fabest_summary")),
        figure_offsets(family, batches, offsets),
        figure_interaction(family, batches, cells),
        figure_transfer_bias(family, batches, penalty),
        figure_designer(family, batches, designer),
    ]
    return [p for p in made if p is not None]
