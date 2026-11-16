"""Homepage figure: what happens when the constants come from another batch of the same steel.

Data: ``examples/data/batch_comparison/42CrMo4/transfer_penalty.csv`` (issue #71), one row per
method and design -> target pair, with the signed load-equivalent error of the target's assessed
tests predicted with the target's own constants (``bias_own``) and with the design batch's
(``bias_transfer``). Each marker is the median over the methods, the bar their interquartile
range. e_L = ln(1/s*) > 0 is conservative (the model reaches its curve at a lower load than the
test); the load is misjudged by exp(e_L) - 1. Regenerate with ``examples/batch_comparison.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import (  # noqa: E402
    DOCS_IMG_DIR,
    INK,
    SERIES,
    apply_homepage_style,
    homepage_title,
)

DATA = REPO_ROOT / "examples" / "data" / "batch_comparison" / "42CrMo4" / "transfer_penalty.csv"
BATCHES = ("42CrMo4_Santus", "42CrMo4_Reis", "42CrMo4_FABEST")  # weaker lot first
LABEL = {b: b.split("_")[1] for b in BATCHES}
COLOUR = dict(zip(BATCHES, SERIES))
MARKER = dict(zip(BATCHES, ("o", "s", "^")))


def transfer_summary(path: Path = DATA) -> pd.DataFrame:
    """Median and quartiles over the methods of the signed load error, per design -> target."""
    df = pd.read_csv(path)
    rows = []
    for target in BATCHES:
        sub = df[df["Target"] == target]
        own = sub.groupby("Method")["bias_own"].mean()  # the same tests, whatever the design
        rows.append(dict(target=target, design=target, values=own.to_numpy()))
        for design in BATCHES:
            pair = sub[sub["Design"] == design]
            if design != target and len(pair):
                rows.append(
                    dict(target=target, design=design, values=pair["bias_transfer"].to_numpy())
                )
    out = pd.DataFrame(rows)
    out["n_methods"] = out["values"].map(len)
    out["median"] = out["values"].map(np.median)
    out["q1"] = out["values"].map(lambda v: np.quantile(v, 0.25))
    out["q3"] = out["values"].map(lambda v: np.quantile(v, 0.75))
    return out.drop(columns="values")


def plot_batch_transfer(out_path: Path = DOCS_IMG_DIR / "batch_transfer.png") -> Path:
    import matplotlib.pyplot as plt

    apply_homepage_style()
    summary = transfer_summary()
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.72, bottom=0.14)
    for row, target in enumerate(BATCHES):
        y0 = len(BATCHES) - 1 - row
        others = [b for b in BATCHES if b != target]
        offsets = {target: 0.24, others[0]: 0.0, others[1]: -0.24}  # own constants on top
        for _, r in summary[summary["target"] == target].iterrows():
            y = y0 + offsets[r["design"]]
            own = r["design"] == target
            colour = INK["muted"] if own else COLOUR[r["design"]]
            ax.plot([r["q1"], r["q3"]], [y, y], color=colour, lw=2.2, solid_capstyle="round")
            ax.plot(r["median"], y, MARKER[r["design"]], ms=8.5, color=colour,
                    markerfacecolor="white" if own else colour, markeredgewidth=1.6)  # fmt: skip
            text = "own constants" if own else f"from {LABEL[r['design']]}"
            pct = 100 * (np.exp(r["median"]) - 1)
            ax.annotate(f"{text}  {pct:+.0f} %", (r["q3"], y), xytext=(6, 0),
                        textcoords="offset points", va="center", fontsize=8.8,
                        color=INK["secondary"] if own else colour,
                        bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))  # fmt: skip
    ax.axvline(0.0, color=INK["secondary"], lw=1.0)
    ax.set_yticks(range(len(BATCHES)))
    ax.set_yticklabels([f"{LABEL[b]} tests" for b in reversed(BATCHES)], fontsize=10.5)
    ax.set_ylim(-0.6, len(BATCHES) - 0.4)
    ax.grid(axis="y", visible=False)
    lo, hi = summary["q1"].min(), summary["q3"].max()
    ax.set_xlim(min(-0.3, lo - 0.05), max(0.45, hi + 0.22))
    ax.set_xlabel("load-equivalent error $e_L$ = ln(1/s*), median over methods (bar: quartiles)")
    ax.text(0.99, 1.02, "conservative →", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color=INK["secondary"])  # fmt: skip
    ax.text(0.01, 1.02, "← non-conservative", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, color=INK["secondary"])  # fmt: skip
    n = int(summary["n_methods"].max())
    homepage_title(
        fig,
        "Constants from another batch of the same steel",
        f"42CrMo4, three laboratories; up to {n} methods, each calibrated on one batch and used on "
        "another.\nPercentages: median misjudgement of the load.",
    )
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(plot_batch_transfer())
