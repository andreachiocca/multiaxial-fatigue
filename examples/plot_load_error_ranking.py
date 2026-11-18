"""Ranking of every method by its load-equivalent error, across materials (issue #37).

Rows: methods, ordered by their mean rank over the materials they ran on (best on top). A material
on which fewer than MIN_RANKED methods are scored (Ti6Al4V_Meneghetti: only MOREL, issue #74) shows its cells
but is left out of the ranks, since a rank among one method carries no information.
Columns: materials. Cells: root mean square of Error_load = ln(1/s*) on the assessed failures,
which combines bias and scatter; lower is better. A cell of 0.10 means a typical load
misjudgement of about 10 %. The same scale holds for every
method, whatever its damage parameter, S-N slope or life law. Data: examples/data/load_error_ranking.csv
(see build_load_error_ranking.py).
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import (
    DATA_DIR, DOCS_IMG_DIR, INK, SURFACE, apply_homepage_style, blue_cmap, homepage_title,
)  # fmt: skip

MIN_RANKED = 2  # methods scored on a material for it to enter the ranks


GROUPS = (
    ("fatigue_life", "failures (life datasets)"),
    ("fatigue_limit", "limit states (fatigue-limit datasets)"),
)


def _ranks(sd: pd.DataFrame) -> pd.DataFrame:
    ranks = sd.rank(axis=0, method="average")
    ranks.loc[:, sd.notna().sum(axis=0) < MIN_RANKED] = np.nan  # no comparison on that material
    return ranks


def plot_load_error_ranking(
    out_path: Path | None = None, data_path: Path = DATA_DIR / "load_error_ranking.csv"
) -> Path:
    apply_homepage_style()
    out_path = out_path or DOCS_IMG_DIR / "load_error_ranking.png"
    table = pd.read_csv(data_path)
    if "Assessment" not in table:
        table["Assessment"] = "fatigue_life"
    # one block of columns per kind of observation, each with its own mean rank: a limit dataset
    # rests on 3-40 limit states, so its columns are not pooled with the life datasets
    # a material on which fewer than MIN_RANKED methods are scored compares nothing: left out of
    # the figure and named in a note (Ti6Al4V_Meneghetti: only MOREL, #74)
    counts = table.groupby("Material")["Method"].nunique()
    hidden = sorted(counts[counts < MIN_RANKED].index)
    table = table[~table["Material"].isin(hidden)]
    blocks, block_ranks = [], {}
    for kind, _ in GROUPS:
        part = table[table["Assessment"] == kind]
        if part.empty:
            continue
        sd_k = part.pivot(index="Method", columns="Material", values="Error_load_rms")
        blocks.append((kind, list(sd_k.columns)))
        block_ranks[kind] = _ranks(sd_k)
    sd = table.pivot(index="Method", columns="Material", values="Error_load_rms")
    materials = [m for _, cols in blocks for m in cols]
    first = blocks[0][0]
    order = block_ranks[first].mean(axis=1).reindex(sd.index).sort_values(na_position="last").index
    sd = sd.loc[order, materials]

    width, height = 1.25 * len(materials) + 3.6, 0.30 * len(sd) + 3.0
    fig, ax = plt.subplots(figsize=(width, height), dpi=300)
    fig.subplots_adjust(left=2.9 / width, right=1 - 1.3 / width, top=1 - 2.3 / height,
                        bottom=0.7 / height)  # fmt: skip
    ax.grid(False)
    values = sd.to_numpy(dtype=float)
    shown = np.ma.masked_invalid(values)
    cmap = blue_cmap().copy()
    cmap.set_bad("#f4f3f0")
    image = ax.imshow(shown, cmap=cmap, aspect="auto", vmin=0.0, vmax=np.nanpercentile(values, 95))
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            text = "–" if np.isnan(values[i, j]) else f"{values[i, j]:.2f}"
            dark = not np.isnan(values[i, j]) and values[i, j] > 0.6 * image.norm.vmax
            ax.text(j, i, text, ha="center", va="center", fontsize=7.5,
                    color="white" if dark else INK["primary"])  # fmt: skip

    def rank_text(name: str, kind: str) -> str:
        r = block_ranks[kind].reindex(sd.index)
        n_ranked = int(r.notna().any(axis=0).sum())
        covered = int(r.loc[name].notna().sum())
        mean = r.loc[name].mean()
        return "–" if covered == 0 else f"{mean:.1f}, {covered}/{n_ranked}"

    labels = [
        f"{name}  (" + " | ".join(rank_text(name, kind) for kind, _ in blocks) + ")"
        for name in sd.index
    ]
    ax.set_yticks(range(len(sd)), labels, fontsize=8.5)
    ax.set_xticks(range(len(materials)), [m.replace("_", " ") for m in materials], fontsize=8.5)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=30, ha="left")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    start = 0
    for k, (kind, cols) in enumerate(blocks):
        if k:
            ax.axvline(start - 0.5, color=SURFACE, lw=5)
        title = dict(GROUPS)[kind]
        ax.text(start + (len(cols) - 1) / 2, len(sd) - 0.2, title, ha="center", va="top",
                fontsize=10, fontweight="semibold", color=INK["secondary"], clip_on=False)  # fmt: skip
        start += len(cols)
    bar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.015)
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=0)
    bar.set_label(r"RMS of $\ln(1/s^*)$ on assessed tests (lower is better)", fontsize=10,
                  color=INK["secondary"])  # fmt: skip
    homepage_title(
        fig,
        "Which criterion misjudges the load least?",
        r"Load-equivalent error $\ln(1/s^*)$ of every method on every dataset ($s^*$: load factor "
        r"reproducing $N_{f,\exp}$ or reaching the limit). Row label: mean rank and materials "
        "ranked, for life | limit datasets",
        y=1 - 0.25 / height,
    )
    if hidden:
        fig.text(0.01, 0.005, "Not shown (fewer than %d methods scored): %s" % (MIN_RANKED, ", ".join(
            m.replace("_", " ") for m in hidden)), fontsize=9, color=INK["secondary"])  # fmt: skip
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(plot_load_error_ranking())
