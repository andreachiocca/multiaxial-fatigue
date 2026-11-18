"""Shared academic plotting style and color palettes for multifat examples."""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_IMG_DIR = REPO_ROOT / "docs" / "images"
DATA_DIR = REPO_ROOT / "examples" / "data"

# Academic palette
PALETTE = {
    "navy": "#1B365D",
    "crimson": "#C0392B",
    "cerulean": "#2980B9",
    "emerald": "#27AE60",
    "purple": "#8E44AD",
    "amber": "#E67E22",
    "charcoal": "#2B2B2B",
    "grid": "#888888",
    "band_f2": "#2C3E50",
    "safe_fill": "#EBF5FB",
}


# Categorical palette of the batch figures (issue #71): blue, orange, aqua, violet, validated for
# colour-vision deficiencies on a light surface; the markers carry the identity without colour.
BATCH_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
BATCH_MARKERS = ["o", "s", "^", "D"]
# Diverging pair with a neutral midpoint: red = non-conservative, blue = conservative.
DIVERGING = ["#c8343b", "#f0efec", "#2a6fc4"]


def apply_academic_style() -> None:
    """Configure matplotlib with Times New Roman and Computer Modern math."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.edgecolor": PALETTE["charcoal"],
            "axes.linewidth": 1.0,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9.2,
            "legend.framealpha": 0.95,
            "figure.titlesize": 13.5,
            "figure.dpi": 300,
        }
    )


# ---------------------------------------------------------------------------------------------
# Homepage (README) visual system: one palette, one type family, one title layout for every
# figure. Categorical slots are the validated order of BATCH_COLOURS (adjacent CVD dE >= 9.2,
# normal-vision dE >= 27.6 on a light surface); aqua is below 3:1 contrast, so every series is
# also labelled directly. Magnitude uses one blue ramp, light (small) to dark (large).
# ---------------------------------------------------------------------------------------------
SERIES = BATCH_COLOURS  # blue, orange, aqua, violet
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#8a8984", "grid": "#e7e6e2"}
SURFACE = "#ffffff"
BLUE_RAMP = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
    "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]  # fmt: skip


def blue_cmap(start: str = "#f4f8fd"):
    """Sequential colormap on the blue ramp, from a near-surface tint to the darkest step."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("multifat_blue", [start, *BLUE_RAMP])


def apply_homepage_style() -> None:
    """Inter (DejaVu Sans fallback), recessive axes and grid, white surface."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Inter", "DejaVu Sans"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Inter",
            "mathtext.it": "Inter:italic",
            "mathtext.bf": "Inter:bold",
            "text.color": INK["primary"],
            "axes.labelcolor": INK["secondary"],
            "axes.edgecolor": INK["muted"],
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": INK["grid"],
            "grid.linewidth": 0.8,
            "axes.axisbelow": True,
            "axes.labelsize": 11,
            "xtick.color": INK["secondary"],
            "ytick.color": INK["secondary"],
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.5,
            "legend.frameon": False,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "figure.dpi": 220,
        }
    )


def homepage_title(fig, title: str, subtitle: str = "", y: float = 0.985) -> None:
    """Left-aligned bold title with a secondary subtitle, the same on every homepage figure."""
    fig.text(0.02, y, title, ha="left", va="top", fontsize=15, fontweight="bold",
             color=INK["primary"])  # fmt: skip
    if subtitle:  # a fixed 0.33 in below the title, whatever the figure height
        from matplotlib.transforms import ScaledTranslation

        below = fig.transFigure + ScaledTranslation(0.0, -0.33, fig.dpi_scale_trans)
        fig.text(0.02, y, subtitle, ha="left", va="top", fontsize=10.5,
                 color=INK["secondary"], transform=below)  # fmt: skip
