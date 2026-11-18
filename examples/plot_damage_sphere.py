"""Findley damage parameter over every material plane, drawn on the unit sphere of plane normals.

Each point of the sphere is the normal n of a candidate plane, coloured by the Findley damage
parameter DP(n) = tau_a + k_FI sigma_n,max on that plane; the critical plane (largest DP) is
marked. Same illustrative loading as the polar gallery figure: sigma_a = 240 MPa axial (x) with
tau_a = 150 MPa shear (xy), in phase and 90 deg out of phase, k_FI = 0.28.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from multifat.eval.evaluator import evaluate_planes, prepare_planes
from multifat.eval.plane_search import frames_from_normals
from multifat.models.registry import get_models
from examples.style import (
    DOCS_IMG_DIR,
    INK,
    SERIES,
    apply_homepage_style,
    blue_cmap,
    homepage_title,
)

SIGMA_A, TAU_A, K_FI = 240.0, 150.0, 0.28


def findley_on_sphere(S_sin, S_cos, n_polar: int = 121, n_azimuth: int = 241):
    """Unit normals on a (polar x azimuth) grid and the Findley DP of each plane."""
    polar = np.linspace(0.0, np.pi, n_polar)
    azimuth = np.linspace(0.0, 2.0 * np.pi, n_azimuth)
    P, A = np.meshgrid(polar, azimuth, indexing="ij")
    normals = np.stack([np.sin(P) * np.cos(A), np.sin(P) * np.sin(A), np.cos(P)], axis=-1)
    flat = normals.reshape(-1, 3)
    zero = np.zeros((3, 3))
    planes = prepare_planes(
        S0=S_sin, E0=zero, S1=-S_sin, E1=zero, R_list=frames_from_normals(flat),
        S_sin=S_sin, S_cos=S_cos,
    )  # fmt: skip
    damage, _ = evaluate_planes(get_models(["FIN"])[0]._m, planes, {"k_FI": K_FI})
    return normals, damage.reshape(P.shape)


def plot_damage_sphere(out_path: Path | None = None) -> Path:
    """Two spheres (in phase, 90 deg out of phase) coloured by the Findley damage parameter."""
    apply_homepage_style()
    out_path = out_path or DOCS_IMG_DIR / "critical_plane_sphere.png"
    shear = np.array([[0.0, TAU_A, 0.0], [TAU_A, 0.0, 0.0], [0.0, 0.0, 0.0]])
    axial = np.diag([SIGMA_A, 0.0, 0.0])
    cases = [
        (r"In phase, $\delta = 0^\circ$", axial + shear, np.zeros((3, 3))),
        (r"Out of phase, $\delta = 90^\circ$", axial, shear),
    ]
    results = [(title, *findley_on_sphere(s, c)) for title, s, c in cases]
    vmin = min(float(d.min()) for _, _, d in results)
    vmax = max(float(d.max()) for _, _, d in results)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = blue_cmap()

    fig = plt.figure(figsize=(11.0, 5.6), dpi=300)
    boxes = ([0.0, 0.13, 0.5, 0.62], [0.5, 0.13, 0.5, 0.62])  # 3D axes placed by hand
    for i, (title, normals, damage) in enumerate(results):
        ax = fig.add_axes(boxes[i], projection="3d")
        X, Y, Z = normals[..., 0], normals[..., 1], normals[..., 2]
        ax.plot_surface(
            X, Y, Z, facecolors=cmap(norm(damage)), rstride=1, cstride=1,
            linewidth=0, antialiased=False, shade=False,
        )  # fmt: skip
        k = np.unravel_index(int(np.argmax(damage)), damage.shape)
        n_cp = normals[k] if normals[k][0] >= 0 else -normals[k]  # n and -n: the same plane
        ax.computed_zorder = False  # draw in the order given: sphere, then axes and normal
        # Face the critical-plane normal towards the viewer
        azim, elev = -55.0, 22.0  # same viewpoint for both panels
        if n_cp @ np.array([np.cos(np.radians(azim)), np.sin(np.radians(azim)), 0.4]) < 0:
            n_cp = -n_cp  # n and -n are the same plane: show the one facing the viewer
        camera = np.array(
            [
                np.cos(np.radians(elev)) * np.cos(np.radians(azim)),
                np.cos(np.radians(elev)) * np.sin(np.radians(azim)),
                np.sin(np.radians(elev)),
            ]
        )
        for vec, label in (((1, 0, 0), r"$x$ (axial)"), ((0, 1, 0), r"$y$"), ((0, 0, 1), r"$z$")):
            vec = np.asarray(vec, dtype=float)
            vec = vec if vec @ camera >= 0 else -vec  # the half-axis facing the viewer
            sign = "" if vec.sum() > 0 else r"$-$"
            ax.plot(*zip(1.0 * vec, 1.35 * vec), color=INK["muted"], linewidth=1.0, zorder=5)
            ax.text(*(1.45 * vec), sign + label, fontsize=9, color=INK["secondary"], zorder=6)
        tip = 1.5 * n_cp
        ax.plot(*zip(n_cp, tip), color=SERIES[1], linewidth=2.6, zorder=10)
        ax.scatter(
            *tip,
            color=SERIES[1],
            s=34,
            depthshade=False,
            zorder=11,
            edgecolors="white",
            linewidths=1.2,
        )
        ax.text(
            *(1.68 * n_cp), r"$\mathbf{n}^*$", color=INK["primary"],
            fontsize=13, ha="center", va="center", zorder=12,
        )  # fmt: skip
        x0 = boxes[i][0] + boxes[i][2] / 2
        fig.text(x0, 0.825, title, ha="center", va="bottom", fontsize=11.5, color=INK["primary"])
        fig.text(x0, 0.82, f"max DP = {damage.max():.0f} MPa", ha="center", va="top",
                 fontsize=10, color=INK["secondary"])  # fmt: skip
        ax.view_init(elev=elev, azim=azim)
        ax.set_box_aspect((1, 1, 1), zoom=1.95)
        ax.set_axis_off()
        lim = 1.7
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
    bar = fig.colorbar(
        cm.ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([0.3, 0.1, 0.4, 0.025]),
        orientation="horizontal",
    )  # fmt: skip
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=0)
    bar.set_label(
        r"Findley DP $= \tau_a + k\,\sigma_{n,\max}$ on the plane of normal $\mathbf{n}$ [MPa]",
        color=INK["secondary"],
    )
    homepage_title(
        fig,
        "Every material plane, one damage value",
        rf"Findley criterion over all plane normals, $\sigma_a = {SIGMA_A:.0f}$ MPa, "
        rf"$\tau_a = {TAU_A:.0f}$ MPa, $k_{{FI}} = {K_FI}$ (illustrative loading); "
        r"$\mathbf{n}^*$ is the critical plane",
    )
    fig.savefig(out_path, facecolor="#FFFFFF")  # no tight box: it clips the 3D labels
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    print(plot_damage_sphere())
