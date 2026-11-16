"""Homepage figures of the README, drawn in the one homepage style of ``examples/style.py``.

1. Critical-plane damage on the sphere of plane normals (header)
2. Biaxial fatigue limit envelopes in normalized stress space
3. Method ranking by the load-equivalent error across materials
4. Cross-batch transfer of the constants (42CrMo4)
"""

import sys
from pathlib import Path
from typing import List

# Ensure repository root is on sys.path for direct script execution
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.plot_batch_transfer import plot_batch_transfer  # noqa: E402
from examples.plot_damage_sphere import plot_damage_sphere  # noqa: E402
from examples.plot_fatigue_limit import plot_fatigue_limit_envelopes  # noqa: E402
from examples.plot_load_error_ranking import plot_load_error_ranking  # noqa: E402

DOCS_IMG_DIR = REPO_ROOT / "docs" / "images"
FIGURES = (
    ("critical_plane_sphere.png", plot_damage_sphere, "Critical-plane damage sphere"),
    ("fatigue_limit_envelopes.png", plot_fatigue_limit_envelopes, "Fatigue limit envelopes"),
    ("load_error_ranking.png", plot_load_error_ranking, "Load-equivalent error ranking"),
    ("batch_transfer.png", plot_batch_transfer, "Cross-batch transfer"),
)


def generate_all_gallery_plots(target_dir: Path = DOCS_IMG_DIR) -> List[Path]:
    """Draw every homepage figure into target_dir."""
    target_dir.mkdir(parents=True, exist_ok=True)
    generated: List[Path] = []
    for k, (name, plot, label) in enumerate(FIGURES, start=1):
        path = plot(out_path=target_dir / name)
        print(f"  [{k}/{len(FIGURES)}] {label}: {path}")
        generated.append(path)
    return generated


if __name__ == "__main__":
    generate_all_gallery_plots()
