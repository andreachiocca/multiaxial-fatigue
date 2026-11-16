"""README figures: the damage sphere agrees with the production critical-plane search."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from examples import homepage_data  # noqa: E402
from examples.plot_damage_sphere import K_FI, SIGMA_A, TAU_A, findley_on_sphere  # noqa: E402
from multifat.eval.evaluator import best_cp_for_case  # noqa: E402
from multifat.models.cp_methods.fin import Findley  # noqa: E402
from multifat.planes.sampling import build_R_list  # noqa: E402


@pytest.mark.parametrize("phase", [0, 90])
def test_sphere_maximum_matches_the_refined_search(phase):
    shear = np.array([[0.0, TAU_A, 0.0], [TAU_A, 0.0, 0.0], [0.0, 0.0, 0.0]])
    axial = np.diag([SIGMA_A, 0.0, 0.0])
    S_sin, S_cos = (axial + shear, np.zeros((3, 3))) if phase == 0 else (axial, shear)
    _, damage = findley_on_sphere(S_sin, S_cos)
    zero = np.zeros((3, 3))
    best = best_cp_for_case(
        S0=S_sin, E0=zero, S1=-S_sin, E1=zero, R_list=build_R_list(), method=Findley(),
        params={"k_FI": K_FI}, S_sin=S_sin, S_cos=S_cos,
    )  # fmt: skip
    # the 1.5-degree sphere grid reaches the refined maximum to a few parts in ten thousand
    assert damage.max() <= best.ext * (1 + 1e-9)
    assert damage.max() == pytest.approx(best.ext, rel=2e-3)


def test_ranking_data_is_consistent():
    table = pd.read_csv(
        Path(__file__).resolve().parents[1] / "examples" / "data" / "load_error_ranking.csv"
    )
    assert {"Material", "Method", "n", "Error_load_mean", "Error_load_rms"} <= set(table.columns)
    rms2 = table.Error_load_rms**2
    var = table.Error_load_sd**2 * (table.n - 1) / table.n + table.Error_load_mean**2
    assert np.allclose(rms2, var, rtol=1e-4)  # RMS^2 = mean^2 + population variance
    assert not table.Method.str.endswith("_ext").any()


def test_committed_homepage_data_is_current():
    """A new method, changed revision or changed workbook needs homepage CSVs regenerated."""
    problems = homepage_data.check()
    assert not problems, (
        "; ".join(problems)
        + " -> regenerate examples/data/load_error_ranking.csv with build_load_error_ranking.py"
    )


def test_homepage_manifest_detects_staleness(tmp_path, monkeypatch):
    csv_path = tmp_path / "sample.csv"
    csv_path.write_text("Material,Method\nM,FIN\n", encoding="utf-8")
    monkeypatch.setattr(homepage_data, "registered_models", lambda: ["FIN", "SWT"])
    monkeypatch.setattr(
        homepage_data, "dataset_fingerprints", lambda mats: {m: f"fp-{m}" for m in mats}
    )
    monkeypatch.setattr(homepage_data, "_implementation_revision", lambda: "2026-11-11")
    mpath = homepage_data.write_manifest(csv_path, ["M"], ["FIN"])

    assert homepage_data.check_manifest(csv_path, ["M"]) == []
    assert any(
        "models added ['CROSSLAND']" in p
        for p in homepage_data.check_manifest(
            csv_path, ["M"], current_models=["CROSSLAND", "FIN", "SWT"]
        )
    )
    assert any(
        "implementation revision" in p
        for p in homepage_data.check_manifest(csv_path, ["M"], current_revision="2026-12-01")
    )
    assert any(
        "datasets changed ['M']" in p
        for p in homepage_data.check_manifest(
            csv_path, ["M"], current_fingerprints={"M": "fp-changed"}
        )
    )

    csv_path.write_text("Material,Method\nM,SWT\n", encoding="utf-8")
    assert any("sha256" in p for p in homepage_data.check_manifest(csv_path, ["M"]))
    mpath.unlink()
    assert any("missing manifest" in p for p in homepage_data.check_manifest(csv_path, ["M"]))


def test_readme_batch_transfer_numbers_follow_the_data():
    """The percentages quoted in the README section are those of the committed batch data."""
    from examples.plot_batch_transfer import transfer_summary

    readme = (Path(__file__).resolve().parents[1] / "README.md").read_text(encoding="utf-8")
    s = transfer_summary().set_index(["design", "target"])["median"]

    def pct(design, target):
        return round(100 * (np.exp(s[(f"42CrMo4_{design}", f"42CrMo4_{target}")]) - 1))

    lo, hi = sorted((pct("Santus", "Reis"), pct("Santus", "FABEST")))
    assert (f"+{lo} %" if lo == hi else f"+{lo}–{hi} %") in readme
    assert f"−{-pct('Reis', 'Santus')} % from Reis" in readme
    assert f"−{-pct('FABEST', 'Santus')} % from FABEST" in readme
    assert f"(−{-pct('FABEST', 'Reis')} % and +{pct('Reis', 'FABEST')} %)" in readme
