"""34CrNiMo6_Pallares: fatigue-limit dataset from Pallares-Santasmartas et al. 2018 (decision d120).

Sources: [A] Metals 8 (2018) 213, doi:10.3390/met8040213, Tables 2, 3 and 5; [T] Int J Fatigue
113 (2018) 54-68, doi:10.1016/j.ijfatigue.2018.04.008, Tables 3, 4 and 5. Both transcribed in
Experimental_Data/sources/Pallares2018/34crnimo6_limits.csv; tools/build_pallares_34crnimo6.py
documents the scope and builds the workbook.

The transcription is checked against numbers the papers derive from it: the Crossland column of
[A] Table 5, the "no effect hypothesis" column of [T] Table 4, and the constant c of [T] Table 5.
"""

import math
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from build_pallares_34crnimo6 import SOURCE, build, read_transcription, stress_ratio  # noqa: E402

from multifat.imp.classification import classify_dataset  # noqa: E402
from multifat.imp.experimental_data import load_material_workbook  # noqa: E402
from multifat.models.registry import get_models  # noqa: E402

WORKBOOK = ROOT / "Experimental_Data" / "20_34CrNiMo6_Pallares.xlsx"
PROPS, LIMITS = read_transcription()
SIGM1, TAUM1 = 615.0, 432.5


def _rows(path, sheet):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return [tuple(r) for r in wb[sheet].iter_rows(values_only=True)]


def test_committed_workbook_is_the_build(tmp_path):
    info = build(SOURCE, tmp_path / "pallares.xlsx")
    for sheet in ("Summary", "Dataset"):
        assert _rows(tmp_path / "pallares.xlsx", sheet) == _rows(WORKBOOK, sheet)
    assert (info["Sigm1"], info["Taum1"], info["rows"]) == (SIGM1, TAUM1, 9)


def test_printed_stress_ratios_follow_from_the_stresses():
    """[A] Table 3 prints R next to sigma_m and sigma_a."""
    for lim in LIMITS:
        if lim["R"] is not None:
            assert stress_ratio(lim) == pytest.approx(lim["R"], abs=0.005)


def test_constant_c_of_the_proposed_method():
    """[T] Eq. (18c) and Table 5: c = sigma_-1 / (3 - sqrt(3) kappa) = 1145.1 MPa."""
    kappa = SIGM1 / TAUM1
    assert SIGM1 / (3.0 - math.sqrt(3.0) * kappa) == pytest.approx(1145.1, abs=0.05)


def _crossland_amplitude(sigma_m):
    """Axial amplitude at the Crossland limit for a fixed mean, with the registered CROSSLAND."""
    (model,) = get_models(["CROSSLAND"])
    params = {"Sigm1": SIGM1, "Taum1": TAUM1}

    def dp(sa):
        S_mean, S_sin = np.zeros((3, 3)), np.zeros((3, 3))
        S_mean[1, 1], S_sin[1, 1] = sigma_m, sa
        h = {"S_mean": S_mean, "S_sin": S_sin, "S_cos": np.zeros((3, 3))}
        return model.evaluate_case(S0=None, E0=None, S1=None, E1=None, params=params,
                                   harmonics=h).values[""]  # fmt: skip

    lo, hi = 1.0, 2000.0
    for _ in range(80):  # the index is monotonic in the amplitude
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dp(mid) < 1.0 else (lo, mid)
    return 0.5 * (lo + hi)


@pytest.mark.parametrize("sigma_m, printed", [(-216.0, -1.0), (181.0, -7.5), (522.0, -10.5)])
def test_crossland_column_of_metals_table5(sigma_m, printed):
    """[A] Table 5: relative error (sigma_a,exp - sigma_a,Crossland) / sigma_a,exp in %."""
    sa = next(lim["sa"] for lim in LIMITS if lim["load"] == "axial" and lim["sm"] == sigma_m)
    error = 100.0 * (sa - _crossland_amplitude(sigma_m)) / sa
    assert error == pytest.approx(printed, abs=0.05)


@pytest.mark.parametrize(
    "tau_m, printed",
    [
        (150.0, -0.9),
        (250.0, -5.7),
        (350.0, -9.7),
        pytest.param(
            500.0, -28.2,
            marks=pytest.mark.xfail(strict=True, reason="the printed tau_a = 337.1 and "
                                    "tau_-1 = 432.5 MPa give -28.30 %, printed -28.2"),
        ),  # fmt: skip
    ],
)
def test_no_effect_column_of_ijf_table4(tau_m, printed):
    """[T] Table 4, "no effect hypothesis": (tau_a,exp - tau_-1) / tau_a,exp in %."""
    ta = next(lim["ta"] for lim in LIMITS if lim["load"] == "torsion" and lim["tm"] == tau_m)
    assert 100.0 * (ta - TAUM1) / ta == pytest.approx(printed, abs=0.05)


def test_imported_as_a_fatigue_limit_dataset():
    block = load_material_workbook(WORKBOOK)["34CrNiMo6_Pallares"]
    m = block["material"]
    assert (m["E_MPa"], m["nu"]) == (None, None)  # neither paper gives them
    assert (m["Tensile_strength_MPa"], m["Yield_strength_MPa"]) == (1210.0, 1084.0)
    assert (m["Tensile_fatigue_limit_MPa"], m["Torsional_fatigue_limit_MPa"]) == (SIGM1, TAUM1)
    (sigm1,) = block["fatigue_limit_sources"]["Sigm1"]
    assert sigm1["loading"] == "push-pull" and sigm1["reference_cycles"] == 2e6
    assert not block["identity_issues"]
    assert (len(block["cases"]), len(block["runouts"])) == (0, 9)
    assert classify_dataset(block, "34CrNiMo6_Pallares")["workflow"] == "fatigue_limit"
