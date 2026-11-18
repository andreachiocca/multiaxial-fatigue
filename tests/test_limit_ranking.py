"""Fatigue-limit datasets in the load-error ranking (decision d110).

A limit state enters with Error_load = ln(1 + FI_load); the pure R = -1 axial limit states are
the reference and are left out, as the design failures of a life dataset are.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from examples import build_load_error_ranking as br  # noqa: E402
from multifat.cli import run as runmod  # noqa: E402
from multifat.report.subsets import dataset_fingerprint  # noqa: E402


def _case(mean=0.0, axial=0.0, shear=0.0, cos_shear=0.0):
    z = np.zeros((3, 3))
    s_mean, s_sin, s_cos = z.copy(), z.copy(), z.copy()
    s_mean[1, 1] = mean
    s_sin[1, 1] = axial
    s_sin[0, 1] = s_sin[1, 0] = shear
    s_cos[0, 1] = s_cos[1, 0] = cos_shear
    return {"harmonics": {"S_mean": s_mean, "S_sin": s_sin, "S_cos": s_cos}}


@pytest.mark.parametrize(
    "case, expected",
    [
        (_case(axial=652.0), True),
        (_case(mean=300.0, axial=630.0), False),  # axial limit with a mean stress: assessed
        (_case(shear=411.0), False),  # torsion: assessed, as torsion failures are
        (_case(axial=442.0, shear=255.0), False),
        (_case(axial=567.0, cos_shear=328.0), False),  # 90 deg
    ],
)
def test_only_pure_reversed_axial_limits_are_the_reference(case, expected):
    assert br.is_axial_reference(case) is expected


def test_limit_error_is_the_log_of_the_load_factor(tmp_path):
    fi = [0.0, -0.2, 0.25, 1.0]  # FI_load = 1/s - 1
    path = tmp_path / "M_X_FI.csv"
    pd.DataFrame(dict(Material="X", Method="M", Specimen_id=list("abcd"), FI_load=fi)).to_csv(
        path, index=False
    )
    errors = br.limit_errors(path, [True, False, False, False], list("abcd"))
    s = 1.0 / (1.0 + np.array(fi[1:]))
    assert errors == pytest.approx(np.log(1.0 / s))
    with pytest.raises(ValueError, match="do not match"):
        br.limit_errors(path, [True, False, False], list("abc"))


def test_delahay_reference_is_the_plane_bending_limit():
    data, mat = runmod.load_material("Ti6Al4V_Delahay")
    flags = [br.is_axial_reference(c) for c in data[mat]["runouts"].values()]
    assert sum(flags) == 1  # plane bending 652 MPa; torsion and the two combined limits assessed


def test_runout_fingerprint_is_opt_in():
    data, mat = runmod.load_material("Ti6Al4V_Delahay")
    block = data[mat]
    assert dataset_fingerprint(block) == dataset_fingerprint({**block, "runouts": {}})
    assert dataset_fingerprint(block, runouts=True) != dataset_fingerprint(
        {**block, "runouts": {}}, runouts=True
    )


def test_committed_ranking_keeps_the_life_rows_and_leaves_out_data_under_review():
    table = pd.read_csv(br.OUT)
    assert set(table["Assessment"]) == {"fatigue_life", "fatigue_limit"}
    assert not set(br.UNDER_REVIEW) & set(table["Material"])
    assert (table["n"] >= 3).all()


def test_runout_fingerprint_ignores_last_bit_differences():
    """Windows and Linux libm differ in the last bits of cos/sin (cos 90 deg ~ 6e-17)."""
    import copy

    data, mat = runmod.load_material("30NCD16")
    block = data[mat]
    other = copy.deepcopy(block)
    for case in other["runouts"].values():
        for key in ("S_sin", "S_cos"):
            case["harmonics"][key] = np.asarray(case["harmonics"][key]) * (1 + 1e-12) + 1e-15
    assert dataset_fingerprint(other, runouts=True) == dataset_fingerprint(block, runouts=True)
    case = next(iter(other["runouts"].values()))
    case["harmonics"]["S_sin"] = np.asarray(case["harmonics"]["S_sin"]) + 0.01  # a real change
    assert dataset_fingerprint(other, runouts=True) != dataset_fingerprint(block, runouts=True)
