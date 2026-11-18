"""One selector of the fully reversed uniaxial design tests (#77).

The master curve (``extract_uniaxial_tension``), the Basquin constants
(``reversed_uniaxial_tests``) and the plane-constant calibration (``_is_axial_tension_case``)
must pick the same tests; before #77 they disagreed on these perturbed copies of a real case.
"""

import copy
from pathlib import Path

import numpy as np
import pytest

from multifat.cli.calibrate import INITIAL_RTOL_R, _is_axial_tension_case
from multifat.eval.material_fits import reversed_uniaxial_tests
from multifat.imp.experimental_data import (
    DESIGN_RTOL_R,
    extract_uniaxial_tension,
    load_material_workbook,
)

DATA = Path(__file__).resolve().parents[1] / "Experimental_Data"


@pytest.fixture(scope="module")
def tension_case():
    data = load_material_workbook(DATA / "6_42CrMo4_Santus.xlsx")
    block = data["42CrMo4_Santus"]
    case_id = extract_uniaxial_tension(data)["42CrMo4_Santus"]["points"][0]["case_id"]
    return block, case_id


def _mean_stress(case, fraction):
    a = case[1][1]["stress"][1, 1]
    for step in (1, 2):
        case[1][step]["stress"] = case[1][step]["stress"].copy()
        case[1][step]["stress"][1, 1] += fraction * a
    case["harmonics"]["S_mean"] = case["harmonics"]["S_mean"].copy()
    case["harmonics"]["S_mean"][1, 1] += fraction * a


def _swap_steps(case, _):
    case[1][1], case[1][2] = case[1][2], case[1][1]


def _parasitic_shear(case, fraction):
    a = case[1][1]["stress"][1, 1]
    for step, sign in ((1, 1.0), (2, -1.0)):
        S = case[1][step]["stress"].copy()
        S[0, 1] = S[1, 0] = sign * fraction * a
        case[1][step]["stress"] = S
    for key in ("S_sin",):
        H = case["harmonics"][key].copy()
        H[0, 1] = H[1, 0] = fraction * a
        case["harmonics"][key] = H


@pytest.mark.parametrize(
    "perturb, amount, expected",
    [
        (_mean_stress, 0.005, False),  # R = -0.990: outside 1e-3 everywhere
        (_swap_steps, None, True),  # load-step order does not matter
        (_parasitic_shear, 0.0005, True),  # 0.05 % shear: within the 1e-3 uniaxial tolerance
        (_parasitic_shear, 0.005, False),
    ],
)
def test_the_three_selectors_agree(tension_case, perturb, amount, expected):
    block, case_id = tension_case
    case = copy.deepcopy(block["cases"][case_id])
    perturb(case, amount)
    data = {"M": {**block, "cases": {case_id: case}}}
    master = bool(extract_uniaxial_tension(data).get("M", {}).get("points"))
    basquin = case_id in reversed_uniaxial_tests(data, "M")["tension"][0]
    calibration = _is_axial_tension_case(case, rtol_R=INITIAL_RTOL_R)[0]
    assert master == basquin == calibration == expected


def test_calibration_starts_from_the_design_set_tolerance():
    assert INITIAL_RTOL_R == DESIGN_RTOL_R
    assert np.isfinite(DESIGN_RTOL_R) and DESIGN_RTOL_R <= 1e-3
