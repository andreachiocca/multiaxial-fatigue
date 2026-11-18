"""Calibration on tests outside R = -1 is accepted but flagged (issue #77)."""

import numpy as np
import pytest

from multifat.cli import calibrate
from multifat.cli.pipeline import Outcome, _annotate_calibration, _calibration_state
from multifat.models.registry import get_models

R_MINUS_1 = dict(fallback=False, tension_rtol_R=calibrate.INITIAL_RTOL_R,
                 torsion_rtol_R=calibrate.INITIAL_RTOL_R,
                 tension_any_R=False, torsion_any_R=False)  # fmt: skip


def test_no_flag_on_r_minus_one_tests():
    assert calibrate.r_fallback_text(R_MINUS_1) == ""
    assert "r_fallback" not in _calibration_state({"selection": R_MINUS_1}, set())


@pytest.mark.parametrize(
    "selection, expected",
    [
        (dict(R_MINUS_1, fallback=True, torsion_rtol_R=0.1),
         "calibrated on torsion tests with |R + 1| up to 0.1"),
        (dict(R_MINUS_1, fallback=True, tension_any_R=True, tension_rtol_R=None),
         "calibrated on tension tests of any stress ratio"),
    ],
)  # fmt: skip
def test_fallback_is_described(selection, expected):
    assert calibrate.r_fallback_text(selection) == expected


def test_model_using_a_fallback_constant_is_flagged():
    state = _calibration_state(
        {"selection": {"r_fallback": "calibrated on tension tests of any stress ratio"}}, set()
    )
    model = get_models(["FS"])[0]
    flagged = _annotate_calibration(Outcome("FS", "ok", calibrated="k_FS"), model, state)
    assert flagged.status == "ok (R fallback)"
    assert "k_FS calibrated on tension tests of any stress ratio" in flagged.notes
    # a model that did not adopt a plane constant keeps its status
    other = _annotate_calibration(Outcome("SWT", "ok", calibrated="sigma_f,b_fat"), model, state)
    assert other.status == "ok"


def test_selection_relaxes_and_records_the_tolerance():
    """No test within the initial tolerance: the selection relaxes it and says so."""

    def selector(case, *, r_target=-1.0, rtol_R=0.02, **_):
        R = case["R"]
        if r_target is None:
            return True, R
        return bool(np.isclose(R, r_target, rtol=rtol_R)), R

    cases = {0: {"R": -0.92}, 1: {"R": -0.5}}
    with pytest.warns(RuntimeWarning, match="relaxed rtol_R=0.1"):
        ids, info = calibrate._pick_cases_with_relaxation(
            cases, selector, selector_kwargs=dict(r_target=-1.0, rtol_R=0.02)
        )
    assert ids == [0] and info["fallback_used"] and info["rtol_R_used"] == 0.1
    with pytest.warns(RuntimeWarning, match="ANY R"):
        ids, info = calibrate._pick_cases_with_relaxation(
            {1: {"R": 0.1}}, selector, selector_kwargs=dict(r_target=-1.0, rtol_R=0.02)
        )
    assert ids == [1] and info["any_R_used"]
