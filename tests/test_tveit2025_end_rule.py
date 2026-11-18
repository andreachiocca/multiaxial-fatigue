"""End of the TVEIT2025 block integration: shakedown, cyclic steady state and the failure fraction.

The integration is checked against 50,000-block histories run without any stopping rule
(tests/data/tveit2025_long_histories.csv, stage-1 constants). Every load that shakes down must be
found infinite at the block where the damage increment first vanishes, and every load in a steady
cycle must give the life of the long history. Comparing only two consecutive blocks, the earlier
rule stopped during these shakedowns and returned 9e7 to 2.4e9 cycles.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multifat.models.invariant_methods.tveit2025 import N_POINTS, Tveit2025, tveit2025_life

ROOT = Path(__file__).resolve().parents[1]
HISTORIES = pd.read_csv(ROOT / "tests" / "data" / "tveit2025_long_histories.csv", comment="#")


def _constants(material):
    params = json.loads(
        (ROOT / "Calibrations" / material / "TVEIT2025_stage1.json").read_text(encoding="utf-8")
    )["hand_entered_constants"]
    keys = ("Sigm1_TVEIT", "Taum1_TVEIT", "A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT")
    return [params[k] for k in keys]


def _cycle(loading, R, amplitude):
    S_mean, S_sin = np.zeros((3, 3)), np.zeros((3, 3))
    i, j = (1, 1) if loading == "axial" else (0, 1)
    S_mean[i, j] = S_mean[j, i] = amplitude * (1 + R) / (1 - R)
    S_sin[i, j] = S_sin[j, i] = amplitude
    return S_mean, S_sin, np.zeros((3, 3))


@pytest.mark.parametrize(
    "row",
    list(HISTORIES.itertuples()),
    ids=lambda r: f"{r.material}-{r.loading}-R{r.R:g}-{r.amplitude:g}",
)
def test_long_history_is_reproduced(row):
    life = tveit2025_life(*_cycle(row.loading, row.R, row.amplitude), *_constants(row.material))
    if not math.isnan(row.shakedown_block):
        assert life.termination == "shakedown" and math.isinf(life.Nf)
        # The test needs two blocks; a few blocks of tolerance cover rounding across platforms
        assert abs(life.blocks - max(int(row.shakedown_block), 2)) <= 2
    else:
        expected = 50000 + (1 - row.D_50000) / row.dD_50000
        assert life.termination == "cyclic_steady_state"
        # 2 %: at 50,000 blocks the damage per block of Ti6Al4V at 510 MPa still drifts by 1.4 %
        assert life.Nf == pytest.approx(expected, rel=0.02)


def test_creeping_back_stress_is_not_taken_for_a_steady_state():
    """Al7075_T6 stage-1 constants, torsion at its limit with a small constant out-of-plane shear.

    The back-stress swings in every block (its net drift over 200 blocks is 0.9 % of its path
    length) but creeps one way until the damage stops at block 2446. The steady-state test of #67
    took block 400 for a steady state and extrapolated 4.6e9 cycles (found by an automated review).
    """
    S_mean, S_sin, S_cos = _cycle("torsion", -1.0, 99.996)
    S_mean[1, 2] = S_mean[2, 1] = 0.895
    life = tveit2025_life(S_mean, S_sin, S_cos, *_constants("Al7075_T6"))
    assert life.termination == "shakedown" and math.isinf(life.Nf)
    assert abs(life.blocks - 2446) <= 2


def test_unresolved_decay_ends_at_the_block_limit_with_a_finite_life():
    """A damage decay not resolved within max_blocks is extrapolated, conservatively, with a warning.

    Ti6Al4V at 495 MPa shakes down at block 7219 (long history); stopped at 3000 blocks, it is still
    decaying, and the life is extrapolated from the mean damage of the last 200 blocks. The model
    (the case-level evaluation the pipeline calls) warns about it.
    """
    S_mean, S_sin, S_cos = _cycle("axial", -1.0, 495.0)
    life = tveit2025_life(S_mean, S_sin, S_cos, *_constants("Ti6Al4V_Meneghetti"), max_blocks=3000)
    assert life.termination == "max_block_extrapolation" and life.blocks == 3000
    assert 3000 < life.Nf < math.inf
    keys = ("Sigm1_TVEIT", "Taum1_TVEIT", "A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT")
    params = dict(zip(keys, _constants("Ti6Al4V_Meneghetti")), max_blocks_TVEIT=3000)
    with pytest.warns(UserWarning, match="extrapolated after max_blocks_TVEIT blocks"):
        result = Tveit2025().evaluate_case(
            S0=S_mean + S_sin,
            E0=None,
            S1=S_mean - S_sin,
            E1=None,
            params=params,
            harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos},
        )
    assert result.aux["Nf"] == life.Nf


def test_failure_inside_the_first_increment_counts_less_than_one_increment():
    """The life ends inside the increment that crosses D = 1, not one increment later."""
    S_mean, S_sin, S_cos = _cycle("axial", -1.0, 1e5)  # damaging from the first increment on
    life = tveit2025_life(S_mean, S_sin, S_cos, 100.0, 80.0, 0.2, 0.1, 0.1, 50.0)
    assert life.termination == "explicit_failure" and life.blocks == 1
    assert 0.0 < life.Nf < 1.0 / (N_POINTS - 1)
