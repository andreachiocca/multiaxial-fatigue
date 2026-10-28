"""Tveit, Ottosen, Erduran & Reyes (2025), Mech. Mater. 209:105452, doi:10.1016/j.mechmat.2025.105452.

Published model S-N curves of Figs. 5-6 with the Table 1 parameters (p.5-7), the torsional limit of
Eq. (24), infinite life inside the endurance surface and after shakedown under mean shear stress, and
the fatigue limits and numerical settings read from a constants file. Every case uses the Table 1
parameters.

The published curves are reproduced within 1 % in stress: at each digitized vertex (N, S) the published
life lies between the lives at 0.99 S and 1.01 S.
"""

import math

import numpy as np
import pytest

from _paper import load_table

from multifat.cli.pipeline import _constants_calibration_can_supply
from multifat.models.invariant_methods.tveit2025 import (
    MIN_BLOCKS,
    SETTINGS,
    Tveit2025,
    fatigue_index,
    tveit2025_life,
)

SQRT3 = math.sqrt(3.0)
TABLE1 = load_table("tveit2025_table1").set_index("material")
CURVES = load_table("tveit2025_fig5_fig6")


def _cycle(sigma_m=0.0, sigma_a=0.0, tau_m=0.0, tau_a=0.0, phase_shift=0.0):
    """S_mean, S_sin, S_cos of sigma_xx = sigma_m + sigma_a sin wt, tau_xy = tau_m + tau_a sin(wt + phase)."""
    S_mean, S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
    S_mean[0, 0] = sigma_m
    S_mean[0, 1] = S_mean[1, 0] = tau_m
    S_sin[0, 0] = sigma_a
    S_sin[0, 1] = S_sin[1, 0] = tau_a * math.cos(math.radians(phase_shift))
    S_cos[0, 1] = S_cos[1, 0] = tau_a * math.sin(math.radians(phase_shift))
    return S_mean, S_sin, S_cos


def _table1(material):
    row = TABLE1.loc[material]
    t_m1 = (row["lambda"] + row["f_m1"]) / SQRT3  # Eq. (24)
    return row["f_m1"], t_m1, row["A"], row["C"], row["K"], row["L"]


@pytest.mark.parametrize("row", list(CURVES.itertuples()), ids=lambda r: f"{r.panel}-{r.loading}")
def test_published_model_curves_within_one_percent_in_stress(row):
    constants = _table1(row.material)

    def life(S):
        cycle = _cycle(sigma_a=S) if row.loading == "bending" else _cycle(tau_a=S)
        return tveit2025_life(*cycle, *constants).Nf

    assert life(1.01 * row.S) < row.N < life(0.99 * row.S)


def test_eq24_torsional_limit_and_uniaxial_index():
    f_m1, t_m1, A = 199.0, 160.0, 0.2
    assert fatigue_index(*_cycle(tau_a=t_m1), f_m1, t_m1, A) == pytest.approx(1.0, rel=1e-12)
    index = fatigue_index(*_cycle(sigma_a=150.0), f_m1, t_m1, A)
    assert index == pytest.approx((1.0 + A) * 150.0 / f_m1, rel=1e-12)


def test_inside_endurance_surface_is_infinite():
    constants = _table1("Inoculated cast iron")
    cycle = _cycle(sigma_a=150.0)
    assert fatigue_index(*cycle, *constants[:3]) < 1.0
    life = tveit2025_life(*cycle, *constants)
    assert math.isinf(life.Nf)
    assert life.termination == "stabilized_extrapolation" and life.blocks == MIN_BLOCKS


def test_shakedown_under_mean_shear_is_infinite():
    constants = _table1("Inoculated cast iron")
    cycle = _cycle(tau_m=150.0, tau_a=50.0)
    assert fatigue_index(*cycle, *constants[:3]) > 1.0
    life = tveit2025_life(*cycle, *constants)
    assert math.isinf(life.Nf)
    assert life.termination == "stabilized_extrapolation" and life.blocks > MIN_BLOCKS


def test_case_model_returns_index_and_native_life():
    model = Tveit2025()
    constants = _table1("Inoculated cast iron")
    f_m1, t_m1, A, C, K, L = constants
    params = {
        "Sigm1": 2.0 * f_m1,  # workbook limits, replaced by the model's own
        "Taum1": 2.0 * t_m1,
        "Sigm1_TVEIT": f_m1,
        "Taum1_TVEIT": t_m1,
        "A_TVEIT": A,
        "C_TVEIT": C,
        "K_TVEIT": K,
        "L_TVEIT": L,
    }
    assert model.required_params() == {"A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT"}
    S_mean, S_sin, S_cos = _cycle(sigma_a=280.0, tau_a=170.0, phase_shift=90.0)
    harmonics = {"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos}
    result = model.evaluate_case(
        S0=S_mean + S_sin, E0=None, S1=S_mean - S_sin, E1=None, params=params, harmonics=harmonics
    )
    expected = tveit2025_life(S_mean, S_sin, S_cos, *constants)
    assert math.isfinite(expected.Nf)
    assert result.values[""] == fatigue_index(S_mean, S_sin, S_cos, f_m1, t_m1, A)
    assert result.aux["Nf"] == expected.Nf
    assert model.native_life(result.aux, params, model.life_fit([result.aux], [1.0], params)) == (
        expected.Nf
    )
    without_limits = {k: v for k, v in params.items() if k not in ("Sigm1", "Sigm1_TVEIT")}
    with pytest.raises(ValueError, match="Sigm1_TVEIT in the constants file or Sigm1"):
        model.evaluate_case(
            S0=S_mean + S_sin,
            E0=None,
            S1=S_mean - S_sin,
            E1=None,
            params=without_limits,
            harmonics=harmonics,
        )


def test_numerical_settings_from_the_constants_file():
    supplied = _constants_calibration_can_supply(Tveit2025())
    assert SETTINGS.keys() | {"Sigm1_TVEIT", "Taum1_TVEIT"} <= supplied
    constants = _table1("Inoculated cast iron")
    f_m1, t_m1, A, C, K, L = constants
    params = {
        "Sigm1": f_m1,
        "Taum1": t_m1,
        "A_TVEIT": A,
        "C_TVEIT": C,
        "K_TVEIT": K,
        "L_TVEIT": L,
        "max_blocks_TVEIT": 20.0,
    }
    S_mean, S_sin, S_cos = _cycle(sigma_a=280.0, tau_a=170.0, phase_shift=90.0)
    harmonics = {"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos}
    with pytest.warns(UserWarning, match="max_blocks_TVEIT"):
        result = Tveit2025().evaluate_case(
            S0=S_mean + S_sin,
            E0=None,
            S1=S_mean - S_sin,
            E1=None,
            params=params,
            harmonics=harmonics,
        )
    expected = tveit2025_life(S_mean, S_sin, S_cos, *constants, max_blocks=20)
    assert expected.termination == "max_block_extrapolation" and expected.blocks == 20
    assert result.aux["Nf"] == expected.Nf


def test_invalid_constants_name_the_symbol():
    with pytest.raises(ValueError, match="f_m1 > 0"):
        tveit2025_life(*_cycle(sigma_a=100.0), 0.0, 100.0, 0.2, 0.1, 1e-5, 2.0)
    with pytest.raises(ValueError, match="K > 0"):
        tveit2025_life(*_cycle(sigma_a=100.0), 150.0, 100.0, 0.2, 0.1, 0.0, 2.0)
    with pytest.raises(ValueError, match="n_points >= 3"):
        tveit2025_life(*_cycle(sigma_a=100.0), 150.0, 100.0, 0.2, 0.1, 1e-5, 2.0, n_points=2)
    with pytest.raises(ValueError, match="max_blocks >= 8"):
        tveit2025_life(*_cycle(sigma_a=100.0), 150.0, 100.0, 0.2, 0.1, 1e-5, 2.0, max_blocks=5)
