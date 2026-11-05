"""Master-curve (S-N) models of issue #38: Kohout & Vechet (2001) form (H2) and Stromeyer."""

from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd
import pytest

from multifat.eval.fitting import fit_power_law_with_survival_std
from multifat.eval.sn_curves import curve_dp, curve_life, fit_sn_curve

N = np.logspace(4.0, 7.5, 25)


def _kv(Nf, sigma_inf=300.0, C=2.0e6, b=-0.12):
    """Kohout & Vechet (2001) Eq. (H2): sigma = sigma_inf (N / (N + C))^b."""
    return sigma_inf * (Nf / (Nf + C)) ** b


def test_basquin_is_the_established_fit_bit_for_bit():
    DP = 500.0 * N**-0.1 * np.exp(np.random.default_rng(0).normal(0, 0.05, N.size))
    fit = fit_sn_curve(N, DP, stdnum=1.3)
    reference = fit_power_law_with_survival_std(N, DP, stdnum=1.3)
    assert {k: fit[k] for k in reference} == reference
    assert fit["sn_model"] == "basquin" and fit["valid"]


def test_kohout_vechet_recovers_its_parameters_and_geometry():
    fit = fit_sn_curve(N, _kv(N), model="kohout_vechet")
    assert fit["valid"], fit["reason"]
    p = fit["params"]
    assert p["sigma_inf"] == pytest.approx(300.0, rel=1e-6)
    assert p["C"] == pytest.approx(2.0e6, rel=1e-5)
    assert p["b"] == pytest.approx(-0.12, rel=1e-6)
    # Eq. (10): sigma_inf = a C^b, so the Basquin asymptote is a N^b for N << C
    assert fit["A"] == pytest.approx(300.0 * 2.0e6**0.12, rel=1e-5)
    small = np.array([10.0, 100.0])
    assert curve_dp(fit, small) == pytest.approx(fit["A"] * small**-0.12, rel=1e-4)
    # Log-log slope b/2 at N = C (paper, Section 5)
    h = 1e-4
    slope = (math.log(_kv(2.0e6 * math.exp(h))) - math.log(_kv(2.0e6 * math.exp(-h)))) / (2 * h)
    assert slope == pytest.approx(-0.06, rel=1e-6)


@pytest.mark.parametrize("model", ["kohout_vechet", "stromeyer"])
def test_inverse_round_trip_and_infinite_life_below_the_limit(model):
    DP = _kv(N) * np.exp(np.random.default_rng(1).normal(0.0, 0.02, N.size))
    fit = fit_sn_curve(N, DP, model=model)
    assert fit["valid"], fit["reason"]
    lives = np.array([3e4, 2e5, 5e6])
    assert curve_life(fit, curve_dp(fit, lives)) == pytest.approx(lives, rel=1e-9)
    limit = fit["params"]["sigma_inf"] * fit["shift"]
    assert np.all(np.isinf(curve_life(fit, np.array([0.5 * limit, limit]))))


def test_survival_shift_is_a_vertical_shift_in_ln_dp():
    DP = _kv(N) * np.exp(np.random.default_rng(2).normal(0.0, 0.03, N.size))
    median = fit_sn_curve(N, DP, model="kohout_vechet")
    shifted = fit_sn_curve(N, DP, stdnum=1.2815516, model="kohout_vechet")
    ratio = curve_dp(shifted, N) / curve_dp(median, N)
    assert ratio == pytest.approx(np.exp(-1.2815516 * median["STD_ln_DP"]), rel=1e-9)


def test_unidentifiable_fits_are_rejected_with_a_reason():
    straight = 500.0 * N**-0.1  # no knee: C runs to its bound
    fit = fit_sn_curve(N, straight, model="kohout_vechet")
    assert not fit["valid"] and "not identifiable" in fit["reason"]
    few = fit_sn_curve(N[:4], _kv(N[:4]), model="kohout_vechet")
    assert not few["valid"] and "at least 6" in few["reason"]
    with pytest.raises(ValueError, match="unknown S-N model"):
        fit_sn_curve(N, _kv(N), model="palmgren")


def test_pipeline_records_the_model_and_compare_refuses_mixtures(tmp_path):
    from multifat.cli.compare import collect_all_results
    from multifat.cli.pipeline import run_workbook
    from workbook_factory import write_workbook

    wb = write_workbook(tmp_path / "Sn.xlsx")
    run_workbook(wb, models=["CROSSLAND"], results_dir=tmp_path / "r" / "a")
    outcomes = run_workbook(
        wb, models=["CROSSLAND"], results_dir=tmp_path / "r" / "b", sn_model="stromeyer"
    )
    side = json.loads((tmp_path / "r" / "b" / "Sn_calibration.json").read_text(encoding="utf-8"))
    assert side["provenance"]["settings"]["master_curve_model"] == "stromeyer"
    if outcomes[0].status == "ok":
        df = pd.read_csv(tmp_path / "r" / "b" / "CROSSLAND_Sn.csv")
        assert (df.SN_model == "stromeyer").all()
    else:  # an unidentifiable fit is skipped with its reason, never replaced by Basquin
        assert "S-N model not usable" in outcomes[0].reason
    base = pd.read_csv(tmp_path / "r" / "a" / "CROSSLAND_Sn.csv")
    assert (base.SN_model == "basquin").all()
    if outcomes[0].status == "ok":
        with pytest.raises(ValueError, match="mixes master-curve models"):
            collect_all_results(tmp_path / "r")
