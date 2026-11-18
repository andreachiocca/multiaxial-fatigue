"""User selection of the plane shear amplitude (longest chord / MCC / MCE) for every model."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.eval.shear_amplitude_option import (
    SELECTIONS,
    apply_shear_amplitude,
    plane_rule,
    resolved_algorithm,
)
from multifat.models.registry import available_cp_methods, get_models


def test_selections_map_to_two_plane_rules():
    assert plane_rule("method_default") is None
    assert {plane_rule(s) for s in SELECTIONS[1:]} == {"longest_chord", "rss"}
    with pytest.raises(ValueError):
        plane_rule("mcc")


@pytest.mark.parametrize("name", sorted(available_cp_methods()))
def test_every_critical_plane_method_is_selectable(name):
    for selection, expected in (
        ("longest_chord", "mcc_harmonic"),
        ("mce_harmonic", "mce_harmonic"),
    ):
        model = get_models([name])[0]
        assert apply_shear_amplitude(model, selection)
        assert resolved_algorithm(model) == expected


def test_selection_does_not_leak_into_new_instances():
    model = get_models(["FIN"])[0]
    apply_shear_amplitude(model, "longest_chord")
    assert resolved_algorithm(get_models(["FIN"])[0]) == "mce_harmonic"  # Findley's own rule


def test_sih_follows_the_selection_and_fixed_models_are_reported():
    assert resolved_algorithm(get_models(["SIH"])[0]) == "mcc_harmonic"
    sih = get_models(["SIH"])[0]
    assert apply_shear_amplitude(sih, "rss") and resolved_algorithm(sih) == "mce_harmonic"
    for fixed in ("BP", "CROSSLAND", "WANG2018"):
        assert not apply_shear_amplitude(get_models([fixed])[0], "rss")


def test_the_selection_changes_an_out_of_phase_findley_value_only():
    """90 deg tension-torsion: MCC and MCE differ out of phase, coincide in phase."""
    from multifat.cli.run import build_R_list

    R_list = build_R_list(n1=6, n2=12)
    params = {"k_FI": 0.3}
    values = {}
    for delta in (0.0, 90.0):
        d = np.radians(delta)
        S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3))
        S_sin[1, 1] = 200.0
        S_sin[0, 1] = S_sin[1, 0] = 100.0 * np.cos(d)
        S_cos[0, 1] = S_cos[1, 0] = 100.0 * np.sin(d)
        h = {"S_mean": np.zeros((3, 3)), "S_sin": S_sin, "S_cos": S_cos}
        for selection in ("longest_chord", "mce_harmonic"):
            model = get_models(["FIN"])[0]
            apply_shear_amplitude(model, selection)
            S0, S1 = S_sin, -S_sin
            values[delta, selection] = model.evaluate_case(
                S0=S0, E0=S0 * 0, S1=S1, E1=S1 * 0, params=params, R_list=R_list, harmonics=h
            ).values[""]
    assert values[0.0, "longest_chord"] == pytest.approx(values[0.0, "mce_harmonic"], rel=1e-9)
    assert values[90.0, "longest_chord"] < values[90.0, "mce_harmonic"]
