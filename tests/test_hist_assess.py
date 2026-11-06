"""End-to-end assessment of a history: events, damage parameters, accumulation, record."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pytest

from multifat.history.accumulate import Endurance, SNCurve
from multifat.history.assess import HistoryAssessment, assess_history
from multifat.history.eligibility import IneligibleHistoryError
from multifat.history.signal import TimeHistory
from multifat.models.registry import _PROBE_PARAMS, declared_param_names, get_models
from multifat.planes.sampling import build_R_list

E_MOD, NU = 210000.0, 0.3
R_LIST = build_R_list(n1=3, n2=6)


def _params(model):
    names = set(declared_param_names()) | set(model.required_params())
    out = {k: _PROBE_PARAMS.get(k, 1.0) for k in names}
    out.update(E=E_MOD, nu=NU)
    return out


def _hooke(S):
    trace = np.trace(S, axis1=-2, axis2=-1)[..., None, None]
    return ((1 + NU) * S - NU * trace * np.eye(3)) / E_MOD


def _block(levels, *, strain=True, mean=30.0, shear=0.4, per_segment=4):
    """Proportional tension-torsion history through the given load-factor levels."""
    lam = np.concatenate(
        [np.linspace(a, b, per_segment, endpoint=False) for a, b in zip(levels[:-1], levels[1:])]
        + [[levels[-1]]]
    )
    S = np.zeros((len(lam), 3, 3))
    S[:, 2, 2] = mean + 100.0 * lam
    S[:, 0, 2] = S[:, 2, 0] = 100.0 * shear * lam
    return TimeHistory.from_arrays(
        np.arange(len(lam), dtype=float), S, _hooke(S) if strain else None, stress_unit="MPa"
    )


def _curve(name, variant="", **kw):
    return SNCurve(A=2000.0, b=-0.12, method=name, variant=variant, **kw)


def _direct_dp(model, history, pair, variant=""):
    """Damage parameter of one event, evaluated the way a two-step constant-amplitude case is."""
    a, b = pair
    S0, S1 = history.S_hist[a], history.S_hist[b]
    E0, E1 = history.E_hist[a], history.E_hist[b]
    harmonics = {
        "S_mean": 0.5 * (S0 + S1), "S_sin": 0.5 * (S0 - S1), "S_cos": np.zeros((3, 3)),
        "E_mean": 0.5 * (E0 + E1), "E_sin": 0.5 * (E0 - E1), "E_cos": np.zeros((3, 3)),
    }  # fmt: skip
    result = model.evaluate_case(
        S0=S0, E0=E0, S1=S1, E1=E1, params=_params(model), R_list=R_LIST, harmonics=harmonics
    )
    return result.values[variant]


def _assess(history, name="FS", **kw):
    (model,) = get_models([name])
    kw.setdefault("variant", "")
    return assess_history(
        history, model, _params(model), _curve(name, kw["variant"]), R_list=R_LIST, **kw
    )


@pytest.mark.parametrize("stress_scale", [1e-12, 1.0, 1e12])
def test_strain_methods_need_one_shared_stress_strain_load_factor(stress_scale):
    S = np.zeros((5, 3, 3))
    E = np.zeros_like(S)
    S[:, 2, 2] = np.array([0, 100, 0, -100, 0]) * stress_scale
    E[:, 2, 2] = [0, 0.001, 0, 0.001, 0]
    history = TimeHistory.from_arrays(np.arange(5.0), S, E, stress_unit="MPa")
    with pytest.raises(IneligibleHistoryError, match="joint|shared|stress.strain"):
        _assess(history)
    result = _assess(history, allow_nonproportional_approximation=True)
    assert result.approximation is not None
    assert not result.settings["proportionality"]["joint"]["proportional"]
    assert _assess(history, name="FIN").approximation is None


def test_varying_strain_cannot_be_counted_from_constant_stress(tmp_path):
    S = np.zeros((5, 3, 3))
    E = np.zeros_like(S)
    E[:, 0, 2] = E[:, 2, 0] = [0, 0.001, 0, -0.001, 0]
    history = TimeHistory.from_arrays(np.arange(5.0), S, E, stress_unit="MPa")
    with pytest.raises(IneligibleHistoryError, match="joint|shared|stress.strain"):
        _assess(history)


def test_constant_amplitude_history_reproduces_the_ca_value_and_miner_sum():
    m = 5
    h = _block([-1.0, 1.0] * m + [-1.0])
    (model,) = get_models(["FS"])
    for residue in ("half", "rearrange"):
        result = _assess(h, residue=residue)
        assert isinstance(result, HistoryAssessment) and len(result.events) == 1
        (event,) = result.events
        assert event["count"] == pytest.approx(m)
        assert event["dp"] == pytest.approx(_direct_dp(model, h, tuple(event["pair"])), rel=1e-12)
        n_f = _curve("FS").cycles_to_failure(np.array([event["dp"]]))[0]
        assert result.damage_per_block == pytest.approx(m / n_f, rel=1e-12)
        assert result.blocks_to_failure == pytest.approx(1.0 / result.damage_per_block, rel=1e-12)
        assert result.cycles_per_block == pytest.approx(m)
        assert result.cycles_to_failure == pytest.approx(m * result.blocks_to_failure, rel=1e-12)


def test_variable_amplitude_history_sums_event_damages_and_each_dp_matches_direct_evaluation():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    (model,) = get_models(["FS"])
    result = _assess(h)
    assert len(result.events) >= 3
    curve = _curve("FS")
    total = 0.0
    for event in result.events:
        dp = _direct_dp(model, h, tuple(event["pair"]))
        assert event["dp"] == pytest.approx(dp, rel=1e-12)
        total += event["count"] / curve.cycles_to_failure(np.array([dp]))[0]
    assert result.damage_per_block == pytest.approx(total, rel=1e-12)


def test_residue_half_and_rearrange_are_recorded_and_differ_as_documented():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    half, rearranged = _assess(h, residue="half"), _assess(h, residue="rearrange")
    assert half.settings["residue"] == "half" and rearranged.settings["residue"] == "rearrange"
    assert half.damage_per_block != rearranged.damage_per_block
    assert (
        half.cycles_per_block != rearranged.cycles_per_block
        or half.damage_per_block != rearranged.damage_per_block
    )
    assert (
        half.history_sha256 == rearranged.history_sha256
    )  # same sampled history, different policy


def test_longest_chord_counting_gives_one_event():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    result = _assess(h, counting="longest_chord")
    assert len(result.events) == 1 and result.events[0]["count"] == 1.0
    assert result.settings["counting"] == "longest_chord"


def _elliptical(strain=False):
    t = np.linspace(0.0, 4.0 * np.pi, 160)
    S = np.zeros((len(t), 3, 3))
    S[:, 0, 2] = S[:, 2, 0] = 80.0 * np.sin(t)
    S[:, 1, 2] = S[:, 2, 1] = 40.0 * np.cos(t)
    return TimeHistory.from_arrays(t, S, _hooke(S) if strain else None, stress_unit="MPa")


def test_nonproportional_history_refused_then_allowed_with_recorded_note():
    h = _elliptical()
    with pytest.raises(IneligibleHistoryError, match="non-proportional"):
        _assess(h, name="FIN")
    result = _assess(h, name="FIN", allow_nonproportional_approximation=True)
    assert result.approximation and "approximation" in result.approximation
    assert result.settings["allow_nonproportional_approximation"] is True
    assert result.settings["proportionality"]["proportional"] is False


def test_native_life_method_refused():
    with pytest.raises(IneligibleHistoryError, match="native"):
        _assess(_block([0.0, 1.0, -1.0, 0.0]), name="MDC")


def test_strain_method_refused_without_strain():
    with pytest.raises(IneligibleHistoryError, match="strain"):
        _assess(_block([0.0, 1.0, -1.0, 0.0], strain=False), name="FS")


def test_loading_scope_error_reported_with_method_and_event():
    h = _block([0.0, 1.0, -1.0, 0.0], strain=False, mean=30.0)  # not fully reversed
    with pytest.raises(IneligibleHistoryError, match="GOLOS.*scope|scope.*GOLOS"):
        _assess(h, name="GOLOS")


def test_missing_variant_rejected():
    h = _block([0.0, 1.0, -1.0, 0.0], strain=False)
    (model,) = get_models(["CROSSLAND"])
    with pytest.raises(IneligibleHistoryError, match="variant"):
        assess_history(
            h, model, _params(model), _curve("CROSSLAND", "_ext"), R_list=R_LIST, variant="_ext"
        )


def test_curve_for_another_method_rejected():
    h = _block([0.0, 1.0, -1.0, 0.0])
    (model,) = get_models(["FS"])
    with pytest.raises(IneligibleHistoryError, match="curve"):
        assess_history(h, model, _params(model), _curve("FIN"), R_list=R_LIST)


def test_zero_damage_history_has_infinite_life_without_warnings():
    t = np.arange(6, dtype=float)
    S = np.tile(np.diag([0.0, 0.0, 50.0]), (6, 1, 1))
    h = TimeHistory.from_arrays(t, S, _hooke(S), stress_unit="MPa")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = _assess(h)
    assert result.events == [] and result.damage_per_block == 0.0
    assert np.isinf(result.blocks_to_failure) and np.isinf(result.cycles_to_failure)
    assert result.to_dict()["blocks_to_failure"] == "inf"


def test_endurance_limit_is_applied_per_event_and_recorded():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    full = _assess(h)
    dps = sorted(e["dp"] for e in full.events)
    limit = 0.5 * (dps[0] + dps[1])
    cut = _assess(h, endurance=Endurance.ignore_below(limit))
    assert cut.damage_per_block < full.damage_per_block
    assert cut.settings["endurance"] == {"mode": "ignore_below", "dp_limit": limit}
    assert [e["damage"] for e in cut.events if e["dp"] < limit] == [0.0]


def test_result_json_is_reproducible_and_carries_the_history_hash(tmp_path):
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    a, b = _assess(h), _assess(h)
    assert a.to_dict() == b.to_dict() and a.history_sha256 == h.sha256()
    path = a.to_json(tmp_path / "result.json")
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded["history_sha256"] == h.sha256() and loaded["method"] == "FS"
    assert loaded["curve"]["b"] == -0.12 and loaded["settings"]["counting"] == "rainflow"


def test_event_count_guard():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    with pytest.raises(ValueError, match="limit"):
        _assess(h, max_events=1)


# ---- review fixes ---------------------------------------------------------------------------


def _noisy(h, rel, seed=0):
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=h.S_hist.shape)
    noise = 0.5 * (noise + np.swapaxes(noise, -1, -2))
    S = h.S_hist + rel * np.max(np.abs(h.S_hist)) * noise
    return TimeHistory.from_arrays(h.time, S, None, stress_unit="MPa")


def test_proportionality_tolerance_is_a_parameter_and_recorded():
    h = _noisy(_block([0.0, 1.0, -1.0, 0.5, 0.0], strain=False), 1e-4)
    with pytest.raises(IneligibleHistoryError, match="non-proportional"):
        _assess(h, name="FIN")
    result = _assess(h, name="FIN", proportionality_tol=1e-3)
    assert result.approximation is None
    assert result.settings["proportionality"]["tol"] == 1e-3
    assert result.settings["proportionality"]["proportional"]


def test_default_residue_is_rearrange_and_half_carries_a_note():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3])  # does not close: the end differs from the start
    default, half = _assess(h), _assess(h, residue="half")
    assert default.settings["residue"] == "rearrange" and not default.notes
    assert any("consecutive blocks" in n for n in half.notes)


def test_rearranged_block_damage_equals_the_damage_of_one_more_tiled_block():
    """Oracle for a repeating block: tile it K times and count with the standard half rule; one more
    block in the periodic regime adds the rearranged block damage."""
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3])

    def tiled(k):
        n = h.n_samples
        S = np.concatenate([h.S_hist] * k)
        E = np.concatenate([h.E_hist] * k)
        return TimeHistory.from_arrays(np.arange(n * k, dtype=float), S, E, stress_unit="MPa")

    d20 = _assess(tiled(20), residue="half").damage_per_block
    d21 = _assess(tiled(21), residue="half").damage_per_block
    assert _assess(h, residue="rearrange").damage_per_block == pytest.approx(d21 - d20, rel=1e-9)


def test_event_damage_parameters_available_for_native_life_methods():
    from multifat.history.assess import event_damage_parameters

    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0], strain=False)
    (model,) = get_models(["MDC"])
    rows = event_damage_parameters(h, model, _params(model), R_list=R_LIST)
    assert rows and all(np.isfinite(r["dp"]) and r["count"] > 0.0 for r in rows)
    with pytest.raises(IneligibleHistoryError, match="native"):
        _assess(h, name="MDC")


def test_event_damage_parameters_still_refuse_missing_strain_and_nonproportional():
    from multifat.history.assess import event_damage_parameters

    (fs,) = get_models(["FS"])
    with pytest.raises(IneligibleHistoryError, match="strain"):
        event_damage_parameters(
            _block([0.0, 1.0, -1.0], strain=False), fs, _params(fs), R_list=R_LIST
        )
    (fin,) = get_models(["FIN"])
    with pytest.raises(IneligibleHistoryError, match="non-proportional"):
        event_damage_parameters(_elliptical(), fin, _params(fin), R_list=R_LIST)


def test_variant_is_validated_even_when_there_are_no_events():
    t = np.arange(6, dtype=float)
    S = np.tile(np.diag([0.0, 0.0, 50.0]), (6, 1, 1))
    h = TimeHistory.from_arrays(t, S, _hooke(S), stress_unit="MPa")
    (model,) = get_models(["CROSSLAND"])
    with pytest.raises(IneligibleHistoryError, match="variant"):
        assess_history(
            h, model, _params(model), _curve("CROSSLAND", "_ext"), R_list=R_LIST, variant="_ext"
        )


def test_life_below_one_cycle_is_noted():
    h = _block([0.0, 1.0, -1.0, 0.0], strain=False)
    (model,) = get_models(["FIN"])
    steep = SNCurve(A=1.0, b=-0.12, method="FIN")  # every damage parameter above A: N_f < 1
    result = assess_history(h, model, _params(model), steep, R_list=R_LIST)
    assert any("N_f < 1" in n for n in result.notes)


def test_nonproportional_strain_is_refused_for_strain_methods_only():
    t = np.linspace(0.0, 4.0 * np.pi, 160)
    S = np.zeros((len(t), 3, 3))
    S[:, 2, 2] = 100.0 * np.sin(t)  # proportional stress
    E = np.zeros((len(t), 3, 3))
    E[:, 2, 2] = 5e-4 * np.sin(t)
    E[:, 0, 0] = 2e-4 * np.cos(t)  # the strain path leaves the line
    h = TimeHistory.from_arrays(t, S, E, stress_unit="MPa")
    with pytest.raises(IneligibleHistoryError, match="strain"):
        _assess(h, name="FS")
    assert _assess(h, name="FIN").events  # a stress-only method never reads the strain


def test_plane_grid_hash_and_chord_settings_are_recorded():
    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0])
    result = _assess(
        h, counting="longest_chord", chord_metric="principal_spread", chord_on="strain"
    )
    assert result.settings["plane_grid"]["n_planes"] == len(R_LIST)
    assert len(result.settings["plane_grid"]["sha256"]) == 64
    assert result.settings["chord"] == {"metric": "principal_spread", "on": "strain"}


def test_every_eligible_method_assesses_a_proportional_history_end_to_end():
    from multifat.history.eligibility import method_matrix

    h = _block([0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0], mean=30.0)
    ran, refused = 0, 0
    for row in method_matrix():
        if row["accumulation"] != "yes":
            continue
        (model,) = get_models([row["method"]])
        params = _params(model)
        if (
            row["method"] == "MMP"
        ):  # its limits must satisfy Sigm1 < Sig0 < 2 Sigm1 (probe values do not)
            params.update(Sigm1=300.0, Taum1=200.0, Sig0=450.0, Tau0=300.0)
        try:
            result = assess_history(h, model, params, _curve(row["method"]), R_list=R_LIST)
        except IneligibleHistoryError as exc:
            assert row[
                "scope_restricted"
            ], f"{row['method']} refused without a scope restriction: {exc}"
            refused += 1
            continue
        assert np.isfinite(result.damage_per_block) and result.damage_per_block >= 0.0, row[
            "method"
        ]
        ran += 1
    assert ran >= 20 and refused >= 1
