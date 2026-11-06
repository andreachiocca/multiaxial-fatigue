"""Bridge from FE load-case histories to the case-level models, and the eligibility check."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.fe.bridge import CaseInputs, IneligibleMethodError, check_eligibility, to_case_inputs
from multifat.fe.model import HistoryBlock, HistoryMeta
from multifat.fe.ranges import ExplicitPair, LongestChord
from multifat.imp.experimental_data import extract_tensor
from multifat.models.registry import (
    _available_cp_methods,
    _available_direct_models,
    get_models,
    stress_only_models,
)
from multifat.planes.sampling import build_R_list

E_MOD, NU = 210000.0, 0.3
PARAMS = {"Sy": 700.0, "k_FS": 1.0, "k_FI": 0.3, "E": E_MOD, "nu": NU}


def _sym(rng, shape):
    a = rng.normal(size=shape + (3, 3)) * 100.0
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def _hooke(S):
    return (
        (1 + NU) * S - NU * np.trace(S, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
    ) / E_MOD


def _meta(**kw):
    base = dict(
        material_id="m", solver="s", solver_version="1", analysis_step="a", adapter="t",
        adapter_version="0", source_files=(), source_units={}, source_frame="cartesian",
        averaging="none", location_kind="node",
    )  # fmt: skip
    base.update(kw)
    return HistoryMeta(**base)


def _two_state_block(strain=True):
    rng = np.random.default_rng(0)
    S = _sym(rng, (1, 2))
    return HistoryBlock(np.array([7]), np.zeros((1, 3)), S, _hooke(S) if strain else None)


def test_bridge_equals_extract_tensor_convention():
    block = _two_state_block()
    (case,) = to_case_inputs(block, LongestChord())
    assert isinstance(case, CaseInputs) and case.location_id == 7 and case.pair == (0, 1)
    S, E = block.S_hist[0], block.E_hist[0]
    data = {"mat": {"cases": {0: {1: {1: {"stress": S[0], "strain": E[0]},
                                       2: {"stress": S[1], "strain": E[1]}}}}}}  # fmt: skip
    S0, E0, S1, E1, harm = extract_tensor(data, "mat", 0, include_harmonics=True)
    for got, want in ((case.S0, S0), (case.E0, E0), (case.S1, S1), (case.E1, E1)):
        assert np.array_equal(got, want)
    assert set(case.harmonics) == set(harm)
    for key in harm:
        assert np.array_equal(case.harmonics[key], harm[key]), key


def test_two_state_history_runs_existing_model_identically():
    block = _two_state_block()
    (case,) = to_case_inputs(block, ExplicitPair(0, 1))
    S, E = block.S_hist[0], block.E_hist[0]
    data = {"mat": {"cases": {0: {1: {1: {"stress": S[0], "strain": E[0]},
                                       2: {"stress": S[1], "strain": E[1]}}}}}}  # fmt: skip
    S0, E0, S1, E1, harm = extract_tensor(data, "mat", 0, include_harmonics=True)
    R_list = build_R_list(n1=3, n2=6)
    model = get_models(["FS"])[0]
    ref = model.evaluate_case(
        S0=S0, E0=E0, S1=S1, E1=E1, params=PARAMS, R_list=R_list, harmonics=harm
    )
    got = model.evaluate_case(
        S0=case.S0, E0=case.E0, S1=case.S1, E1=case.E1, params=PARAMS, R_list=R_list,
        harmonics=case.harmonics,
    )  # fmt: skip
    assert got.values == ref.values and np.isfinite(got.values[""])


def test_stress_only_model_runs_without_strain_and_strain_model_fails_clearly():
    (case,) = to_case_inputs(_two_state_block(strain=False), LongestChord())
    assert np.isnan(case.E0).all() and np.isnan(case.harmonics["E_sin"]).all()
    R_list = build_R_list(n1=3, n2=6)

    def run(name):
        return get_models([name])[0].evaluate_case(
            S0=case.S0, E0=case.E0, S1=case.S1, E1=case.E1, params=PARAMS, R_list=R_list,
            harmonics=case.harmonics,
        )  # fmt: skip

    assert "FIN" in stress_only_models() and np.isfinite(run("FIN").values[""])
    assert "FS" not in stress_only_models()
    with pytest.raises(ValueError, match="nonfinite"):
        run("FS")


class _NonlocalStub:
    name = "STUB"
    requires_nonlocal = True


def test_stub_requires_nonlocal_refused_with_reason():
    block = _two_state_block()
    with pytest.raises(IneligibleMethodError, match="nonlocal"):
        to_case_inputs(block, LongestChord(), meta=_meta(), method=_NonlocalStub())
    with pytest.raises(ValueError, match="meta"):
        to_case_inputs(block, LongestChord(), method=_NonlocalStub())
    # a history that is not local-only would accept it
    check_eligibility(_meta(local_only=False), _NonlocalStub())


def test_all_registered_methods_pass_eligibility():
    names = list(_available_cp_methods()) + list(_available_direct_models())
    checked = 0
    for name in names:
        try:
            (model,) = get_models([name])
        except NotImplementedError:  # disabled by the numerical audit
            continue
        check_eligibility(_meta(), model)
        checked += 1
    assert checked >= 30


def test_nonproportional_reduction_is_documented_lossy():
    t = np.arange(360) * (2.0 * np.pi / 360)
    S = np.zeros((1, 360, 3, 3))
    S[0, :, 0, 2] = S[0, :, 2, 0] = 100.0 * np.sin(t)  # xz sine
    S[0, :, 1, 2] = S[0, :, 2, 1] = 100.0 * np.cos(t)  # yz cosine: a circular shear path
    block = HistoryBlock(np.array([1]), np.zeros((1, 3)), S, None)
    (case,) = to_case_inputs(block, LongestChord())
    # the two-state reduction is a straight line through the circle: no quadrature component
    assert not np.any(case.harmonics["S_cos"])
    # while the history itself has one (the true cosine coefficient of the yz shear)
    true_cos = 2.0 * np.mean(S[0, :, 1, 2] * np.cos(t))
    assert abs(true_cos) == pytest.approx(100.0, rel=1e-9)


@pytest.mark.parametrize("name", ["FS", "EK", "SWT", "MKBM"])
def test_strain_methods_refused_when_history_has_no_strain(name):
    block = _two_state_block(strain=False)
    (model,) = get_models([name])
    assert name not in stress_only_models()
    with pytest.raises(IneligibleMethodError, match="strain"):
        to_case_inputs(block, LongestChord(), meta=_meta(), method=model)


def test_stress_only_methods_accepted_without_strain_and_all_accepted_with_it():
    (fin,) = get_models(["FIN"])
    assert (
        len(
            to_case_inputs(_two_state_block(strain=False), LongestChord(), meta=_meta(), method=fin)
        )
        == 1
    )
    (fs,) = get_models(["FS"])
    assert (
        len(to_case_inputs(_two_state_block(strain=True), LongestChord(), meta=_meta(), method=fs))
        == 1
    )
