"""Eligibility of methods and histories for the history bridge."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.history.accumulate import SNCurve
from multifat.history.eligibility import (
    IneligibleHistoryError,
    check_accumulation_eligibility,
    check_history_kind,
    is_native_life,
    main,
    method_matrix,
)
from multifat.history.proportionality import Proportionality, proportionality
from multifat.models.references import DISABLED_METHODS
from multifat.models.registry import (
    _available_cp_methods,
    _available_direct_models,
    get_models,
    stress_only_models,
)


def _curve(name, variant=""):
    return SNCurve(A=1000.0, b=-0.1, method=name, variant=variant)


def _enabled():
    return [
        n
        for n in list(_available_cp_methods()) + list(_available_direct_models())
        if n not in DISABLED_METHODS
    ]


def test_matrix_has_a_row_for_every_enabled_registered_method():
    rows = method_matrix()
    assert sorted(r["method"] for r in rows) == sorted(_enabled())
    assert {
        "method",
        "strain_needed",
        "native_life",
        "requires_nonlocal",
        "accumulation",
        "reason",
    } <= set(rows[0])


def test_native_life_methods_refused_for_accumulation():
    native = [r["method"] for r in method_matrix() if r["native_life"]]
    assert {"MDC", "MMCCM", "MOREL"} <= set(native)
    for name in native:
        (model,) = get_models([name])
        assert is_native_life(model)
        with pytest.raises(IneligibleHistoryError, match="native"):
            check_accumulation_eligibility(model, _curve(name), variant="", has_strain=True)
    row = next(r for r in method_matrix() if r["method"] == "MDC")
    assert row["accumulation"] == "refused" and "native" in row["reason"]


def test_every_non_native_method_with_strain_is_eligible():
    for row in method_matrix():
        if row["native_life"]:
            continue
        (model,) = get_models([row["method"]])
        check_accumulation_eligibility(model, _curve(row["method"]), variant="", has_strain=True)
        assert row["accumulation"] == "yes"


def test_strain_method_refused_without_strain_and_accepted_with_it():
    (fs,) = get_models(["FS"])
    assert "FS" not in stress_only_models()
    with pytest.raises(IneligibleHistoryError, match="strain"):
        check_accumulation_eligibility(fs, _curve("FS"), variant="", has_strain=False)
    check_accumulation_eligibility(fs, _curve("FS"), variant="", has_strain=True)
    (fin,) = get_models(["FIN"])
    assert "FIN" in stress_only_models()
    check_accumulation_eligibility(fin, _curve("FIN"), variant="", has_strain=False)


def test_curve_method_and_variant_mismatch_rejected():
    (fs,) = get_models(["FS"])
    with pytest.raises(IneligibleHistoryError, match="curve"):
        check_accumulation_eligibility(fs, _curve("FIN"), variant="", has_strain=True)
    with pytest.raises(IneligibleHistoryError, match="variant"):
        check_accumulation_eligibility(fs, _curve("FS", "_ext"), variant="", has_strain=True)
    check_accumulation_eligibility(fs, _curve("FS", "_ext"), variant="_ext", has_strain=True)


class _NonlocalStub:
    name = "STUB"
    requires_nonlocal = True


def test_stub_requires_nonlocal_refused():
    with pytest.raises(IneligibleHistoryError, match="nonlocal"):
        check_accumulation_eligibility(_NonlocalStub(), _curve("STUB"), variant="", has_strain=True)


def test_nonproportional_refused_unless_allowed_and_note_returned():
    s = np.zeros((200, 3, 3))
    t = np.linspace(0.0, 2.0 * np.pi, 200)
    s[:, 0, 2] = s[:, 2, 0] = 80.0 * np.sin(t)
    s[:, 1, 2] = s[:, 2, 1] = 40.0 * np.cos(t)
    verdict = proportionality(s)
    assert isinstance(verdict, Proportionality) and not verdict.proportional
    with pytest.raises(IneligibleHistoryError, match="non-proportional"):
        check_history_kind(verdict, allow_approximation=False)
    note = check_history_kind(verdict, allow_approximation=True)
    assert note and "approximation" in note
    assert check_history_kind(Proportionality(0.0, True, 1e-6), allow_approximation=False) is None


def test_markdown_table_lists_every_method(capsys):
    assert main([]) == 0
    out = capsys.readouterr().out
    for name in _enabled():
        assert f"`{name}`" in out
    assert out.lstrip().startswith("| method")


def test_matrix_flags_loading_scope_restricted_methods():
    rows = {r["method"]: r for r in method_matrix()}
    assert rows["GOLOS"]["scope_restricted"] is True
    assert rows["FS"]["scope_restricted"] is False and rows["FIN"]["scope_restricted"] is False
