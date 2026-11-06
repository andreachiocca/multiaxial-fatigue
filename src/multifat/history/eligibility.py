"""Eligibility of methods and histories for the history bridge, with explicit diagnostics.

Nothing is skipped silently: a method or history that the bridge cannot treat on a justified
formulation raises ``IneligibleHistoryError`` naming the reason. ``method_matrix`` lists every enabled
registered method with what it needs and whether accumulation through a Basquin curve of its damage
parameter is allowed; ``python -m multifat.history.eligibility`` prints it as a markdown table.
"""

from __future__ import annotations

import functools
from typing import Sequence

import numpy as np

from multifat.models.errors import LoadingScopeError
from multifat.models.references import DISABLED_METHODS
from multifat.models.registry import (
    _PROBE_PARAMS,
    _available_cp_methods,
    _available_direct_models,
    declared_param_names,
    get_models,
    stress_only_models,
)
from multifat.planes.sampling import build_R_list

from .accumulate import SNCurve
from .proportionality import Proportionality


class IneligibleHistoryError(ValueError):
    """A method, curve or history is outside what the history bridge can treat."""


def is_native_life(model) -> bool:
    """True if the model defines its own life law (the test ``cli/run.py`` applies)."""
    return getattr(model, "native_life", None) is not None


def check_event_eligibility(model, *, has_strain: bool) -> None:
    """Raise ``IneligibleHistoryError`` unless ``model`` can give damage parameters for a history.

    Refused: a model that declares ``requires_nonlocal`` (the history is local) and a strain-needing
    model when the history has no strain. Native-life models pass: their damage parameters per event
    are available even though they cannot accumulate through a Basquin curve.
    """
    name = getattr(model, "name", type(model).__name__)
    if getattr(model, "requires_nonlocal", False):
        raise IneligibleHistoryError(
            f"{name} requires nonlocal data (gradients or volume); a history of one point is local"
        )
    if not has_strain and name not in stress_only_models():
        raise IneligibleHistoryError(
            f"{name} needs strain, but this history carries stress only; supply a strain history "
            "or choose a stress-only method"
        )


def check_accumulation_eligibility(
    model, curve: SNCurve, *, variant: str, has_strain: bool
) -> None:
    """Raise ``IneligibleHistoryError`` unless ``model`` can accumulate damage on ``curve``.

    Refused: a curve of another method or variant; a native-life model (a single Basquin curve of its
    damage parameter is not its formulation); and what ``check_event_eligibility`` refuses.
    """
    name = getattr(model, "name", type(model).__name__)
    if curve.method != name:
        raise IneligibleHistoryError(f"the curve is for method {curve.method!r}, not {name!r}")
    if curve.variant != variant:
        raise IneligibleHistoryError(
            f"the curve is for variant {curve.variant!r} but variant {variant!r} was requested"
        )
    if is_native_life(model):
        raise IneligibleHistoryError(
            f"{name} defines its own native life law; accumulation through a single Basquin curve "
            "of its damage parameter is not its formulation (damage parameters per event remain available)"
        )
    check_event_eligibility(model, has_strain=has_strain)


def check_history_kind(
    verdict: Proportionality, *, allow_approximation: bool, what: str = "stress"
) -> str | None:
    """Raise unless the history is proportional or the approximation is explicitly allowed.

    ``what`` names the path that was measured (``"stress"`` or ``"strain"``). Returns the
    approximation note to record when a non-proportional history is allowed, else ``None``.
    """
    if verdict.proportional:
        return None
    detail = f"{what} path rank ratio {verdict.rank_ratio:.3g} > tolerance {verdict.tol:.3g}"
    if not allow_approximation:
        raise IneligibleHistoryError(
            f"non-proportional {what} history ({detail}): scalar counting with two-state events is "
            "exact only for proportional loading; pass allow_nonproportional_approximation=True to "
            "accept an approximation (it is recorded in the result), or proportionality_tol to "
            "loosen the numerical tolerance"
        )
    return (
        f"non-proportional {what} history ({detail}): events are two-state paths between counted "
        "reversals, an approximation"
    )


def _enabled_names() -> list[str]:
    names = list(_available_cp_methods()) + list(_available_direct_models())
    return [n for n in names if n.removesuffix("_ext") not in DISABLED_METHODS]


@functools.lru_cache(maxsize=1)
def scope_restricted_methods() -> frozenset[str]:
    """Methods that reject a generic load with a mean stress through ``LoadingScopeError``.

    Found by evaluation, like ``stress_only_models``: each model sees one non-reversed
    tension-torsion case; a criterion defined only for a restricted loading (for example GOLOS, fully
    reversed in-phase plane stress) raises ``LoadingScopeError``, and any other outcome means the
    probe load is inside its scope. The restriction is then a property to expect on general events.
    """
    import warnings

    nu, modulus = 0.3, 210000.0

    def hooke(S):
        return ((1 + nu) * S - nu * np.trace(S) * np.eye(3)) / modulus

    S0 = np.array([[0.0, 0.0, 40.0], [0.0, 0.0, 0.0], [40.0, 0.0, 200.0]])
    S1 = np.array([[0.0, 0.0, -10.0], [0.0, 0.0, 0.0], [-10.0, 0.0, -90.0]])
    E0, E1 = hooke(S0), hooke(S1)
    zero = np.zeros((3, 3))
    harmonics = {
        "S_mean": 0.5 * (S0 + S1),
        "S_sin": 0.5 * (S0 - S1),
        "S_cos": zero,
        "E_mean": 0.5 * (E0 + E1),
        "E_sin": 0.5 * (E0 - E1),
        "E_cos": zero,
    }
    R_list = build_R_list(n1=3, n2=6)
    out = set()
    for name in _enabled_names():
        (model,) = get_models([name])
        params = {k: _PROBE_PARAMS.get(k, 1.0) for k in declared_param_names()}
        params.update(E=modulus, nu=nu)
        try:
            with warnings.catch_warnings(), np.errstate(all="ignore"):
                warnings.simplefilter("ignore")
                model.evaluate_case(
                    S0=S0, E0=E0, S1=S1, E1=E1, params=params, R_list=R_list, harmonics=harmonics
                )
        except LoadingScopeError:
            out.add(name)
        except Exception:
            pass
    return frozenset(out)


def method_matrix() -> list[dict]:
    """One row per enabled registered method.

    Keys: ``method``, ``strain_needed``, ``native_life``, ``requires_nonlocal``, ``scope_restricted``
    (the criterion rejects general loads, see ``scope_restricted_methods``), ``accumulation``
    (``"yes"`` or ``"refused"``) and ``reason``.
    """
    stress_only = stress_only_models()
    restricted = scope_restricted_methods()
    rows = []
    for name in _enabled_names():
        (model,) = get_models([name])
        native = is_native_life(model)
        nonlocal_ = bool(getattr(model, "requires_nonlocal", False))
        reasons = []
        if native:
            reasons.append("native life law")
        if nonlocal_:
            reasons.append("needs nonlocal data")
        rows.append(
            {
                "method": name,
                "strain_needed": name not in stress_only,
                "native_life": native,
                "requires_nonlocal": nonlocal_,
                "scope_restricted": name in restricted,
                "accumulation": "refused" if reasons else "yes",
                "reason": ", ".join(reasons),
            }
        )
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    """Print the method matrix as a markdown table."""
    lines = [
        "| method | strain needed | native life | scope restricted | accumulation | reason |",
        "|---|---|---|---|---|---|",
    ]
    for r in method_matrix():
        lines.append(
            f"| `{r['method']}` | {'yes' if r['strain_needed'] else 'no'} | "
            f"{'yes' if r['native_life'] else 'no'} | {'yes' if r['scope_restricted'] else 'no'} | "
            f"{r['accumulation']} | {r['reason']} |"
        )
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
