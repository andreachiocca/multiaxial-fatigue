"""User selection of the plane shear amplitude for every model that exposes one.

A model exposes its plane shear amplitude through ``np_shear_rule``: every critical-plane method
(through the plane preparation of ``eval/evaluator.py``) and the case-level models that declare
the attribute (e.g. SIH). The selection replaces the method's own rule for one run:

- ``method_default``: the rule of the method's source (or the repository choice it documents);
- ``longest_chord`` and ``mcc_harmonic``: half the longest chord of the shear path, which for a
  single-harmonic ellipse equals the minimum circumscribed circle (MCC) radius, its major semiaxis;
- ``mce_harmonic`` and ``rss``: the minimum circumscribed ellipse (MCE) amplitude, the quadratic
  combination sqrt(a^2 + b^2) of the two semiaxes.

For one harmonic the longest chord and the MCC coincide, so the five selections give at most two
numerical results per method. Models whose shear amplitude is fixed by their definition (for
example BP, RFB, PAPADOPOULOS, LIUMAHADEVAN) or that use none do not expose ``np_shear_rule`` and
are left unchanged.
"""

from __future__ import annotations

SELECTIONS = ("method_default", "longest_chord", "mcc_harmonic", "mce_harmonic", "rss")
_PLANE_RULE = {
    "longest_chord": "longest_chord",
    "mcc_harmonic": "longest_chord",
    "mce_harmonic": "rss",
    "rss": "rss",
}


def plane_rule(selection: str) -> str | None:
    """The ``np_shear_rule`` a selection imposes, or None for ``method_default``."""
    if selection not in SELECTIONS:
        raise ValueError(f"Unknown shear amplitude selection {selection!r}; choose {SELECTIONS}")
    return _PLANE_RULE.get(selection)


def amplitude_target(model):
    """The object whose ``np_shear_rule`` sets the model's plane shear amplitude, or None."""
    plane_method = getattr(model, "_m", None)  # CPModelAdapter wraps the plane-level method
    if plane_method is not None:
        return plane_method
    if hasattr(model, "evaluate_on_planes") or hasattr(model, "np_shear_rule"):
        return model  # a plane-level method, or a case-level model declaring its rule
    return None


def apply_shear_amplitude(model, selection: str) -> bool:
    """Impose a selection on a freshly built model; False when the model exposes no amplitude."""
    rule = plane_rule(selection)
    target = amplitude_target(model)
    if target is None:
        return False
    if rule is not None:
        target.np_shear_rule = rule  # instance attribute: the class default is untouched
    return True


def resolved_algorithm(model) -> str | None:
    """Canonical amplitude a built model evaluates: mcc_harmonic, mce_harmonic, or None."""
    target = amplitude_target(model)
    if target is None:
        return None
    rule = getattr(target, "np_shear_rule", "rss")  # the plane preparation's default
    return "mcc_harmonic" if rule == "longest_chord" else "mce_harmonic"
