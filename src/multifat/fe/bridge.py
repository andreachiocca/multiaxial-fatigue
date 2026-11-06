"""Bridge from FE load-case histories to the case-level models, with an eligibility check.

A range strategy selects a pair of load cases per location; the bridge turns the pair into the
two-load-step inputs ``S0, E0, S1, E1`` and the harmonic dictionary of the existing models, with
the convention of ``imp.experimental_data.extract_tensor``. The reduction is exact only for a
two-state path. For curved or phase-shifted histories it discards information (the harmonic
dictionary gets no quadrature component), and the maximum normal stress over all load cases does
not reach the two-step interface; extending the criteria to those histories is issue #20.

Every registered model takes the tensors of one material point, so all are eligible on local FE
data. A method that needs more (gradients, volume) declares ``requires_nonlocal = True`` and is
refused with a diagnostic; no model file is edited for this.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import HistoryBlock, HistoryMeta
from .ranges import RangeStrategy


class IneligibleMethodError(ValueError):
    """The method needs information that a local-only history does not carry."""


@dataclass(frozen=True)
class CaseInputs:
    """Two-load-step inputs of one location for ``FatigueModel.evaluate_case``."""

    location_id: int
    pair: tuple[int, int]
    S0: np.ndarray
    E0: np.ndarray
    S1: np.ndarray
    E1: np.ndarray
    harmonics: dict[str, np.ndarray]


def check_eligibility(meta: HistoryMeta, method) -> None:
    """Raise :class:`IneligibleMethodError` if ``method`` needs more than ``meta`` provides."""
    if meta.local_only and getattr(method, "requires_nonlocal", False):
        name = getattr(method, "name", type(method).__name__)
        raise IneligibleMethodError(
            f"{name} requires nonlocal data (gradients or volume), but this FE history holds "
            "only local tensors; it cannot be evaluated on it"
        )


def _require_stress_only(method) -> None:
    """Refuse a method that needs strain when the history has none.

    The registry knows which registered models run without strain (it finds out by evaluating
    them); a method it does not list is assumed to need strain. Imported here, on first use, so
    the rest of the package does not depend on the models.
    """
    from multifat.models.registry import stress_only_models

    name = getattr(method, "name", type(method).__name__)
    if name not in stress_only_models():
        raise IneligibleMethodError(
            f"{name} needs strain, but this FE history carries stress only; supply a strain "
            "history (for example from the elastic constants) or choose a stress-only method"
        )


def to_case_inputs(
    block: HistoryBlock,
    strategy: RangeStrategy,
    *,
    meta: HistoryMeta | None = None,
    method=None,
) -> list[CaseInputs]:
    """One :class:`CaseInputs` per location of ``block``, from the pair ``strategy`` selects.

    ``S0`` and ``S1`` are the raw tensors at load cases ``a`` and ``b``; ``S_mean = (S0 + S1) / 2``,
    ``S_sin = (S0 - S1) / 2``, ``S_cos = 0``, and likewise for strain. Without strain, ``E0``,
    ``E1`` and every ``E_*`` are NaN, the convention the registry probes stress-only models with.
    When ``method`` is given, ``meta`` is required to check its eligibility, and a history without
    strain refuses every method that is not stress-only (``IneligibleMethodError``).
    """
    if method is not None:
        if meta is None:
            raise ValueError("meta is required to check the eligibility of a method")
        check_eligibility(meta, method)
        if block.E_hist is None:
            _require_stress_only(method)
    nan = np.full((3, 3), np.nan)
    zero = np.zeros((3, 3))
    out = []
    for m in range(block.n_locations):
        E_hist = None if block.E_hist is None else block.E_hist[m]
        a, b = strategy(block.S_hist[m], E_hist).pair
        S0, S1 = block.S_hist[m, a].copy(), block.S_hist[m, b].copy()
        E0, E1 = (
            (nan.copy(), nan.copy()) if E_hist is None else (E_hist[a].copy(), E_hist[b].copy())
        )
        harmonics = {
            "S_mean": 0.5 * (S0 + S1),
            "S_sin": 0.5 * (S0 - S1),
            "S_cos": zero.copy(),
            "E_mean": 0.5 * (E0 + E1),
            "E_sin": 0.5 * (E0 - E1),
            "E_cos": zero.copy() if E_hist is not None else nan.copy(),
        }
        out.append(CaseInputs(block.location_id[m].item(), (a, b), S0, E0, S1, E1, harmonics))
    return out
