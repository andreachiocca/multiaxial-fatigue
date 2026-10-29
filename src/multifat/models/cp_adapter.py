# fatigue_cp/models/cp_adapter.py
from __future__ import annotations
from typing import Mapping

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.base import CaseResult, CriticalPlaneMethod


class CPModelAdapter:
    """
    Wraps a CriticalPlaneMethod (plane-level plugin) into a FatigueModel (case-level model).
    Produces two variants:
      ""     -> by_metric (legacy NAME)
      "_ext" -> ext       (legacy NAME_ext)
    """

    def __init__(self, method: CriticalPlaneMethod):
        self._m = method
        self.name = method.name
        self.optional_params = getattr(method, "optional_params", {})
        if hasattr(method, "comparison_capabilities"):
            self.comparison_capabilities = method.comparison_capabilities
        # Native life law of the plane method, if any (see multifat.cli.run.write_variant).
        for attr in ("native_life", "life_fit", "expected_dp_at_life"):
            if hasattr(method, attr):
                setattr(self, attr, getattr(method, attr))

    def required_params(self) -> set[str]:
        return self._m.required_params()

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        if R_list is None:
            raise ValueError(f"{self.name}: R_list is required for CP models")

        S_sin = S_cos = E_sin = E_cos = None
        if harmonics is not None:
            S_sin = harmonics.get("S_sin")
            S_cos = harmonics.get("S_cos")
            E_sin = harmonics.get("E_sin")
            E_cos = harmonics.get("E_cos")

        best = best_cp_for_case(
            S0=S0,
            E0=E0,
            S1=S1,
            E1=E1,
            R_list=R_list,
            method=self._m,
            params=params,
            S_sin=S_sin,
            S_cos=S_cos,
            E_sin=E_sin,
            E_cos=E_cos,
        )
        if not getattr(self._m, "ext_variant", True):
            return CaseResult(values={"": float(best.by_metric)})
        return CaseResult(
            values={
                "": float(best.by_metric),
                "_ext": float(best.ext),
            }
        )
