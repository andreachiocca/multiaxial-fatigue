"""Optional calibration switch for the Findley shear amplitude (MCC or MCE).

Adapted from the upstream implementation (A. Chiocca, 36ff0a1, 0105b69); the configuration
and its validation are retained, while the versioned calibration store and record
compatibility are omitted. The setting chooses the plane shear rule used by the Findley
criterion: MCE denotes the conventional root-sum-square rule (``rss``), and MCC denotes the
longest chord of the shear path (``longest_chord``); both run through the conventional plane
search. The standalone harmonic evaluator (``models/cp_methods/fin_harmonic.py``) provides an
independent check of this equivalence. Production Findley remains unchanged and uses MCE.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from multifat.models.cp_methods.fin import Findley
from multifat.models.utils.shear_amplitude import ALGORITHMS, DEFINITION_VERSION


class FindleyCalibrationUnavailable(ValueError):
    """The requested load history cannot support explicit amplitude calibration."""


class FindleyMCC(Findley):
    """Findley criterion using the MCC amplitude: half of the longest plane shear-path chord."""

    np_shear_rule = "longest_chord"


def resolve_findley_algorithm(selection: str = "method_default") -> str:
    """Canonical amplitude identifier: ``mce_harmonic`` for the default and the RSS alias."""
    if not isinstance(selection, str) or selection not in (
        *ALGORITHMS,
        "rss",
        "longest_chord",
        "method_default",
    ):
        raise ValueError(f"Unknown Findley shear amplitude algorithm {selection!r}")
    if selection == "longest_chord":  # half the longest chord = MCC radius for one harmonic
        return "mcc_harmonic"
    return "mce_harmonic" if selection in {"method_default", "rss"} else selection


def amplitude_identity(selection: str) -> dict[str, str]:
    """Resolved definition together with the numerical path actually employed."""
    return {
        "algorithm": resolve_findley_algorithm(selection),
        "definition_version": DEFINITION_VERSION,
        "evaluation_path": "legacy" if selection == "method_default" else "harmonic",
    }


def findley_method(selection: str = "method_default") -> Findley:
    """Findley criterion for a calibration: production class or the MCC variant."""
    return FindleyMCC() if resolve_findley_algorithm(selection) == "mcc_harmonic" else Findley()


def require_stress_harmonics(case: Mapping[str, Any], material: str, case_id: int) -> None:
    """Validate stored coefficients before the extractor can construct a fallback from endpoints."""
    coefficient_block = case.get("harmonics")
    for tensor_name in ("S_sin", "S_cos"):
        entry = (
            coefficient_block.get(tensor_name) if isinstance(coefficient_block, Mapping) else None
        )
        matrix = np.asarray(entry)
        if (
            matrix.shape != (3, 3)
            or matrix.dtype.kind not in "biuf"
            or not np.isfinite(matrix).all()
        ):
            raise FindleyCalibrationUnavailable(
                f"{material}/FIN case {case_id}: missing or invalid {tensor_name}; "
                "explicit amplitude calibration requires stored S_sin and S_cos tensors"
            )
