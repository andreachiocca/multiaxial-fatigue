"""Correspondence between multifat nomenclature and FatPy nomenclature (v01).

Implements parameter mapping and unit conversions between multifat (MPa-based)
and FatPy (Young's modulus in GPa, stresses in MPa) as defined in docs/notation.md.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Any, Mapping

MAPPING_FILE = Path(__file__).resolve().parent / "fatpy_keys.json"


@functools.lru_cache(maxsize=1)
def load_fatpy_mapping() -> dict[str, dict[str, Any]]:
    """Load the JSON correspondence table between multifat and FatPy keys."""
    return json.loads(MAPPING_FILE.read_text(encoding="utf-8"))


def to_fatpy(params: Mapping[str, Any]) -> dict[str, Any]:
    """Convert a dictionary of multifat parameters to FatPy nomenclature.

    Unmapped or ambiguous quantities (such as Sig0 and Tau0, pending clarification)
    are excluded from the output. Young's modulus E is converted from MPa to GPa.

    Parameters
    ----------
    params:
        Mapping of multifat parameter names to their numeric values.

    Returns
    -------
    dict:
        Dictionary keyed by FatPy parameter names with converted values.
    """
    mapping = load_fatpy_mapping()
    out: dict[str, Any] = {}
    for key, value in params.items():
        entry = mapping.get(key)
        if entry is None:
            continue
        fatpy_key = entry.get("fatpy_key")
        factor = entry.get("factor_to_fatpy")
        if fatpy_key is not None and factor is not None:
            if value is None:
                out[fatpy_key] = None
            else:
                out[fatpy_key] = float(value) * factor
    return out


def from_fatpy(fatpy_params: Mapping[str, Any]) -> dict[str, Any]:
    """Convert a dictionary of FatPy parameters to multifat nomenclature.

    Inverts FatPy naming and unit conventions (converting Young's modulus from GPa
    to MPa). Supports recognized FatPy alternatives such as ``fat_lim_ben_m1`` for Sigm1.

    Parameters
    ----------
    fatpy_params:
        Mapping of FatPy parameter names to their numeric values.

    Returns
    -------
    dict:
        Dictionary keyed by multifat parameter names with converted values.
    """
    mapping = load_fatpy_mapping()
    # Reverse lookup: fatpy_key -> (multifat_name, factor_from_fatpy)
    reverse: dict[str, tuple[str, float]] = {}
    for name, entry in mapping.items():
        fkey = entry.get("fatpy_key")
        factor = entry.get("factor_to_fatpy")
        if fkey is not None and factor is not None and factor != 0:
            inv = 1.0 / factor
            reverse[fkey] = (name, inv)
            for alt in entry.get("fatpy_alternatives", []):
                reverse[alt] = (name, inv)

    out: dict[str, Any] = {}
    for fkey, value in fatpy_params.items():
        if fkey in reverse:
            mname, inv_factor = reverse[fkey]
            if value is None:
                out[mname] = None
            else:
                out[mname] = float(value) * inv_factor
    return out
