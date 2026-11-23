"""Method taxonomy of the registered models (issue #61): approach, quantity, non-proportionality.

The approach is the step that reduces the multiaxial history to a scalar, after the review
article's families: ``invariant`` (i), ``critical_plane`` (ii), ``integral`` (iii),
``incremental`` (iv) and ``path`` (v). Hybrids take the class of that final reduction.
``quantity`` is what the damage parameter is built from (stress, strain or energy).
``np_sensitivity`` is ``explicit`` when a term or a calibrated curve acts only on non-proportional
loading, ``emergent`` otherwise. Mesoscopic arguments and native life laws are tags.

The catalogue (``docs/catalogue/methods.json``) holds the same fields for every record, with the
reason of each non-obvious case; ``tests/test_taxonomy.py`` keeps the two in agreement.
"""

from __future__ import annotations

from typing import NamedTuple

APPROACHES = ("invariant", "critical_plane", "integral", "incremental", "path")
QUANTITIES = ("stress", "strain", "energy")
APPROACH_LABELS = {
    "invariant": "Invariant",
    "critical_plane": "Critical plane",
    "integral": "Integral",
    "incremental": "Incremental",
    "path": "Load path",
}


class Taxonomy(NamedTuple):
    """Classification of one registered method."""

    approach: str
    quantity: str
    np_sensitivity: str
    tags: tuple[str, ...]


# fmt: off
TAXONOMY: dict[str, Taxonomy] = {
    "BOLCHOUN": Taxonomy("integral", "stress", "explicit", ('native_life',)),
    "BP": Taxonomy("integral", "stress", "emergent", ()),
    "CAIM": Taxonomy("critical_plane", "stress", "emergent", ()),
    "COIN": Taxonomy("invariant", "stress", "emergent", ()),
    "CROSSLAND": Taxonomy("invariant", "stress", "emergent", ()),
    "DANGVAN": Taxonomy("critical_plane", "stress", "emergent", ('mesoscopic',)),
    "EK": Taxonomy("invariant", "energy", "emergent", ()),
    "FIN": Taxonomy("critical_plane", "stress", "emergent", ()),
    "FS": Taxonomy("critical_plane", "strain", "emergent", ()),
    "GAM": Taxonomy("invariant", "stress", "emergent", ()),
    "GOLOS": Taxonomy("invariant", "stress", "emergent", ()),
    "GSA": Taxonomy("critical_plane", "strain", "emergent", ()),
    "GSE": Taxonomy("critical_plane", "energy", "emergent", ()),
    "KAKUNO": Taxonomy("invariant", "stress", "emergent", ()),
    "KHALIJ": Taxonomy("invariant", "stress", "emergent", ()),
    "KUJAWSKI2021": Taxonomy("critical_plane", "energy", "emergent", ()),
    "LI": Taxonomy("critical_plane", "energy", "emergent", ()),
    "LIU2021": Taxonomy("critical_plane", "strain", "emergent", ()),
    "LIU2025": Taxonomy("critical_plane", "stress", "emergent", ()),
    "LIUMAHADEVAN": Taxonomy("critical_plane", "stress", "emergent", ('native_life',)),
    "MATAKE": Taxonomy("critical_plane", "stress", "emergent", ()),
    "MATSUBARA": Taxonomy("invariant", "stress", "emergent", ()),
    "MCD": Taxonomy("critical_plane", "stress", "emergent", ()),
    "MDC": Taxonomy("critical_plane", "stress", "emergent", ('native_life',)),
    "MGSE_YU": Taxonomy("critical_plane", "energy", "emergent", ()),
    "MGSE_ZHU": Taxonomy("critical_plane", "energy", "emergent", ()),
    "MKBM": Taxonomy("critical_plane", "strain", "emergent", ()),
    "MMCCM": Taxonomy("critical_plane", "strain", "emergent", ('native_life',)),
    "MMMH": Taxonomy("invariant", "stress", "emergent", ()),
    "MMP": Taxonomy("invariant", "stress", "emergent", ()),
    "MOREL": Taxonomy("critical_plane", "stress", "emergent", ('mesoscopic', 'native_life')),
    "MSWT": Taxonomy("critical_plane", "energy", "emergent", ()),
    "MWCM": Taxonomy("critical_plane", "stress", "emergent", ()),
    "MWHS": Taxonomy("critical_plane", "strain", "emergent", ()),
    "NINIC": Taxonomy("critical_plane", "stress", "emergent", ()),
    "OTT": Taxonomy("invariant", "stress", "emergent", ()),
    "PAPADOPOULOS": Taxonomy("critical_plane", "stress", "emergent", ('mesoscopic',)),
    "PBP": Taxonomy("invariant", "stress", "emergent", ()),
    "PCR": Taxonomy("critical_plane", "stress", "emergent", ()),
    "PDS_FS": Taxonomy("integral", "strain", "emergent", ()),
    "PDS_SWT": Taxonomy("integral", "energy", "emergent", ()),
    "PI": Taxonomy("integral", "stress", "emergent", ()),
    "RFB": Taxonomy("integral", "stress", "emergent", ()),
    "SIH": Taxonomy("integral", "stress", "emergent", ()),
    "SINES": Taxonomy("invariant", "stress", "emergent", ()),
    "SSF": Taxonomy("invariant", "stress", "emergent", ()),
    "SWT": Taxonomy("critical_plane", "energy", "emergent", ()),
    "SWTD": Taxonomy("invariant", "energy", "emergent", ()),
    "TVEIT2025": Taxonomy("incremental", "stress", "emergent", ('native_life',)),
    "VHN": Taxonomy("invariant", "stress", "emergent", ()),
    "VSE": Taxonomy("critical_plane", "energy", "emergent", ()),
    "WANG2018": Taxonomy("critical_plane", "stress", "emergent", ()),
    "WANG2022": Taxonomy("critical_plane", "stress", "emergent", ()),
    "ZHANGYAO": Taxonomy("critical_plane", "stress", "emergent", ()),
}
# fmt: on


def taxonomy(method: str) -> Taxonomy:
    """Classification of a registered method; the ``_ext`` variant shares its base's entry."""
    return TAXONOMY[method.removesuffix("_ext")]


def approach_label(method: str) -> str:
    """Display name of the method's approach, or ``"Unknown"`` for an unregistered name."""
    base = method.removesuffix("_ext")
    return APPROACH_LABELS[TAXONOMY[base].approach] if base in TAXONOMY else "Unknown"
