"""Resolve user-facing fatigue-method names and explicit plane-selection variants.

The normal (unsuffixed) method name selects the method's own default plane:
FS uses maximum shear-strain amplitude; FIN uses maximum full Findley damage.
The _ext suffix requests maximisation of the entire damage parameter.
These are selection conventions, not calibration-mode overrides.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from multifat.models.references import ensure_enabled
from multifat.models.registry import _available_cp_methods, _available_direct_models


@dataclass(frozen=True)
class MethodSelection:
    """One named method and one explicitly resolved output variant."""

    base_name: str
    suffix: str = ""

    @property
    def tag(self) -> str:
        return self.base_name + self.suffix


def resolve_method_tag(tag: str) -> MethodSelection:
    """Validate a public method name and its optional _ext selection."""

    if not isinstance(tag, str) or not tag or tag != tag.strip():
        raise ValueError(f"Invalid method name: {tag!r}")
    suffix = "_ext" if tag.endswith("_ext") else ""
    base_name = tag[: -len(suffix)] if suffix else tag

    critical_plane = _available_cp_methods()
    direct = _available_direct_models()
    if base_name in critical_plane:
        if suffix and not getattr(critical_plane[base_name], "ext_variant", True):
            raise ValueError(f"{base_name} has no _ext variant")
    elif base_name in direct:
        if suffix:
            raise ValueError(f"{base_name} has no _ext variant")
    else:
        raise KeyError(f"Unknown method {tag!r}")

    ensure_enabled(base_name)
    return MethodSelection(base_name, suffix)


def select_method_variants(
    names: Iterable[str] | None = None, *, all_variants: bool = False
) -> list[MethodSelection]:
    """Choose one default variant per method unless explicitly requested otherwise.

    An explicit suffix selects only that variant. With all_variants=True, a
    bare critical-plane name expands to its default and available _ext
    variant. Duplicate requests are coalesced in their original order.
    """

    if names is None:
        requested = list(_available_cp_methods()) + list(_available_direct_models())
    else:
        requested = list(names)
        if not requested:
            raise ValueError("At least one method must be requested")

    selected: list[MethodSelection] = []
    seen: set[str] = set()
    critical_plane = _available_cp_methods()
    for name in requested:
        entry = resolve_method_tag(name)
        candidates = [entry]
        if (
            all_variants
            and not entry.suffix
            and entry.base_name in critical_plane
            and getattr(critical_plane[entry.base_name], "ext_variant", True)
        ):
            candidates.append(MethodSelection(entry.base_name, "_ext"))
        for candidate in candidates:
            if candidate.tag not in seen:
                seen.add(candidate.tag)
                selected.append(candidate)
    return selected
