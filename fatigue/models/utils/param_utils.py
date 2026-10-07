from __future__ import annotations

"""Small helpers for accessing material parameters.

The project frequently works with partially-specified material property
dictionaries coming from Excel. For research-grade methods we prefer to:

1) Use the provided value if present.
2) Otherwise fall back to a reasonable estimate (when possible).
3) Otherwise keep running but emit a warning.

This file centralizes the warning style so all methods behave consistently.
"""

from typing import Mapping, Optional
import warnings


def get_param(
    params: Mapping[str, float],
    key: str,
    *,
    default: Optional[float] = None,
    warn: bool = True,
    context: str = "",
) -> Optional[float]:
    """Return params[key] if present and finite, else `default`.

    If `warn=True`, emits a runtime warning when falling back.
    """
    v = params.get(key, None)
    try:
        if v is None:
            raise KeyError
        v = float(v)
        if v != v:  # NaN
            raise ValueError
        return v
    except Exception:
        if warn:
            msg = f"Material property '{key}' was not present in the material data file"
            if context:
                msg += f" (needed by {context})"
            if default is not None:
                msg += f". Using default/estimated value: {default}."
            else:
                msg += "."
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return default
