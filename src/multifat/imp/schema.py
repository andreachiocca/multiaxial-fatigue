"""Workbook input contract: column rules, value parsing and validation report.

A material workbook has a ``Summary`` sheet (one header row, material values on the row whose
first cells contain "Average") and a ``Dataset`` sheet (one row per test). See
``docs/workbook-format.md`` for the full contract and the per-model parameter table.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


class WorkbookError(ValueError):
    """A workbook violates the input contract; the message names file, sheet and location."""

    def __init__(self, file: str, sheet: str, where: str, expected: str):
        self.file, self.sheet, self.where, self.expected = file, sheet, where, expected
        loc = f"{sheet}!{where}" if where else sheet
        super().__init__(f"[{file}] {loc}: {expected}")


# Summary columns read by header substring (legacy rule). Only E and nu are required at load;
# the others are required per model through ``required_params()``.
SUMMARY_FIELDS = {
    "E_MPa": ("Young's modulus", True),
    "nu": ("Poisson's ratio", True),
    "Yield_strength_MPa": ("Yield strength", False),
    "Tensile_strength_MPa": ("Tensile strength", False),
    "Tensile_fatigue_limit_MPa": ("Tensile fatigue limit", False),
    "Torsional_fatigue_limit_MPa": ("Torsional fatigue limit", False),
    "k_FS": ("k FS", False),
    "k_FI": ("k FI", False),
}

# Optional Summary constants, matched on the whole header (units suffix and '*' ignored).
# Key -> accepted header names (lower case).
SUMMARY_CONSTANTS = {
    "Sig0": (
        "sig0",
        "tensile fatigue strength r=0",
        "tensile fatigue limit r=0",
        "tensile fatigue strength r = 0",
        "tensile fatigue limit r = 0",
    ),
    "Tau0": (
        "tau0",
        "torsional fatigue strength r=0",
        "torsional fatigue limit r=0",
        "torsional fatigue strength r = 0",
        "torsional fatigue limit r = 0",
    ),
    "sigma_f": ("sigma_f", "fatigue strength coefficient"),
    "epsilon_f": ("epsilon_f", "fatigue ductility coefficient"),
    "b_fat": ("b_fat", "fatigue strength exponent"),
    "c_fat": ("c_fat", "fatigue ductility exponent"),
    "tau_f": ("tau_f", "shear fatigue strength coefficient"),
    "gamma_f": ("gamma_f", "shear fatigue ductility coefficient"),
    "b0": ("b0", "shear fatigue strength exponent"),
    "c0": ("c0", "shear fatigue ductility exponent"),
    "k_MGSE": ("k_mgse", "k mgse"),
    "a_MSWT": ("a_mswt", "a mswt"),
    "alpha_w": ("alpha_w", "alpha w"),
    "xi_ZHU": ("xi_zhu", "xi zhu"),
    "L_LI": ("l_li", "li weight"),
}

# Dataset columns: key -> tokens that must all appear in the header (legacy rule), or a tuple of
# such token lists (alternatives, first match wins). Forces and moments are never stresses.
DATASET_REQUIRED = {
    "Nf": ["number", "cycles", "failure"],
    "notched": ["notched"],
    "runout": ["run", "out"],
    "mean_rad": ["mean", "radial", "stress"],
    "amp_rad": ["amplit", "radial", "stress"],
    "mean_ax": (["mean", "axial", "stress"], ["mean", "bending", "stress"]),
    "amp_ax": (["amplit", "axial", "stress"], ["amplit", "bending", "stress"]),
    "mean_sh": ["mean", "shear", "stress"],
    "amp_sh": ["amplit", "shear", "stress"],
}
DATASET_OPTIONAL = {
    "phase_rad": ["phase", "radial", "load"],
    "phase_tor": ["phase", "torsional", "load"],
    # Revised format (October 2026 collection); a Phase AXIAL column marks the format.
    "phase_ax": ["phase", "axial", "load"],
    "freq_rad": ["freq", "radial"],
    "freq_ax": ["freq", "axial"],
    "freq_tor": ["freq", "torsional"],
    "mean_p": ["mean", "pressure"],
    "amp_p": ["amplit", "pressure"],
    "specimen": ["specimen"],
    # Decision d113 (issue #18): failure, run-out (censored life) or fatigue limit (limit state)
    "observation_type": ["observation", "type"],
}

# Values of the Observation type column (decision d113). A blank cell leaves the row undeclared:
# a run-out is then a limit state only in a workbook without failures (the rule before d113).
OBSERVATION_TYPES = ("failure", "run-out", "fatigue limit")


def workbook_files(folder) -> list:
    """The workbooks of a folder (not its subfolders), without Excel lock files (``~$*.xlsx``)."""
    from pathlib import Path

    return sorted(p for p in Path(folder).glob("*.xlsx") if not p.name.startswith("~$"))


def material_id(stem: str) -> str:
    """Material ID of a workbook: the file stem without an ordering prefix such as ``6_``."""
    return re.sub(r"^\d+_", "", stem)


def tokens_text(tokens) -> str:
    """Header rule as text, for docs and messages."""
    alts = tokens if isinstance(tokens, tuple) else (tokens,)
    return " or ".join(" + ".join(t) for t in alts)


# -----------------------------
# Value parsing (unchanged legacy rules)
# -----------------------------
def norm_text(s: Any) -> str:
    return " ".join(str(s).strip().split()).lower()


def header_key(s: Any) -> str:
    """Normalised header without a trailing units group '(...)' or footnote '*'."""
    t = norm_text(s)
    t = re.sub(r"\s*\([^()]*\)\s*$", "", t)
    return t.rstrip("* ").strip()


def is_nullish(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, float) and np.isnan(x):
        return True
    if isinstance(x, (np.floating,)):
        return bool(np.isnan(float(x)))
    if isinstance(x, str):
        t = x.strip().lower()
        return t in {"", "null", "nan", "na", "n/a", "-", "—"}
    return False


def _parse_numeric_string(s: str) -> Tuple[str, bool]:
    s = s.strip()
    if s == "":
        return "", False

    if s.startswith((">", "<", "≥", "≤")):
        s = s[1:].strip()

    is_percent = False
    if "%" in s:
        s = s.replace("%", "")
        is_percent = True

    s = s.replace(" ", "").replace(" ", "")

    if "," in s and "." in s:
        last_comma = s.rfind(",")
        last_dot = s.rfind(".")
        if last_comma > last_dot:
            s = s.replace(".", "")
            s = s.replace(",", ".")
        else:
            s = s.replace(",", "")
    elif "," in s and "." not in s:
        parts = s.split(",")
        if len(parts) == 2:
            left, right = parts[0], parts[1]
            if "e" in right.lower():
                s = left + "." + right
            else:
                if left == "0" or len(right) != 3:
                    s = left + "." + right
                else:
                    s = left + right
        else:
            s = s.replace(",", "")

    return s, is_percent


def to_float_or_none(x: Any) -> Optional[float]:
    if is_nullish(x):
        return None
    if isinstance(x, (int, float, np.integer, np.floating)) and not (
        isinstance(x, float) and np.isnan(x)
    ):
        return float(x)

    s_norm, is_percent = _parse_numeric_string(str(x))
    if s_norm == "":
        return None

    try:
        val = float(s_norm)
    except Exception:
        v = pd.to_numeric(pd.Series([s_norm]), errors="coerce").iloc[0]
        if isinstance(v, float) and np.isnan(v):
            return None
        val = float(v)

    if is_percent:
        val /= 100.0
    return float(val)


def to_float_or_zero(x: Any) -> float:
    v = to_float_or_none(x)
    return 0.0 if v is None else float(v)


def find_col(columns: List[str], must_contain) -> str:
    ncols = [(c, norm_text(c)) for c in columns]
    alts = must_contain if isinstance(must_contain, tuple) else (must_contain,)
    for tokens in alts:
        for original, normalized in ncols:
            if all(token.lower() in normalized for token in tokens):
                return original
    raise KeyError(f"Could not find column tokens {must_contain}. Available: {list(columns)}")


def is_no_value(x: Any) -> bool:
    if is_nullish(x):
        return False
    if isinstance(x, (bool, np.bool_)):
        return x is False
    if isinstance(x, (int, np.integer)):
        return int(x) == 0
    t = str(x).strip().lower()
    return t.startswith("no") or t in {"false", "0", "failure", "failed"}


def is_yes_value(x: Any) -> bool:
    if is_nullish(x):
        return False
    if isinstance(x, (bool, np.bool_)):
        return x is True
    if isinstance(x, (int, np.integer)):
        return int(x) == 1
    t = str(x).strip().lower()
    return t.startswith("yes") or t in {"true", "1", "runout", "run-out", "run out"}


def parse_observation_type(x: Any) -> str | None:
    """Declared Observation type (one of ``OBSERVATION_TYPES``), None if blank.

    Raises ValueError for any other text: the column is a declaration and is not guessed.
    """
    if is_nullish(x):
        return None
    t = norm_text(x).replace("runout", "run-out").replace("run out", "run-out")
    if t not in OBSERVATION_TYPES:
        raise ValueError(f"Observation type {x!r} is not one of {', '.join(OBSERVATION_TYPES)}")
    return t


def match_constant_columns(headers: Iterable[Any], names: dict[str, Iterable[str]]) -> dict:
    """Map each constant key to the first header column (1-based) whose header_key matches."""
    keys = [header_key(h) if h is not None else "" for h in headers]
    out: dict[str, int] = {}
    for key, aliases in names.items():
        al = {a.lower() for a in aliases}
        for c_idx, k in enumerate(keys, start=1):
            if k and k in al:
                out[key] = c_idx
                break
    return out


# -----------------------------
# Validation report
# -----------------------------
@dataclass
class ValidationReport:
    file: str
    material: str = ""
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    info: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def text(self) -> str:
        lines = [f"Workbook: {self.file}", f"Material: {self.material}"]
        for title, items in (("ERRORS", self.errors), ("WARNINGS", self.warnings)):
            lines.append(f"{title} ({len(items)})")
            lines += [f"  - {s}" for s in items]
        lines.append("INFO")
        lines += [f"  - {s}" for s in self.info]
        return "\n".join(lines) + "\n"
