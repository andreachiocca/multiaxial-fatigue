"""Confidence of a calibration: findings that question it, and what to try to remove them.

A calibration can run to the end and still rest on weak inputs: an estimated constant, fatigue
limits from another source than the tests, calibration points close together. A procedure reports
each such finding with its severity, what was found, why it matters and what data or choice would
remove it. The level summarises them: "high" without findings, "medium" with minor findings only,
"low" with at least one major finding.

The findings describe the calibration inputs and the fit, never the predictions of the other
tests, so that no calibration is changed in view of its validation results.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, List

MAJOR = "major"
MINOR = "minor"
LEVELS = ("high", "medium", "low")


@dataclass(frozen=True)
class Finding:
    """One reason to doubt a calibration.

    Attributes
    ----------
    code : str
        Stable identifier, e.g. "A_estimated".
    severity : str
        "major" or "minor".
    what : str
        What was found, with the numbers.
    why : str
        Why it matters for the constants or the predictions.
    advice : str
        What data or choice would remove it.
    """

    code: str
    severity: str
    what: str
    why: str
    advice: str

    def __post_init__(self):
        if self.severity not in (MAJOR, MINOR):
            raise ValueError(f"severity must be {MAJOR!r} or {MINOR!r}, got {self.severity!r}")


def level(findings: Iterable[Finding]) -> str:
    """'low' with any major finding, 'medium' with minor findings only, else 'high'."""
    severities = {f.severity for f in findings}
    if MAJOR in severities:
        return "low"
    return "medium" if severities else "high"


def summary(findings: Iterable[Finding]) -> Dict[str, Any]:
    """JSON-ready confidence: level, counts and the findings (major first, then by code)."""
    ordered = sorted(findings, key=lambda f: (f.severity != MAJOR, f.code))
    return {
        "level": level(ordered),
        "major": sum(f.severity == MAJOR for f in ordered),
        "minor": sum(f.severity == MINOR for f in ordered),
        "findings": [asdict(f) for f in ordered],
    }


def headline(confidence: Dict[str, Any]) -> str:
    """One line, e.g. 'confidence low (2 major, 1 minor: A_estimated, ...)'."""
    codes = ", ".join(f["code"] for f in confidence["findings"])
    counts = f"{confidence['major']} major, {confidence['minor']} minor"
    return f"confidence {confidence['level']} ({counts}" + (f": {codes})" if codes else ")")


def report_lines(confidence: Dict[str, Any]) -> List[str]:
    """The headline, then each finding with its reason and advice, for the console."""
    lines = [headline(confidence)]
    for f in confidence["findings"]:
        lines.append(f"[{f['severity']}] {f['what']}")
        lines.append(f"    why: {f['why']}")
        lines.append(f"    try: {f['advice']}")
    return lines
