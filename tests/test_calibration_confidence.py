"""Confidence of a calibration: the level summarises the findings, major first."""

import pytest

from multifat.eval.calibration_confidence import (
    Finding,
    headline,
    level,
    report_lines,
    summary,
)


def _finding(code, severity):
    return Finding(code, severity, f"{code} found", "it matters", "do this")


def test_level_follows_the_worst_finding():
    assert level([]) == "high"
    assert level([_finding("b", "minor")]) == "medium"
    assert level([_finding("b", "minor"), _finding("a", "major")]) == "low"


def test_summary_headline_and_report():
    result = summary([_finding("z_minor", "minor"), _finding("y_major", "major")])
    assert (result["level"], result["major"], result["minor"]) == ("low", 1, 1)
    assert [f["code"] for f in result["findings"]] == ["y_major", "z_minor"]
    assert headline(result) == "confidence low (1 major, 1 minor: y_major, z_minor)"
    assert report_lines(result)[1:4] == [
        "[major] y_major found",
        "    why: it matters",
        "    try: do this",
    ]
    assert headline(summary([])) == "confidence high (0 major, 0 minor)"


def test_severity_is_checked():
    with pytest.raises(ValueError, match="severity"):
        _finding("a", "serious")
