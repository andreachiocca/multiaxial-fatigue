"""Tests for FatPy correspondence, round-tripping, and notation integrity."""

from __future__ import annotations

import ast
from pathlib import Path
import pytest

from multifat.data.fatpy import from_fatpy, load_fatpy_mapping, to_fatpy

DOCS_NOTATION = Path(__file__).resolve().parents[1] / "docs" / "notation.md"
REPO_ROOT = Path(__file__).resolve().parents[1]


def test_fatpy_mapping_file_loads():
    mapping = load_fatpy_mapping()
    assert isinstance(mapping, dict)
    assert "E" in mapping
    assert mapping["E"]["fatpy_key"] == "elastic_modulus"
    assert mapping["E"]["fatpy_unit"] == "GPa"
    assert mapping["E"]["factor_to_fatpy"] == 0.001


def test_round_trip_mapped_parameters():
    params = {
        "E": 206000.0,
        "nu": 0.29,
        "nu_p": 0.5,
        "Sy": 680.0,
        "Su": 850.0,
        "Sigm1": 390.0,
        "Taum1": 240.0,
        "sigma_f": 1150.0,
        "b_fat": -0.085,
        "epsilon_f": 0.45,
        "c_fat": -0.58,
        "tau_f": 670.0,
        "b0": -0.085,
        "gamma_f": 0.85,
        "c0": -0.58,
        "K_prime": 1180.0,
        "n_prime": 0.14,
        "Nf": 150000.0,
        "sigma_a": 320.0,
        "sigma_m": 50.0,
        "phase_shift": 45.0,
    }
    fatpy = to_fatpy(params)
    assert fatpy["elastic_modulus"] == pytest.approx(206.0)
    assert fatpy["yield_strength"] == pytest.approx(680.0)
    assert fatpy["fat_lim_ten_m1"] == pytest.approx(390.0)

    recovered = from_fatpy(fatpy)
    for key, val in params.items():
        assert recovered[key] == pytest.approx(val), f"Mismatch for {key}"


def test_alternative_fatpy_keys():
    # fat_lim_ben_m1 is an accepted bending alternative for Sigm1
    fatpy = {"fat_lim_ben_m1": 420.0}
    recovered = from_fatpy(fatpy)
    assert recovered["Sigm1"] == pytest.approx(420.0)


def test_ambiguous_and_unmapped_keys_excluded():
    # Sig0 and Tau0 are ambiguous (max vs amplitude in FatPy v01); method constants have no keys
    params = {"Sig0": 300.0, "Tau0": 200.0, "Tu": 500.0, "k_FS": 0.5}
    fatpy = to_fatpy(params)
    assert fatpy == {}


def test_nullable_values_preserved():
    params = {"E": None, "nu": 0.3}
    fatpy = to_fatpy(params)
    assert fatpy["elastic_modulus"] is None
    assert fatpy["poissons_ratio"] == pytest.approx(0.3)

    recovered = from_fatpy(fatpy)
    assert recovered["E"] is None
    assert recovered["nu"] == pytest.approx(0.3)


def test_notation_markdown_table_matches_mapping():
    """Verify that docs/notation.md and fatpy_keys.json do not drift."""
    mapping = load_fatpy_mapping()
    text = DOCS_NOTATION.read_text(encoding="utf-8")

    # Find the table in docs/notation.md
    assert "## FatPy nomenclature" in text or "FatPy nomenclature" in text
    table_lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith("|") and not line.strip().startswith("|---")
    ]
    # Filter lines belonging to the FatPy section
    in_fatpy = False
    fatpy_rows = []
    for line in text.splitlines():
        if "FatPy key (v01)" in line:
            in_fatpy = True
            continue
        if in_fatpy:
            if line.startswith("|") and not line.startswith("|---"):
                parts = [c.strip() for c in line.split("|")[1:-1]]
                if len(parts) >= 3:
                    fatpy_rows.append(parts)
            elif not line.startswith("|") and line.strip():
                # End of table
                break

    assert len(fatpy_rows) >= 15, "Failed to parse FatPy table from docs/notation.md"

    for parts in fatpy_rows:
        raw_names = parts[0]
        # Names may be backtick-quoted: e.g. `E` or `sigma_f`, `b_fat`
        names = [n.strip(" `") for n in raw_names.split(",")]
        for name in names:
            if name.endswith("…"):
                continue
            assert (
                name in mapping
            ), f"Multifat parameter '{name}' in docs/notation.md not in fatpy_keys.json"


def test_no_private_cross_package_imports_in_cli_tools_examples():
    """Enforce Issue #81 rule: no 'from ... import _...' across packages."""
    targets = [REPO_ROOT / "src/multifat/cli", REPO_ROOT / "tools", REPO_ROOT / "examples"]
    violations = []

    for base in targets:
        for py_file in base.glob("**/*.py"):
            tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom):
                    mod = node.module or ""
                    # Check imports from multifat across packages
                    for alias in node.names:
                        if alias.name.startswith("_") and not alias.name.startswith("__"):
                            violations.append(
                                f"{py_file.relative_to(REPO_ROOT)}: imports private symbol {alias.name} from {mod}"
                            )

    assert not violations, "Found private cross-package imports:\n" + "\n".join(violations)


def test_no_model_name_branches_in_run_cli():
    """Enforce Issue #81 rule: no model-name literal in cli/run.py outside Findley amplitude."""
    run_file = REPO_ROOT / "src/multifat/cli/run.py"
    text = run_file.read_text(encoding="utf-8")

    assert "LOAD_ERROR_SKIP" not in text
    assert 'in {"BP", "CAIM"}' not in text
    assert 'in ("BP", "CAIM")' not in text
