"""FE histories stay apart from experimental data: import graph and entry-point rejection."""

from __future__ import annotations

import ast
import pathlib

import numpy as np
import pytest

from multifat.fe.model import HistoryBlock

SRC = pathlib.Path(__file__).resolve().parents[1] / "src" / "multifat"
FE = SRC / "fe"
HISTORY = SRC / "history"  # the one intended consumer of fe (docs/load-histories.md)


def _imports(path: pathlib.Path):
    """Absolute module names imported by a file; relative imports resolve inside its package."""
    pkg = ["multifat"] + list(path.relative_to(SRC).parent.parts)
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            yield from (a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = pkg[: len(pkg) - (node.level - 1)]
                yield ".".join(base + ([node.module] if node.module else []))
            elif node.module:
                yield node.module


def test_no_existing_module_imports_fe():
    """Only ``fe`` itself and ``history``, built on it, import ``multifat.fe``."""
    offenders = [
        str(p.relative_to(SRC))
        for p in SRC.rglob("*.py")
        if FE not in p.parents
        and HISTORY not in p.parents
        and any(m == "multifat.fe" or m.startswith("multifat.fe.") for m in _imports(p))
    ]
    assert offenders == []


def test_fe_imports_no_workbook_calibration_model_or_report_modules():
    forbidden = ("imp", "eval", "report", "cli", "models", "planes")
    lazy_in_bridge = ("bridge.py", "multifat.models.registry")  # stress-only check, imported on use
    offenders = {
        str(p.relative_to(SRC)): m
        for p in FE.rglob("*.py")
        for m in _imports(p)
        if any(m == f"multifat.{f}" or m.startswith(f"multifat.{f}.") for f in forbidden)
        and (p.name, m) != lazy_in_bridge
    }
    assert offenders == {}


@pytest.fixture
def block():
    S = np.zeros((1, 2, 3, 3))
    return HistoryBlock(np.array([1]), np.zeros((1, 3)), S, None)


def test_workbook_and_calibration_entry_points_reject_fe_objects(block, tmp_path):
    from multifat.cli.compare import collect_all_results
    from multifat.imp import experimental_data as exp

    h5 = tmp_path / "history.h5"
    h5.write_bytes(b"not a workbook")
    with pytest.raises(Exception, match="xlsx|extension|format") as info:
        exp.load_material_workbook(h5)
    assert not isinstance(info.value, AttributeError)
    with pytest.raises(FileNotFoundError, match="xlsx"):
        exp.import_experimental_data_grouped_by_material(tmp_path)
    with pytest.raises(TypeError):
        exp.load_material_workbook(block)
    with pytest.raises(TypeError):
        exp.extract_tensor(block, "m", 0)
    with pytest.raises(TypeError):
        exp.get_material_properties(block, "m")
    with pytest.raises(FileNotFoundError, match="CSV"):
        collect_all_results(tmp_path)


def test_fe_csv_modules_import_no_workbook_calibration_or_history_code():
    forbidden = ("imp", "eval", "models", "history", "cli", "report", "planes")
    offenders = {
        name: m
        for name in ("csvio.py", "csvwriter.py")
        for m in _imports(FE / name)
        if any(m == f"multifat.{f}" or m.startswith(f"multifat.{f}.") for f in forbidden)
    }
    assert offenders == {}


def test_header_parser_and_writer_work_without_h5py_and_the_writer_without_numpy(tmp_path):
    """The writer is copyable into a solver script; reading a header needs only the parser."""
    import subprocess
    import sys

    from multifat.fe.synthetic import TubeGeometry, TubeLoad, TubeMaterial, write_tube_csv

    csv = write_tube_csv(
        tmp_path / "t.csv", TubeGeometry(6.0, 7.5), TubeMaterial(210000.0, 0.3),
        [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 0.0, 0.0)], n_radii=2, n_theta=2, strain="none",
    )  # fmt: skip
    writer = FE / "csvwriter.py"
    code = (
        "import importlib.util, sys\n"
        "sys.modules['numpy'] = None\n"
        "sys.modules['h5py'] = None\n"
        f"spec = importlib.util.spec_from_file_location('standalone_writer', r'{writer}')\n"
        "module = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(module)\n"
        "print('writer-ok')\n"
    )
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "writer-ok" in out.stdout, out.stderr
    code = (
        "import sys\n"
        "sys.modules['h5py'] = None\n"
        "from multifat.fe.csvio import read_header\n"
        f"header, n = read_header(r'{csv}')\n"
        "print('header-ok', header.solver.split()[0], n)\n"
    )
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and "header-ok analytical" in out.stdout, out.stderr


def test_csv_doc_lists_every_header_key_alias_and_the_example_command():
    from multifat.fe.csvio import COLUMN_ALIASES, REQUIRED_HEADER_KEYS

    text = (SRC.parent.parent / "docs" / "fe-csv-format.md").read_text(encoding="utf-8")
    missing_keys = [k for k in REQUIRED_HEADER_KEYS if f"`{k}`" not in text]
    assert missing_keys == []
    for optional in (
        "strain_convention",
        "strain_unit",
        "cycle",
        "time_unit",
        "n_nodes",
        "n_steps",
    ):
        assert f"`{optional}`" in text, optional
    # every alias, in the solver's own spelling, appears in the alias table
    spelled = [
        "S11",
        "S13",
        "EE12",
        "PE23",
        "E11",
        "SX",
        "SXZ",
        "EPELXY",
        "EPPLYZ",
        "EPTOZ",
        "TEMP",
    ]
    assert [a for a in spelled if f"`{a}`" not in text] == []
    assert all(alias in COLUMN_ALIASES for alias in (a.lower() for a in spelled))
    assert "python -m multifat.fe.synthetic --csv" in text
    for token in ("import_fe_csv", "assess_nodes", "FeCsvWriter", "unit_case"):
        assert token in text, token
