"""The history package stays apart from existing code; its example and documentation are checked."""

from __future__ import annotations

import ast
import json
import pathlib

import numpy as np

from multifat.history.eligibility import method_matrix

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "multifat"
HIST = SRC / "history"


def _imports(path: pathlib.Path):
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


def _is(module: str, package: str) -> bool:
    return module == package or module.startswith(package + ".")


def test_no_existing_module_imports_history():
    offenders = [
        str(p.relative_to(SRC))
        for p in SRC.rglob("*.py")
        if HIST not in p.parents and any(_is(m, "multifat.history") for m in _imports(p))
    ]
    assert offenders == []


def test_history_imports_no_cli_report_or_workbook_code():
    forbidden = ("cli", "report", "imp")
    offenders = {
        str(p.relative_to(SRC)): m
        for p in HIST.rglob("*.py")
        for m in _imports(p)
        if any(_is(m, f"multifat.{f}") for f in forbidden)
    }
    assert offenders == {}


def test_documented_example_runs(tmp_path, capsys):
    from multifat.history.example import main
    from multifat.history.signal import TimeHistory

    assert main(["--out", str(tmp_path)]) == 0
    assert "blocks to failure" in capsys.readouterr().out
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert 0.0 < result["blocks_to_failure"] < float("inf") and result["method"] == "FIN"
    history = TimeHistory.from_record(tmp_path / "history")
    assert history.sha256() == result["history_sha256"]
    assert history.provenance["source"] == "function" and np.isfinite(history.S_hist).all()


def test_doc_lists_every_matrix_method_and_the_example_command():
    text = (ROOT / "docs" / "load-histories.md").read_text(encoding="utf-8")
    missing = [r["method"] for r in method_matrix() if f"`{r['method']}`" not in text]
    assert missing == []
    assert "python -m multifat.history.example" in text
    for source in ("IneligibleHistoryError", "allow_nonproportional_approximation", "SNCurve"):
        assert source in text
