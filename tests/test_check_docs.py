"""Documentation checker: validates fenced tools/*.py commands against tool parsers."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("check_docs", ROOT / "tools/check_docs.py")
check_docs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_docs)

TOOLS = sorted(p.name for p in (ROOT / "tools").glob("*.py"))


@pytest.mark.parametrize("tool", TOOLS)
def test_every_tool_supports_help(tool: str) -> None:
    out = subprocess.run(
        [sys.executable, str(ROOT / "tools" / tool), "--help"],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert out.returncode == 0, out.stderr


def test_fenced_tool_commands_accepts_valid_invocations() -> None:
    md = "\n".join(
        [
            "# Example",
            "",
            "```bash",
            "uv run python tools/generate_method_docs.py --check  # comment",
            "uv run python tools/regression_snapshot.py \\",
            "    run <dir> --jobs 4",
            "uv run python tools/regression_snapshot.py compare /tmp/a /tmp/b",
            "python tools/check_docs.py && python tools/check_scientific_reviews.py",
            "git diff tools/regression_snapshot.py",
            "```",
            "",
            "```python",
            "uv run python tools/regression_snapshot.py",
            "```",
        ]
    )
    count, errors = check_docs.check_tool_commands("example.md", md, ROOT)
    assert count == 5
    assert not errors


@pytest.mark.parametrize(
    "line",
    [
        "uv run python tools/regression_snapshot.py",
        "uv run python tools/generate_method_docs.py --unknown-flag",
        "uv run python tools/check_docs.py --unexpected",
        "uv run python tools/missing_tool.py",
    ],
)
def test_fenced_tool_commands_rejects_invalid_invocations(line: str) -> None:
    md = f"```bash\n{line}\n```\n"
    count, errors = check_docs.check_tool_commands("bad.md", md, ROOT)
    assert count == 1
    assert len(errors) == 1
