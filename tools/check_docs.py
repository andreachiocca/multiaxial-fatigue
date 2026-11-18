"""Check local Markdown destinations and fenced tool commands without network access.

External links, code spans and frozen literature/provenance documents are excluded.
This checks the maintained documentation's file and fragment links, and parses fenced
``tools/*.py`` commands against those tools' argument parsers; it does not certify the
availability or scientific accuracy of external sources.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
from pathlib import Path
import re
import shlex
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
FROZEN = {"docs/research", "docs/catalogue"}
SHELL_FENCES = {"", "bash", "sh", "shell", "console"}
RUNNERS = {"uv", "run", "python", "python3"}


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__.splitlines()[0])


def prose(text):
    text = re.sub(r"^```.*?^```[^\n]*", lambda m: "\n" * m[0].count("\n"), text, flags=re.M | re.S)
    return re.sub(r"`[^`\n]+`", "", text)


def anchors(text):
    """Explicit anchors and GitHub-style heading slugs, including duplicates."""
    found = set(re.findall(r"<a\s+(?:id|name)=[\"\']([^\"\']+)", text))
    counts = {}
    for title in re.findall(r"^#{1,6}\s+(.+?)\s*#*\s*$", text, flags=re.M):
        title = re.sub(r"<[^>]*>", "", title).lower()
        slug = re.sub(r"[^\w\- ]", "", title).replace(" ", "-")
        n = counts.get(slug, 0)
        counts[slug] = n + 1
        found.add(f"{slug}-{n}" if n else slug)
    return found


def _is_runner_prefix(tokens: list[str]) -> bool:
    return all(tok in RUNNERS or ("=" in tok and not tok.startswith("-")) for tok in tokens)


def fenced_tool_commands(text: str) -> list[tuple[int, str, list[str]]]:
    """Return ``(lineno, script_path, argv)`` for ``tools/*.py`` calls in fenced shell blocks."""
    found: list[tuple[int, str, list[str]]] = []
    in_fence = False
    fence_lang = ""
    pending = ""
    pending_line = 0

    def process(cmd_line: str, cmd_lineno: int) -> None:
        if "tools/" not in cmd_line:
            return
        try:
            tokens = shlex.split(cmd_line, comments=True)
        except ValueError:
            return
        groups: list[list[str]] = [[]]
        for tok in tokens:
            if tok in {"&&", ";"}:
                groups.append([])
            else:
                groups[-1].append(tok)
        for group in groups:
            for idx, tok in enumerate(group):
                if tok.startswith("tools/") and tok.endswith(".py"):
                    if _is_runner_prefix(group[:idx]):
                        found.append((cmd_lineno, tok, group[idx + 1 :]))
                    break

    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        stripped = raw_line.strip()
        if stripped.startswith("```"):
            if not in_fence:
                in_fence = True
                fence_lang = stripped[3:].strip().lower()
            else:
                if pending:
                    process(pending.strip(), pending_line)
                    pending = ""
                    pending_line = 0
                in_fence = False
                fence_lang = ""
            continue
        if not in_fence or fence_lang not in SHELL_FENCES:
            continue
        if fence_lang == "console" and stripped.startswith("$ "):
            stripped = stripped[2:].lstrip()
        if not stripped or stripped.startswith("#"):
            if pending:
                process(pending.strip(), pending_line)
                pending = ""
                pending_line = 0
            continue
        if stripped.endswith("\\"):
            if not pending:
                pending_line = lineno
            pending = (pending + " " + stripped[:-1]).strip()
            continue
        if pending:
            process((pending + " " + stripped).strip(), pending_line)
            pending = ""
            pending_line = 0
        else:
            process(stripped, lineno)
    if pending:
        process(pending.strip(), pending_line)
    return found


def _load_parser(
    script: Path, cache: dict[Path, argparse.ArgumentParser]
) -> argparse.ArgumentParser:
    resolved = script.resolve()
    if resolved in cache:
        return cache[resolved]
    if resolved == Path(__file__).resolve():
        parser = build_parser()
        cache[resolved] = parser
        return parser
    spec = importlib.util.spec_from_file_location(f"_check_docs_{script.stem}", resolved)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "build_parser"):
        raise RuntimeError(f"{script.name} does not expose build_parser()")
    parser = module.build_parser()
    cache[resolved] = parser
    return parser


def check_tool_commands(
    name: str,
    text: str,
    root: Path = ROOT,
    parser_cache: dict[Path, argparse.ArgumentParser] | None = None,
) -> tuple[int, list[str]]:
    """Validate fenced ``tools/*.py`` commands in one Markdown document."""
    cache = parser_cache if parser_cache is not None else {}
    errors: list[str] = []
    checked = 0
    for lineno, script_rel, argv in fenced_tool_commands(text):
        checked += 1
        script = root / script_rel
        if not script.is_file():
            errors.append(f"{name}:{lineno}: missing {script_rel}")
            continue
        parser = _load_parser(script, cache)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                parser.parse_args(argv)
        except SystemExit as exc:
            if exc.code != 0:
                lines = [ln for ln in buf.getvalue().splitlines() if ln.strip()]
                detail = lines[-1] if lines else f"exit {exc.code}"
                errors.append(f"{name}:{lineno}: invalid {script_rel} command: {detail}")
    return checked, errors


def main():
    build_parser().parse_args()
    tracked = (
        subprocess.check_output(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=ROOT
        )
        .decode()
        .split("\0")
    )
    # Include untracked work in progress, respecting the repository ignore rules.
    paths = sorted(set(tracked))
    errors, checked, checked_cmds = [], 0, 0
    cache = {}
    parser_cache: dict[Path, argparse.ArgumentParser] = {}
    for name in paths:
        if not name.endswith(".md") or any(name.startswith(p + "/") for p in FROZEN):
            continue
        source = ROOT / name
        if not source.is_file():
            continue
        raw_text = source.read_text(encoding="utf-8")
        cmd_count, cmd_errors = check_tool_commands(name, raw_text, ROOT, parser_cache)
        checked_cmds += cmd_count
        errors.extend(cmd_errors)
        text = prose(raw_text)
        for match in re.finditer(r"\[[^\]\n]*\]\(([^\s]+?)(?:\s+\"[^\"]*\")?\)", text):
            dest = match.group(1)
            url = urlsplit(dest)
            if url.scheme or url.netloc:
                continue
            target = (source.parent / unquote(url.path)).resolve() if url.path else source
            line = text.count("\n", 0, match.start()) + 1
            checked += 1
            if not target.exists():
                errors.append(f"{name}:{line}: missing {dest}")
            elif url.fragment and target.suffix == ".md":
                if target not in cache:
                    cache[target] = anchors(target.read_text(encoding="utf-8"))
                if unquote(url.fragment) not in cache[target]:
                    errors.append(f"{name}:{line}: missing fragment {dest}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(
        f"{checked} local documentation links and {checked_cmds} tool commands OK "
        "(external links and frozen provenance excluded)"
    )


if __name__ == "__main__":
    main()
