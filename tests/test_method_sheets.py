"""Structure checks of the LaTeX method sheets in docs/method-sheets."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SHEETS_DIR = ROOT / "docs" / "method-sheets"

# Registered models whose sheet is still to be written. This list may only shrink:
# a new method must arrive with its sheet.
PENDING = {
    "BP", "CAIM", "COIN", "CROSSLAND", "DANGVAN", "EK", "FIN", "FS", "GAM", "GOLOS", "GSA",
    "GSE", "KUJAWSKI2021", "LI", "LIU2021", "LIU2025", "LIUMAHADEVAN", "MATAKE", "MCD", "MDC",
    "MGSE_YU", "MGSE_ZHU", "MKBM", "MMCCM", "MMMH", "MMP", "MOREL", "MSWT", "MWCM", "MWHS",
    "NINIC", "OTT", "PAPADOPOULOS", "PBP", "PCR", "PI", "RFB", "SINES", "SSF", "SWT", "SWTD",
    "TVEIT2025", "VHN", "VSE", "WANG2022", "ZHANGYAO",
}  # fmt: skip

BLOCKS = (
    "Source",
    "Idea and scope",
    "Inputs",
    "Formulation",
    "Constants and calibration",
    "Life",
    "Repository choices",
    "Worked example",
    "Verification",
)


def _registered() -> set[str]:
    from multifat.models.registry import _available_cp_methods, _available_direct_models

    return set({**_available_cp_methods(), **_available_direct_models()})


def _sheets() -> dict[Path, str]:
    return {p: p.read_text(encoding="utf-8") for p in sorted((SHEETS_DIR / "sheets").glob("*.tex"))}


def _sheet_models() -> dict[str, Path]:
    found: dict[str, Path] = {}
    for path, text in _sheets().items():
        first = text.splitlines()[0]
        match = re.fullmatch(r"% registry: (.+)", first)
        assert match, f"{path.name}: first line must be '% registry: NAME[, NAME]'"
        for name in (n.strip() for n in match.group(1).split(",")):
            assert name not in found, f"{name} has two sheets"
            found[name] = path
    return found


def test_every_registered_model_has_a_sheet_or_is_pending():
    registered = _registered()
    documented = set(_sheet_models())
    assert documented <= registered, f"sheets for unknown models: {documented - registered}"
    assert not (documented & PENDING), f"remove from PENDING: {documented & PENDING}"
    assert PENDING <= registered, f"PENDING lists unknown models: {PENDING - registered}"
    missing = registered - documented - PENDING
    assert not missing, f"registered models without a sheet: {sorted(missing)}"


def test_every_sheet_is_input_by_the_master_file():
    master = (SHEETS_DIR / "methods.tex").read_text(encoding="utf-8")
    inputs = set(re.findall(r"\\input\{([^}]+)\}", master))
    for target in inputs:
        assert (SHEETS_DIR / f"{target}.tex").is_file(), f"\\input{{{target}}} has no file"
    for path in _sheets():
        assert f"sheets/{path.stem}" in inputs, f"{path.name} is not input by methods.tex"


@pytest.mark.parametrize("path", sorted(_sheets()), ids=lambda p: p.stem)
def test_sheet_has_every_block_in_order(path):
    text = path.read_text(encoding="utf-8")
    positions = [text.find(f"\\block{{{name}}}") for name in BLOCKS]
    assert all(
        p >= 0 for p in positions
    ), f"missing blocks: {[b for b, p in zip(BLOCKS, positions) if p < 0]}"
    assert positions == sorted(positions), "blocks out of order"
    assert "\\methodhead{" in text and "\\label{sheet:" in text


def test_citations_resolve_and_braces_balance():
    bib = (SHEETS_DIR / "references.bib").read_text(encoding="utf-8")
    keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    files = [SHEETS_DIR / "methods.tex", SHEETS_DIR / "notation.tex", *_sheets()]
    for path in files:
        text = re.sub(r"(?<!\\)%.*", "", path.read_text(encoding="utf-8"))
        for group in re.findall(r"\\cite\{([^}]+)\}", text):
            for key in group.split(","):
                assert key.strip() in keys, f"{path.name}: \\cite{{{key}}} not in references.bib"
        stripped = text.replace("\\{", "").replace("\\}", "")
        assert stripped.count("{") == stripped.count("}"), f"{path.name}: unbalanced braces"
        for env in set(re.findall(r"\\begin\{(\w+\*?)\}", text)):
            n_begin = text.count(f"\\begin{{{env}}}")
            assert n_begin == text.count(f"\\end{{{env}}}"), f"{path.name}: unbalanced {env}"
