"""Add the column "Observation type" to the Dataset sheet of a workbook (decision d113, issue #18).

The column states what each row is: ``failure``, ``run-out`` (a test stopped without failure: a
censored life) or ``fatigue limit`` (an established limit state, from a staircase or derived from
an S-N curve). Run-out = "No" is always a failure. For Run out = "Yes" the curators' rule of
d113 applies to the workbooks collected before the column existed:

  - a workbook with at least one failure holds S-N curves: its run-outs are interrupted tests,
    ``run-out``;
  - a workbook without any failure holds fatigue limits only: its rows are ``fatigue limit``.

The edit works on the XML of the sheet, so formulas, their cached values, styles, comments and
every other part of the file stay as they are; only the Dataset sheet and, when its Excel table
ends at the last column, that table's definition are rewritten. The new column follows the last
header; an Excel table that ends there is extended to it.

Usage: uv run python tools/add_observation_type.py WORKBOOK [WORKBOOK ...] [--check]
"""

from __future__ import annotations

import argparse
import posixpath
import re
import sys
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

HEADER = "Observation type (failure/run-out/fatigue limit)"
TYPES = ("failure", "run-out", "fatigue limit")

_CELL = re.compile(r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', re.S)
_ROW = re.compile(r'<row r="(\d+)"([^>]*)>(.*?)</row>', re.S)


def col_index(letters: str) -> int:
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n


def col_letters(n: int) -> str:
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    xml = z.read("xl/sharedStrings.xml").decode("utf-8")
    return [
        "".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
        for si in re.findall(r"<si>(.*?)</si>", xml, re.S)
    ]


def _cell_text(attrs: str, body: str | None, shared: list[str]) -> str:
    if not body:
        return ""
    if 't="s"' in attrs:
        v = re.search(r"<v>(.*?)</v>", body, re.S)
        return shared[int(v.group(1))] if v else ""
    if 't="inlineStr"' in attrs:
        return "".join(re.findall(r"<t[^>]*>(.*?)</t>", body, re.S))
    v = re.search(r"<v>(.*?)</v>", body, re.S)
    return v.group(1) if v else ""


def _sheet_path(z: zipfile.ZipFile, name: str) -> str:
    wb = z.read("xl/workbook.xml").decode("utf-8")
    rid = re.search(rf'<sheet [^>]*name="{name}"[^>]*r:id="([^"]+)"', wb).group(1)
    rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8")
    target = next(
        re.search(r'Target="([^"]+)"', rel).group(1)
        for rel in re.findall(r"<Relationship [^>]*/>", rels)
        if f'Id="{rid}"' in rel
    )
    return target.lstrip("/") if target.startswith("/") else posixpath.join("xl", target)


def _table_paths(z: zipfile.ZipFile, sheet: str) -> list[str]:
    rels = posixpath.join(posixpath.dirname(sheet), "_rels", posixpath.basename(sheet) + ".rels")
    if rels not in z.namelist():
        return []
    xml = z.read(rels).decode("utf-8")
    return [
        posixpath.normpath(posixpath.join(posixpath.dirname(sheet), t))
        for t in re.findall(r'Type="[^"]*/table" Target="([^"]+)"', xml)
        + re.findall(r'Target="([^"]+)" [^>]*Type="[^"]*/table"', xml)
    ]


def observation_types(rows: dict[int, str]) -> dict[int, str]:
    """Type of every row from its Run out value (d113): see the module docstring."""
    yes = {r for r, v in rows.items() if v.strip().lower().startswith("yes")}
    no = {r for r, v in rows.items() if v.strip().lower().startswith("no")}
    unclear = set(rows) - yes - no
    if unclear:
        raise ValueError(f"rows {sorted(unclear)}: Run out is neither Yes nor No")
    runout_type = "run-out" if no else "fatigue limit"
    return {**{r: "failure" for r in no}, **{r: runout_type for r in yes}}


def add_column(path: Path, output: Path | None = None) -> dict:
    """Write `path` (or `output`) with the Observation type column; return what was done."""
    output = output or path
    with zipfile.ZipFile(path) as z:
        members = [(info, z.read(info.filename)) for info in z.infolist()]
        shared = _shared_strings(z)
        sheet = _sheet_path(z, "Dataset")
        tables = _table_paths(z, sheet)
    data = dict((info.filename, blob) for info, blob in members)
    xml = data[sheet].decode("utf-8")

    rows = {int(m.group(1)): m for m in _ROW.finditer(xml)}
    header = {
        col_index(c.group(1)): _cell_text(c.group(3), c.group(4), shared)
        for c in _CELL.finditer(rows[1].group(3))
    }
    if any(h.strip().lower().startswith("observation type") for h in header.values()):
        raise ValueError(f"{path.name}: the Observation type column already exists")
    last = max(i for i, h in header.items() if h.strip())
    new = last + 1
    new_l = col_letters(new)
    run_col = next(
        i for i, h in header.items() if "run" in h.lower() and "out" in h.lower() and h.strip()
    )

    runout = {}
    for r, m in rows.items():
        if r == 1:
            continue
        cells = {col_index(c.group(1)): c for c in _CELL.finditer(m.group(3))}
        if new in cells and (cells[new].group(4) or "").strip():
            raise ValueError(f"{path.name}: cell {new_l}{r} is not empty")
        if run_col in cells:
            text = _cell_text(cells[run_col].group(3), cells[run_col].group(4), shared).strip()
            if text:
                runout[r] = text
    types = observation_types(runout)

    def new_row(m: re.Match) -> str:
        r = int(m.group(1))
        value = HEADER if r == 1 else types.get(r)
        if value is None:  # no observation on this row: left as it is
            return m.group(0)
        attrs, body = m.group(2), m.group(3)
        cells = [c for c in _CELL.finditer(body) if col_index(c.group(1)) != new]
        ref_cell = next((c for c in cells if col_index(c.group(1)) == last), None)
        style = re.search(r'\ss="(\d+)"', ref_cell.group(3)) if ref_cell else None
        s_attr = f' s="{style.group(1)}"' if style else ""
        cell = f'<c r="{new_l}{r}"{s_attr} t="inlineStr"><is><t>{escape(value)}</t></is></c>'
        before = [c.group(0) for c in cells if col_index(c.group(1)) < new]
        after = [c.group(0) for c in cells if col_index(c.group(1)) > new]
        spans = re.search(r'spans="(\d+):(\d+)"', attrs)
        if spans and int(spans.group(2)) < new:
            attrs = attrs.replace(spans.group(0), f'spans="{spans.group(1)}:{new}"')
        return f'<row r="{r}"{attrs}>{"".join(before)}{cell}{"".join(after)}</row>'

    xml = _ROW.sub(new_row, xml)
    dim = re.search(r'<dimension ref="([A-Z]+)(\d+):([A-Z]+)(\d+)"', xml)
    if dim and col_index(dim.group(3)) < new:
        xml = xml.replace(
            dim.group(0), f'<dimension ref="{dim.group(1)}{dim.group(2)}:{new_l}{dim.group(4)}"'
        )
    data[sheet] = xml.encode("utf-8")

    extended = []
    for table in tables:
        txml = data[table].decode("utf-8")
        ref = re.search(r'\bref="([A-Z]+)(\d+):([A-Z]+)(\d+)"', txml)
        if col_index(ref.group(3)) != last:
            continue
        new_ref = f"{ref.group(1)}{ref.group(2)}:{new_l}{ref.group(4)}"
        txml = txml.replace(f'ref="{ref.group(1)}{ref.group(2)}:{ref.group(3)}{ref.group(4)}"',
                            f'ref="{new_ref}"')  # fmt: skip
        count = re.search(r'<tableColumns count="(\d+)"', txml)
        n = int(count.group(1))
        ids = [int(i) for i in re.findall(r'<tableColumn id="(\d+)"', txml)]
        txml = txml.replace(count.group(0), f'<tableColumns count="{n + 1}"')
        txml = txml.replace(
            "</tableColumns>",
            f'<tableColumn id="{max(ids) + 1}" name="{escape(HEADER)}"/></tableColumns>',
        )
        data[table] = txml.encode("utf-8")
        extended.append(table)

    tmp = output.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(tmp, "w") as out:
        for info, _ in members:
            out.writestr(info, data[info.filename])
    tmp.replace(output)
    counts = {t: sum(v == t for v in types.values()) for t in TYPES}
    return dict(column=new_l, tables=extended, **counts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workbooks", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.workbooks:
        try:
            print(path.name, add_column(path))
        except ValueError as exc:
            print(f"{path.name}: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
