"""Writer for the FE CSV exchange format, standard library only so solver scripts can copy it.

``FeCsvWriter`` writes the ``# key: value`` header, the column line and the rows; floats are
written with 17 significant digits so a value round-trips exactly. It checks that the required
header keys are present and that rows match the declared columns, and nothing else: units and
conventions are the caller's declaration, validated when the file is imported
(``multifat.fe.csvio``). The format is specified in ``docs/fe-csv-format.md``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

REQUIRED_HEADER_KEYS = (
    "format",
    "solver",
    "solver_version",
    "stress_unit",
    "length_unit",
    "coordinate_system",
    "averaging",
    "analysis_type",
    "step_role",
)


def _text(value: object, what: str) -> str:
    text = str(value)
    if "\n" in text or "\r" in text:
        raise ValueError(f"{what} contains a newline")
    return text


class FeCsvWriter:
    """Write an FE CSV file row by row.

    Parameters
    ----------
    path : str or Path
        Output file (UTF-8, ``\\n`` line ends).
    header : mapping
        Header keys and values; the required keys must be present. ``format`` is written first.
    columns : sequence of str
        Column names, in the order of the file; each row must give exactly these.
    """

    def __init__(self, path: str | Path, header: Mapping[str, object], *, columns: Sequence[str]):
        missing = [k for k in REQUIRED_HEADER_KEYS if k not in header]
        if missing:
            raise ValueError(f"missing required header key(s): {', '.join(missing)}")
        names = [_text(c, "column name") for c in columns]
        if len(set(names)) != len(names):
            raise ValueError(f"duplicate column names in {names}")
        for name in names:
            if "," in name:
                raise ValueError(f"column name {name!r} contains a comma")
        self._columns = tuple(names)
        keys = ["format"] + [k for k in header if k != "format"]
        lines = []
        for key in keys:
            value = _text(header[key], f"header value of {key!r}")
            lines.append(f"# {key}: {value}\n")
        self._handle = open(path, "w", encoding="utf-8", newline="")
        self._handle.write("".join(lines))
        self._handle.write(",".join(self._columns) + "\n")

    def write_row(self, row: Mapping[str, object]) -> None:
        """Append one row; ``row`` gives exactly the declared columns."""
        absent = [c for c in self._columns if c not in row]
        if absent:
            raise ValueError(f"row lacks column(s): {', '.join(absent)}")
        extra = [c for c in row if c not in self._columns]
        if extra:
            raise ValueError(f"row has extra column(s): {', '.join(map(str, extra))}")
        cells = []
        for column in self._columns:
            value = row[column]
            if isinstance(value, str):
                if "," in value:
                    raise ValueError(f"column {column!r}: a comma in {value!r} (no quoting)")
                cells.append(_text(value, f"column {column!r}"))
            elif hasattr(value, "__index__"):  # int, bool and integer types of numpy
                cells.append(str(int(value)))
            else:
                number = float(value)
                if not math.isfinite(number):
                    raise ValueError(f"column {column!r}: value must be finite, got {number}")
                cells.append(format(number, ".17g"))
        self._handle.write(",".join(cells) + "\n")

    def write_rows(self, rows: Iterable[Mapping[str, object]]) -> None:
        """Append several rows."""
        for row in rows:
            self.write_row(row)

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "FeCsvWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
