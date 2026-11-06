"""CSV writer helper (standard library only) and the CSV flavour of the synthetic tube generator."""

from __future__ import annotations

import ast
import csv
from pathlib import Path

import pytest

from multifat.fe import csvio, csvwriter
from multifat.fe.csvio import read_header, resolve_columns
from multifat.fe.csvwriter import REQUIRED_HEADER_KEYS, FeCsvWriter
from multifat.fe.synthetic import (
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    main,
    tube_cartesian,
    tube_points,
    write_tube_csv,
)

GEOM, MAT = TubeGeometry(a=6.0, b=7.5), TubeMaterial(E=210000.0, nu=0.3)
LOADS = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 12000.0, 8e4)]
HEADER = {
    "format": "multifat-fe-csv 1.0",
    "solver": "unit test",
    "solver_version": "1",
    "stress_unit": "MPa",
    "length_unit": "mm",
    "coordinate_system": "global_cartesian",
    "averaging": "none",
    "analysis_type": "linear_elastic",
    "step_role": "time_sample",
}
SIX = ("xx", "yy", "zz", "xy", "yz", "xz")


def _read(path):
    """Rows of a CSV file as dicts (standard library only), after the comment header."""
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [ln for ln in handle if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def test_writer_uses_only_the_standard_library():
    tree = ast.parse(Path(csvwriter.__file__).read_text(encoding="utf-8"))
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0, "no relative imports: the helper must work standalone"
            roots.add((node.module or "").split(".")[0])
    assert roots <= {"__future__", "collections", "math", "pathlib", "typing"}, roots


def test_required_header_keys_match_the_importer():
    assert tuple(REQUIRED_HEADER_KEYS) == tuple(csvio.REQUIRED_HEADER_KEYS)


def test_writer_header_and_rows_exact_text(tmp_path):
    path = tmp_path / "a.csv"
    header = dict(HEADER, notes="two nodes")
    with FeCsvWriter(path, header, columns=["node_id", "step", "S_xx", "material"]) as w:
        w.write_row({"node_id": 7, "step": 1, "S_xx": 0.1, "material": "steel"})
        w.write_rows([{"node_id": 8, "step": 1, "S_xx": -2.5e-3, "material": "alu"}])
    expected = "".join(f"# {k}: {v}\n" for k, v in header.items())
    expected += "node_id,step,S_xx,material\n"
    expected += "7,1,0.10000000000000001,steel\n"
    expected += "8,1,-0.0025000000000000001,alu\n"
    assert path.read_bytes().decode("utf-8") == expected


def test_writer_puts_format_first_whatever_the_order(tmp_path):
    path = tmp_path / "a.csv"
    header = dict(reversed(list(HEADER.items())))
    with FeCsvWriter(path, header, columns=["node_id"]) as w:
        w.write_row({"node_id": 1})
    assert path.read_text(encoding="utf-8").splitlines()[0] == "# format: multifat-fe-csv 1.0"


def test_writer_rejects_missing_required_header_key(tmp_path):
    for key in REQUIRED_HEADER_KEYS:
        header = {k: v for k, v in HEADER.items() if k != key}
        with pytest.raises(ValueError, match=key):
            FeCsvWriter(tmp_path / "a.csv", header, columns=["node_id"])


def test_writer_refuses_bad_rows_and_values(tmp_path):
    with FeCsvWriter(tmp_path / "a.csv", HEADER, columns=["node_id", "S_xx"]) as w:
        with pytest.raises(ValueError, match="S_xx"):
            w.write_row({"node_id": 1})
        with pytest.raises(ValueError, match="extra"):
            w.write_row({"node_id": 1, "S_xx": 1.0, "extra": 2})
        with pytest.raises(ValueError, match="finite"):
            w.write_row({"node_id": 1, "S_xx": float("nan")})
        with pytest.raises(ValueError, match="comma"):
            w.write_row({"node_id": 1, "S_xx": "a,b"})
    with pytest.raises(ValueError, match="newline"):
        FeCsvWriter(tmp_path / "b.csv", dict(HEADER, notes="x\ny"), columns=["node_id"])
    with pytest.raises(ValueError, match="duplicate"):
        FeCsvWriter(tmp_path / "c.csv", HEADER, columns=["node_id", "node_id"])


def test_writer_output_parses_with_the_importer_header(tmp_path):
    path = tmp_path / "a.csv"
    with FeCsvWriter(path, HEADER, columns=["node_id"]) as w:
        w.write_row({"node_id": 1})
    header, n_lines = read_header(path)
    assert header.solver == "unit test" and n_lines == len(HEADER)


def test_write_tube_csv_flavours_use_the_solver_names(tmp_path):
    kw = dict(n_radii=2, n_theta=2, strain="split", with_instance=True, with_material=True)
    cols = {}
    for flavour in ("canonical", "abaqus", "ansys"):
        path = write_tube_csv(tmp_path / f"{flavour}.csv", GEOM, MAT, LOADS, flavour=flavour, **kw)
        with open(path, encoding="utf-8") as handle:
            cols[flavour] = [ln for ln in handle if not ln.startswith("#")][0].strip().split(",")
        resolve_columns(cols[flavour])  # every flavour resolves
        read_header(path)
    assert cols["canonical"] == (
        ["node_id", "instance", "material", "step", "x", "y", "z"]
        + [f"S_{c}" for c in SIX]
        + [f"Eel_{c}" for c in SIX]
        + [f"Epl_{c}" for c in SIX]
    )
    assert cols["abaqus"][-18:] == (
        ["S11", "S22", "S33", "S12", "S23", "S13"]
        + ["EE11", "EE22", "EE33", "EE12", "EE23", "EE13"]
        + ["PE11", "PE22", "PE33", "PE12", "PE23", "PE13"]
    )
    assert cols["ansys"][-18:] == (
        ["SX", "SY", "SZ", "SXY", "SYZ", "SXZ"]
        + ["EPELX", "EPELY", "EPELZ", "EPELXY", "EPELYZ", "EPELXZ"]
        + ["EPPLX", "EPPLY", "EPPLZ", "EPPLXY", "EPPLYZ", "EPPLXZ"]
    )
    total = write_tube_csv(
        tmp_path / "t.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=2, flavour="abaqus", strain="total"
    )
    with open(total, encoding="utf-8") as handle:
        line = [ln for ln in handle if not ln.startswith("#")][0].strip().split(",")
    assert line[-6:] == ["E11", "E22", "E33", "E12", "E23", "E13"]


def test_engineering_shear_doubles_only_strain_columns(tmp_path):
    kw = dict(n_radii=2, n_theta=2, flavour="canonical", strain="total")
    tensor = write_tube_csv(tmp_path / "t.csv", GEOM, MAT, LOADS, shear="tensor", **kw)
    eng = write_tube_csv(tmp_path / "e.csv", GEOM, MAT, LOADS, shear="engineering", **kw)
    assert read_header(tensor)[0].strain_convention == "tensor"
    assert read_header(eng)[0].strain_convention == "engineering"
    ids, xyz, r, theta = tube_points(GEOM, 2, 2)
    for path, factor in ((tensor, 1.0), (eng, 2.0)):
        rows = _read(path)
        assert len(rows) == 8  # 4 nodes x 2 steps, step-major
        for k, load in enumerate(LOADS):
            S, E = tube_cartesian(r, theta, load, GEOM, MAT)
            for n, row in enumerate(rows[4 * k : 4 * k + 4]):
                for c in SIX:
                    i, j = "xyz".index(c[0]), "xyz".index(c[1])
                    assert float(row[f"S_{c}"]) == pytest.approx(
                        S[n, i, j], abs=1e-9
                    )  # never doubled
                    scale = factor if i != j else 1.0
                    assert float(row[f"E_{c}"]) == pytest.approx(scale * E[n, i, j], abs=1e-15)


def test_split_strain_is_hooke_and_a_known_plastic_fraction(tmp_path):
    path = write_tube_csv(
        tmp_path / "s.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=2, analysis_type="elastic_plastic",
        strain="split", plastic_fraction=0.25,
    )  # fmt: skip
    header = read_header(path)[0]
    assert header.analysis_type == "elastic_plastic" and header.cycle == "stabilised"
    ids, xyz, r, theta = tube_points(GEOM, 2, 2)
    S, E = tube_cartesian(r, theta, LOADS[1], GEOM, MAT)
    row = _read(path)[4]  # first node, second step
    assert float(row["Eel_xy"]) == pytest.approx(E[0, 0, 1], abs=1e-15)
    assert float(row["Epl_xy"]) == pytest.approx(0.25 * E[0, 0, 1], abs=1e-15)


def test_units_scaled(tmp_path):
    path = write_tube_csv(
        tmp_path / "u.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=2, strain="none",
        stress_unit="Pa", length_unit="m",
    )  # fmt: skip
    header = read_header(path)[0]
    assert (header.stress_unit, header.length_unit) == ("Pa", "m")
    ids, xyz, r, theta = tube_points(GEOM, 2, 2)
    S, _ = tube_cartesian(r, theta, LOADS[1], GEOM, MAT)
    row = _read(path)[4]
    assert float(row["S_zz"]) == pytest.approx(S[0, 2, 2] * 1e6, rel=1e-12)
    assert float(row["x"]) == pytest.approx(xyz[0, 0] * 1e-3, rel=1e-12)


def test_node_major_order_groups_by_location(tmp_path):
    path = write_tube_csv(
        tmp_path / "n.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=2, strain="none", order="node_major",
        with_time=True,
    )  # fmt: skip
    rows = _read(path)
    assert [(int(r["node_id"]), int(r["step"])) for r in rows] == [
        (n, s) for n in (1, 2, 3, 4) for s in (1, 2)
    ]
    assert [float(r["time"]) for r in rows[:2]] == [0.0, 1.0]
    assert read_header(path)[0].time_unit == "s"


def test_header_counts_and_custom_header_override(tmp_path):
    path = write_tube_csv(
        tmp_path / "h.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=3, strain="none",
        header={"notes": "custom", "solver": "mine"},
    )  # fmt: skip
    header = read_header(path)[0]
    assert (header.n_nodes, header.n_steps) == (6, 2)
    assert header.notes == "custom" and header.solver == "mine"


def test_example_flag_writes_the_two_step_example(tmp_path):
    assert main(["--csv", str(tmp_path)]) == 0
    path = tmp_path / "tube_two_steps.csv"
    header, _ = read_header(path)
    assert (header.n_nodes, header.n_steps) == (4, 2)
    assert len(_read(path)) == 8
    with pytest.raises(SystemExit):
        main([])  # one of --out and --csv is required
