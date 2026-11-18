"""CSV importer: round trips against the tube oracle, conventions, refusals and streaming."""

from __future__ import annotations

import json
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

from multifat.fe.csvio import import_fe_csv
from multifat.fe.csvwriter import FeCsvWriter
from multifat.fe.model import sha256_of
from multifat.fe.select import select_by_ids
from multifat.fe.store import HistoryFile
from multifat.fe.synthetic import (
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    main,
    tube_cartesian,
    tube_points,
    write_tube_csv,
)

pytest.importorskip("h5py")

GEOM, MAT = TubeGeometry(a=6.0, b=7.5), TubeMaterial(E=210000.0, nu=0.3)
LOADS = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 12000.0, 0.0), TubeLoad(50.0, 12000.0, 8e4)]
N_RADII, N_THETA = 2, 3  # three angles: all six Cartesian components are nonzero and distinct

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
S_COLS = [f"S_{c}" for c in SIX]
EEL_COLS = [f"Eel_{c}" for c in SIX]
EPL_COLS = [f"Epl_{c}" for c in SIX]
E_COLS = [f"E_{c}" for c in SIX]


def _oracle(strain, fraction=0.0):
    ids, xyz, r, theta = tube_points(GEOM, N_RADII, N_THETA)
    S = np.empty((len(ids), len(LOADS), 3, 3))
    E = np.empty_like(S)
    for k, load in enumerate(LOADS):
        S[:, k], E[:, k] = tube_cartesian(r, theta, load, GEOM, MAT)
    return ids, xyz, S, (1.0 + fraction) * E, E


def _close(got, want, scale=1.0):
    np.testing.assert_allclose(got, want, rtol=1e-11, atol=1e-13 * scale)


ROUNDTRIP = [
    dict(flavour="canonical", shear="tensor"),
    dict(flavour="abaqus", shear="engineering", with_instance=True, with_material=True,
         with_time=True),
    dict(flavour="ansys", shear="engineering", stress_unit="Pa", length_unit="m",
         order="node_major"),
    dict(flavour="abaqus", shear="tensor", stress_unit="psi", length_unit="in",
         order="node_major", with_instance=True, with_time=True),
    dict(flavour="ansys", shear="tensor", stress_unit="ksi", length_unit="ft", with_material=True),
    dict(flavour="canonical", shear="engineering", stress_unit="N/mm^2", length_unit="cm",
         order="node_major", with_time=True),
]  # fmt: skip


@pytest.mark.parametrize("strain", ["split", "total", "none"])
@pytest.mark.parametrize("options", ROUNDTRIP, ids=lambda o: "-".join(map(str, o.values())))
def test_roundtrip_matches_oracle(tmp_path, options, strain):
    fraction = 0.25 if strain == "split" else 0.0
    kw = dict(options)
    if strain == "none":
        kw.pop("shear")
    path = write_tube_csv(
        tmp_path / "t.csv", GEOM, MAT, LOADS, n_radii=N_RADII, n_theta=N_THETA, strain=strain,
        analysis_type="elastic_plastic" if strain == "split" else "linear_elastic",
        plastic_fraction=fraction if strain == "split" else None, **kw,
    )  # fmt: skip
    meta = import_fe_csv(path, tmp_path / "h.h5")
    ids, xyz, S, E_total, E_elastic = _oracle(strain, fraction)
    with HistoryFile(tmp_path / "h.h5") as f:
        assert f.meta == meta
        block = f.read(slice(None))
        assert list(f.load_cases.steps) == [1, 2, 3]
        assert (f.load_cases.time is not None) == bool(options.get("with_time"))
        if f.load_cases.time is not None:
            assert np.allclose(f.load_cases.time, [0.0, 0.5, 1.0])
        assert (f.instance is not None) == bool(options.get("with_instance"))
        assert (f.material is not None) == bool(options.get("with_material"))
    assert np.array_equal(block.location_id, ids)
    _close(block.coordinates, xyz, scale=10.0)
    _close(block.S_hist, S, scale=1e3)
    if strain == "none":
        assert block.E_hist is None and meta.strain_kind is None
    else:
        _close(block.E_hist, E_total, scale=1.0)
        assert meta.strain_kind == "total_mechanical"
    if strain == "split":
        _close(block.E_elastic_hist, E_elastic)
        _close(block.E_plastic_hist, fraction * E_elastic)
    else:
        assert block.E_elastic_hist is None and block.E_plastic_hist is None


def _write(path, rows, columns, **header):
    with FeCsvWriter(path, {**HEADER, **header}, columns=columns) as writer:
        writer.write_rows(rows)
    return path


def _val(node, step, k):
    return 100.0 * node + 10.0 * step + k + 1.0  # distinct for every component


def _eps(node, step, k):
    return 1e-3 * (node + 0.1 * step + 0.01 * (k + 1))


def _rows(nodes, steps, *, strain=None, names=S_COLS, strain_names=None, extra=None):
    out = []
    for step in steps:
        for node in nodes:
            row = {"node_id": node, "step": step, "x": float(node), "y": 0.0, "z": 0.0}
            row.update({name: _val(node, step, k) for k, name in enumerate(names)})
            if strain is not None:
                row.update({name: _eps(node, step, k) for k, name in enumerate(strain_names)})
            if extra:
                row.update(extra(node, step))
            out.append(row)
    return out


BASE = ["node_id", "step", "x", "y", "z"]


@pytest.mark.parametrize("row_index", [0, 3])
@pytest.mark.parametrize("extra", [True, False])
def test_row_width_checked_before_projection(tmp_path, row_index, extra):
    path = _write(tmp_path / "bad.csv", _rows([1, 2], [1, 2]), BASE + S_COLS)
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if not line.startswith("#")) + 1
    cells = lines[start + row_index].split(",")
    cells = cells + ["999"] if extra else cells[:-1]
    lines[start + row_index] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n")
    output = tmp_path / "existing.h5"
    output.write_bytes(b"previous output")
    with pytest.raises(ValueError, match="column|field|width"):
        import_fe_csv(path, output, chunk_rows=1)
    assert output.read_bytes() == b"previous output"


def test_import_refuses_source_destination_alias(tmp_path):
    path = _write(tmp_path / "source.csv", _rows([1], [1, 2]), BASE + S_COLS)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="source|different|same"):
        import_fe_csv(path, path)
    assert path.read_bytes() == before


def test_writer_values_roundtrip_without_decimal_precision_loss(tmp_path):
    values = [
        1.2345678901234567e-7,
        3.141592653589793e-4,
        2.718281828459045e-6,
        1.4142135623730951e-5,
        1.7320508075688772e-6,
        2.23606797749979e-8,
    ]
    data = _rows([1], [1, 2], strain=True, strain_names=E_COLS)
    for row in data:
        row.update(zip(E_COLS, values))
    path = _write(
        tmp_path / "precision.csv",
        data,
        BASE + S_COLS + E_COLS,
        strain_convention="tensor",
        strain_unit="1",
    )
    block, _, _ = _load(path, tmp_path)
    expected = _tensor(lambda k: values[k])
    np.testing.assert_array_equal(block.E_hist[0, 0], expected)


def test_failed_tensor_pass_preserves_existing_destination(tmp_path):
    data = _rows([1], [1, 2])
    data[1]["S_xx"] = "nan"
    path = _write(tmp_path / "nonfinite.csv", data, BASE + S_COLS)
    output = tmp_path / "existing.h5"
    output.write_bytes(b"previous output")
    with pytest.raises(ValueError, match="finite"):
        import_fe_csv(path, output)
    assert output.read_bytes() == b"previous output"


def _tensor(f, scale=1.0):
    """Symmetric tensor with entries f(0..5) = xx, yy, zz, xy, yz, xz."""
    v = [f(k) * scale for k in range(6)]
    return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]])


def _read(path):
    with HistoryFile(path) as f:
        return f.read(slice(None)), f.meta, f.load_cases


def _load(csv, tmp_path, name="h.h5", **kw):
    """Import a CSV file and read the whole history back: ``(block, meta, load_cases)``."""
    import_fe_csv(csv, tmp_path / name, **kw)
    return _read(tmp_path / name)


@pytest.mark.parametrize(
    "names",
    [
        ["S11", "S22", "S33", "S12", "S23", "S13"],  # Abaqus: S12 -> xy, S23 -> yz, S13 -> xz
        ["SX", "SY", "SZ", "SXY", "SYZ", "SXZ"],  # Ansys
        ["Sxx", "Syy", "Szz", "Sxy", "Syz", "Sxz"],
    ],
)
def test_solver_names_map_each_component_to_its_own_tensor_entry(tmp_path, names):
    rows = _rows([1, 2], [1, 2], names=names)
    path = _write(tmp_path / "a.csv", rows, BASE + names)
    block, meta, _ = _load(path, tmp_path)
    for n, node in enumerate([1, 2]):
        for s, step in enumerate([1, 2]):
            np.testing.assert_allclose(block.S_hist[n, s], _tensor(lambda k: _val(node, step, k)))
    mapping = json.loads(meta.source_conventions["column_aliases_used"])
    assert mapping["S_xz"] == names[5] or names[5].lower() == "sxz" and "S_xz" in mapping


def test_columns_are_matched_by_name_not_position(tmp_path):
    rows = _rows([1, 2, 3], [1, 2])
    shuffled = ["S_xz", "step", "S_yy", "node_id", "z", "S_xx", "y", "S_zz", "x", "S_yz", "S_xy"]
    path = _write(tmp_path / "a.csv", rows, shuffled)
    block, *_ = _load(path, tmp_path)
    np.testing.assert_allclose(block.S_hist[2, 1], _tensor(lambda k: _val(3, 2, k)))


def test_engineering_shear_not_applied_to_stress(tmp_path):
    eps_names = ["E11", "E22", "E33", "E12", "E23", "E13"]
    s_names = ["S11", "S22", "S33", "S12", "S23", "S13"]
    out = {}
    for convention in ("tensor", "engineering"):
        path = _write(
            tmp_path / f"{convention}.csv",
            _rows([1], [1, 2], strain=True, names=s_names, strain_names=eps_names),
            BASE + s_names + eps_names,
            strain_convention=convention, strain_unit="1",
        )  # fmt: skip
        out[convention] = _load(path, tmp_path, f"{convention}.h5")[0]
    np.testing.assert_array_equal(
        out["tensor"].S_hist, out["engineering"].S_hist
    )  # stress untouched
    ref = _tensor(lambda k: _eps(1, 2, k))
    np.testing.assert_allclose(out["tensor"].E_hist[0, 1], ref)
    halved = ref.copy()
    halved[~np.eye(3, dtype=bool)] *= 0.5
    np.testing.assert_allclose(out["engineering"].E_hist[0, 1], halved)


def test_instance_repeated_labels_and_interface_node_import_as_distinct_locations(tmp_path):
    locations = [("A", 1, "steel"), ("B", 1, "steel"), ("A", 2, "steel"), ("A", 2, "alu")]
    rows = []
    for step in (1, 2):
        for n, (inst, node, mat) in enumerate(locations):
            row = {"node_id": node, "instance": inst, "material": mat, "step": step, "x": float(n),
                   "y": 0.0, "z": 0.0}  # fmt: skip
            row.update({c: _val(10 * (n + 1), step, k) for k, c in enumerate(S_COLS)})
            rows.append(row)
    path = _write(
        tmp_path / "a.csv",
        rows,
        ["node_id", "instance", "material", "step", "x", "y", "z"] + S_COLS,
    )
    meta = import_fe_csv(path, tmp_path / "h.h5")
    assert meta.material_id == "alu+steel"
    with HistoryFile(tmp_path / "h.h5") as f:
        assert f.n_locations == 4
        pos_b = select_by_ids(f, [1], instance="B")
        pos_alu = select_by_ids(f, [2], material="alu")
        assert list(select_by_ids(f, [1])) == [0, 1] and list(select_by_ids(f, [2])) == [2, 3]
        b = f.read(pos_b)
        alu = f.read(pos_alu)
    np.testing.assert_allclose(b.S_hist[0, 1], _tensor(lambda k: _val(20, 2, k)))
    np.testing.assert_allclose(alu.S_hist[0, 0], _tensor(lambda k: _val(40, 1, k)))
    assert list(b.instance) == ["B"] and list(alu.material) == ["alu"]


def test_truncated_file_refused_against_n_nodes_and_n_steps(tmp_path):
    rows = _rows([1, 2, 3], [1, 2])
    cols = BASE + S_COLS
    with pytest.raises(ValueError, match="n_nodes"):
        import_fe_csv(
            _write(tmp_path / "a.csv", rows, cols, n_nodes=4, n_steps=2), tmp_path / "h.h5"
        )
    with pytest.raises(ValueError, match="n_steps"):
        import_fe_csv(
            _write(tmp_path / "b.csv", rows, cols, n_nodes=3, n_steps=3), tmp_path / "h.h5"
        )
    assert not (tmp_path / "h.h5").exists()
    ok = _write(tmp_path / "c.csv", rows, cols, n_nodes=3, n_steps=2)
    import_fe_csv(ok, tmp_path / "h.h5")


def test_ragged_grid_refused(tmp_path):
    rows = _rows([1, 2, 3], [1, 2])
    cols = BASE + S_COLS
    last_missing = rows[:-1]  # the last node of the last step is cut off (a truncated export)
    with pytest.raises(ValueError, match="different set of locations"):
        import_fe_csv(_write(tmp_path / "a.csv", last_missing, cols), tmp_path / "h.h5")
    swapped = rows[:-1] + [dict(rows[-1], node_id=99)]
    with pytest.raises(ValueError, match="different set of locations"):
        import_fe_csv(_write(tmp_path / "b.csv", swapped, cols), tmp_path / "h.h5")
    assert not (tmp_path / "h.h5").exists()


def test_duplicate_keys_refused(tmp_path):
    rows = _rows([1, 2, 3], [1, 2])
    rows.append(dict(rows[0]))  # step 1 again at the end: a repeated step block
    with pytest.raises(ValueError, match="step"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, BASE + S_COLS), tmp_path / "h.h5")
    rows = _rows([1, 2, 3], [1, 2])
    rows[1] = dict(rows[1], node_id=1)  # node 1 twice in step 1
    with pytest.raises(ValueError, match="duplicate key"):
        import_fe_csv(_write(tmp_path / "b.csv", rows, BASE + S_COLS), tmp_path / "h.h5")


def test_unsorted_steps_in_step_major_file_refused(tmp_path):
    rows = _rows([1, 2], [2, 1])
    with pytest.raises(ValueError, match="ascend"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, BASE + S_COLS), tmp_path / "h.h5")


def test_node_major_file_is_sorted_by_step_in_memory(tmp_path):
    rows = _rows([1, 2, 3], [1, 2, 3])
    node_major = sorted(rows, key=lambda r: (r["node_id"], -r["step"]))  # steps descending per node
    path = _write(tmp_path / "a.csv", node_major, BASE + S_COLS)
    block, _, cases = _load(path, tmp_path)
    assert list(cases.steps) == [1, 2, 3]
    np.testing.assert_allclose(block.S_hist[1, 0], _tensor(lambda k: _val(2, 1, k)))
    np.testing.assert_allclose(block.S_hist[1, 2], _tensor(lambda k: _val(2, 3, k)))


def test_time_not_increasing_refused_and_never_resorted(tmp_path):
    cols = BASE + ["time"] + S_COLS
    rows = _rows([1, 2], [1, 2, 3], extra=lambda n, s: {"time": {1: 0.0, 2: 2.0, 3: 1.0}[s]})
    with pytest.raises(ValueError, match="time"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, cols, time_unit="s"), tmp_path / "h.h5")
    rows = _rows([1, 2], [1, 2], extra=lambda n, s: {"time": s + 0.1 * n})  # time varies in a step
    with pytest.raises(ValueError, match="one finite value per step"):
        import_fe_csv(_write(tmp_path / "b.csv", rows, cols, time_unit="s"), tmp_path / "h.h5")
    rows = _rows([1, 2], [1, 2], extra=lambda n, s: {"time": 1.5})  # equal times
    with pytest.raises(ValueError, match="strictly"):
        import_fe_csv(_write(tmp_path / "c.csv", rows, cols, time_unit="s"), tmp_path / "h.h5")
    with pytest.raises(ValueError, match="time_unit"):
        rows = _rows([1, 2], [1, 2], extra=lambda n, s: {"time": float(s)})
        import_fe_csv(_write(tmp_path / "d.csv", rows, cols), tmp_path / "h.h5")


def test_non_contiguous_steps_accepted_and_recorded(tmp_path):
    rows = _rows([1, 2], [3, 7, 12])
    block, _, cases = _load(_write(tmp_path / "a.csv", rows, BASE + S_COLS), tmp_path)
    assert list(cases.steps) == [3, 7, 12] and cases.labels == ("step3", "step7", "step12")
    np.testing.assert_allclose(block.S_hist[0, 2], _tensor(lambda k: _val(1, 12, k)))


def _corrupt_last_row(path, column_index, text):
    """Replace one cell of the last data row (the writer itself refuses non-finite values)."""
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    cells = lines[-1].split(",")
    cells[column_index] = text
    lines[-1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.mark.parametrize("bad", ["nan", "inf", "-inf"])
def test_non_finite_values_refused_and_partial_output_removed(tmp_path, bad):
    path = _write(tmp_path / "a.csv", _rows([1, 2, 3], [1, 2, 3]), BASE + S_COLS)
    _corrupt_last_row(path, len(BASE) + 3, bad)  # last step: the output file is already open
    with pytest.raises(ValueError, match="finite|convert"):
        import_fe_csv(path, tmp_path / "h.h5")
    assert not (tmp_path / "h.h5").exists()


def test_blank_numeric_cell_refused(tmp_path):
    path = _write(tmp_path / "a.csv", _rows([1, 2], [1, 2]), BASE + S_COLS)
    _corrupt_last_row(path, len(BASE) + 3, "")
    with pytest.raises(ValueError):
        import_fe_csv(path, tmp_path / "h.h5")
    assert not (tmp_path / "h.h5").exists()


def test_missing_strain_convention_or_strain_unit_refused_when_strain_present(tmp_path):
    rows = _rows([1, 2], [1, 2], strain=True, strain_names=E_COLS)
    cols = BASE + S_COLS + E_COLS
    with pytest.raises(ValueError, match="strain_convention"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, cols, strain_unit="1"), tmp_path / "h.h5")
    with pytest.raises(ValueError, match="strain_unit"):
        import_fe_csv(
            _write(tmp_path / "b.csv", rows, cols, strain_convention="tensor"), tmp_path / "h.h5"
        )


@pytest.mark.parametrize("unit, factor", [("percent", 1e-2), ("microstrain", 1e-6), ("m/m", 1.0)])
def test_strain_in_percent_and_microstrain_converted(tmp_path, unit, factor):
    rows = _rows([1, 2], [1, 2], strain=True, strain_names=E_COLS)
    path = _write(
        tmp_path / "a.csv",
        rows,
        BASE + S_COLS + E_COLS,
        strain_convention="tensor",
        strain_unit=unit,
    )
    block, meta, _ = _load(path, tmp_path)
    np.testing.assert_allclose(block.E_hist[1, 1], _tensor(lambda k: _eps(2, 2, k), factor))
    assert meta.source_units["strain"] == unit


def test_psi_and_inch_file_roundtrips(tmp_path):
    rows = _rows([1, 2], [1, 2])
    rows = [dict(r, x=r["x"] / 25.4) for r in rows]
    path = _write(tmp_path / "a.csv", rows, BASE + S_COLS, stress_unit="psi", length_unit="in")
    block, meta, _ = _load(path, tmp_path)
    np.testing.assert_allclose(
        block.S_hist[1, 1], _tensor(lambda k: _val(2, 2, k), 6.894757293168361e-3)
    )
    np.testing.assert_allclose(block.coordinates[:, 0], [1.0, 2.0])
    assert meta.source_units == {"stress": "psi", "length": "in"}


def test_non_global_frame_refused(tmp_path):
    path = _write(tmp_path / "a.csv", _rows([1], [1, 2]), BASE + S_COLS, coordinate_system="local")
    with pytest.raises(ValueError, match="coordinate_system"):
        import_fe_csv(path, tmp_path / "h.h5")


def test_elastic_plastic_requires_stabilised_cycle(tmp_path):
    rows = _rows([1, 2], [1, 2], strain=True, strain_names=E_COLS)
    cols = BASE + S_COLS + E_COLS
    kw = dict(strain_convention="tensor", strain_unit="1", analysis_type="elastic_plastic")
    with pytest.raises(ValueError, match="cycle"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, cols, **kw), tmp_path / "h.h5")
    meta = import_fe_csv(
        _write(tmp_path / "b.csv", rows, cols, cycle="stabilised", **kw), tmp_path / "h.h5"
    )
    assert meta.analysis_type == "elastic_plastic" and meta.cycle == "stabilised"


def test_total_strain_is_eel_plus_epl(tmp_path):
    def extra(node, step):
        return {c: 0.5 * _eps(node, step, k) for k, c in enumerate(EPL_COLS)}

    rows = _rows([1, 2], [1, 2], strain=True, strain_names=EEL_COLS, extra=extra)
    path = _write(
        tmp_path / "a.csv", rows, BASE + S_COLS + EEL_COLS + EPL_COLS,
        strain_convention="tensor", strain_unit="1",
    )  # fmt: skip
    block, meta, _ = _load(path, tmp_path)
    np.testing.assert_allclose(block.E_elastic_hist[0, 1], _tensor(lambda k: _eps(1, 2, k)))
    np.testing.assert_allclose(block.E_plastic_hist[0, 1], _tensor(lambda k: _eps(1, 2, k), 0.5))
    np.testing.assert_allclose(block.E_hist[0, 1], _tensor(lambda k: _eps(1, 2, k), 1.5))
    assert meta.strain_kind == "total_mechanical"


def test_total_strain_with_split_part_refused(tmp_path):
    rows = _rows(
        [1, 2],
        [1, 2],
        strain=True,
        strain_names=E_COLS,
        extra=lambda n, s: {c: 0.0 for c in EEL_COLS},
    )
    path = _write(
        tmp_path / "a.csv",
        rows,
        BASE + S_COLS + E_COLS + EEL_COLS,
        strain_convention="tensor",
        strain_unit="1",
    )
    with pytest.raises(ValueError, match="total"):
        import_fe_csv(path, tmp_path / "h.h5")


def test_epl_alone_needs_elastic_constants_and_then_derives_the_elastic_part(tmp_path):
    rows = _rows([1, 2], [1, 2], strain=True, strain_names=EPL_COLS)
    path = _write(
        tmp_path / "a.csv", rows, BASE + S_COLS + EPL_COLS, strain_convention="tensor",
        strain_unit="1", analysis_type="elastic_plastic", cycle="stabilised",
    )  # fmt: skip
    with pytest.raises(ValueError, match="elastic_constants"):
        import_fe_csv(path, tmp_path / "h.h5")
    E, nu = 200000.0, 0.3
    meta = import_fe_csv(path, tmp_path / "h.h5", elastic_constants={"": (E, nu)})
    assert meta.strain_kind == "derived_elastic_part"
    block, *_ = _read(tmp_path / "h.h5")
    s = _tensor(lambda k: _val(2, 2, k))
    by_hand = ((1 + nu) * s - nu * np.trace(s) * np.eye(3)) / E
    np.testing.assert_allclose(block.E_elastic_hist[1, 1], by_hand, rtol=1e-12)
    np.testing.assert_allclose(block.E_plastic_hist[1, 1], _tensor(lambda k: _eps(2, 2, k)))
    np.testing.assert_allclose(
        block.E_hist[1, 1], by_hand + _tensor(lambda k: _eps(2, 2, k)), rtol=1e-12
    )


def test_elastic_constants_must_cover_every_material(tmp_path):
    rows = []
    for step in (1, 2):
        for node, mat in ((1, "steel"), (2, "alu")):
            row = {
                "node_id": node,
                "material": mat,
                "step": step,
                "x": float(node),
                "y": 0.0,
                "z": 0.0,
            }
            row.update({c: _val(node, step, k) for k, c in enumerate(S_COLS)})
            row.update({c: _eps(node, step, k) for k, c in enumerate(EPL_COLS)})
            rows.append(row)
    path = _write(
        tmp_path / "a.csv", rows, ["node_id", "material", "step", "x", "y", "z"] + S_COLS + EPL_COLS,
        strain_convention="tensor", strain_unit="1", analysis_type="elastic_plastic", cycle="stabilised",
    )  # fmt: skip
    with pytest.raises(ValueError, match="alu"):
        import_fe_csv(path, tmp_path / "h.h5", elastic_constants={"steel": (210000.0, 0.3)})
    import_fe_csv(
        path,
        tmp_path / "h.h5",
        elastic_constants={"steel": (210000.0, 0.3), "alu": (70000.0, 0.33)},
    )
    block, *_ = _read(tmp_path / "h.h5")
    s = _tensor(lambda k: _val(2, 1, k))
    want = ((1 + 0.33) * s - 0.33 * np.trace(s) * np.eye(3)) / 70000.0
    np.testing.assert_allclose(block.E_elastic_hist[1, 0], want, rtol=1e-12)


def test_eel_alone_is_total_strain_for_linear_elastic_and_refused_under_plasticity(tmp_path):
    rows = _rows([1, 2], [1, 2], strain=True, strain_names=EEL_COLS)
    cols = BASE + S_COLS + EEL_COLS
    ok = _write(tmp_path / "a.csv", rows, cols, strain_convention="tensor", strain_unit="1")
    meta = import_fe_csv(ok, tmp_path / "h.h5")
    block, *_ = _read(tmp_path / "h.h5")
    np.testing.assert_allclose(block.E_hist[0, 0], _tensor(lambda k: _eps(1, 1, k)))
    assert block.E_elastic_hist is None and meta.strain_kind == "total_mechanical"
    bad = _write(
        tmp_path / "b.csv", rows, cols, strain_convention="tensor", strain_unit="1",
        analysis_type="elastic_plastic", cycle="stabilised",
    )  # fmt: skip
    with pytest.raises(ValueError, match="Eel"):
        import_fe_csv(bad, tmp_path / "h2.h5")


def test_linear_elastic_without_strain_derives_only_with_constants_and_is_flagged(tmp_path):
    path = _write(tmp_path / "a.csv", _rows([1, 2], [1, 2]), BASE + S_COLS)
    meta = import_fe_csv(path, tmp_path / "h.h5")
    assert meta.strain_kind is None
    with HistoryFile(tmp_path / "h.h5") as f:
        assert not f.has_strain
    E, nu = 70000.0, 0.33
    meta = import_fe_csv(path, tmp_path / "d.h5", elastic_constants={"": (E, nu)})
    assert meta.strain_kind == "derived_linear_elastic"
    block, *_ = _read(tmp_path / "d.h5")
    s = _tensor(lambda k: _val(1, 2, k))
    np.testing.assert_allclose(
        block.E_hist[0, 1], ((1 + nu) * s - nu * np.trace(s) * np.eye(3)) / E, rtol=1e-12
    )


def test_elastic_plastic_without_any_strain_refused(tmp_path):
    path = _write(
        tmp_path / "a.csv", _rows([1, 2], [1, 2]), BASE + S_COLS, analysis_type="elastic_plastic",
        cycle="stabilised",
    )  # fmt: skip
    with pytest.raises(ValueError, match="strain"):
        import_fe_csv(path, tmp_path / "h.h5", elastic_constants={"": (200000.0, 0.3)})


def test_global_average_across_materials_refused_unless_allowed_and_recorded(tmp_path):
    rows = []
    for step in (1, 2):
        for node, mat in ((1, "steel"), (2, "alu")):
            row = {
                "node_id": node,
                "material": mat,
                "step": step,
                "x": float(node),
                "y": 0.0,
                "z": 0.0,
            }
            row.update({c: _val(node, step, k) for k, c in enumerate(S_COLS)})
            rows.append(row)
    cols = ["node_id", "material", "step", "x", "y", "z"] + S_COLS
    path = _write(tmp_path / "a.csv", rows, cols, averaging="nodal_average_global")
    with pytest.raises(ValueError, match="allow_global_average"):
        import_fe_csv(path, tmp_path / "h.h5")
    meta = import_fe_csv(path, tmp_path / "h.h5", allow_global_average=True)
    assert meta.averaging == "nodal_average_global"
    assert meta.source_conventions["global_average_allowed"] == "true"
    one = _write(
        tmp_path / "b.csv", _rows([1, 2], [1, 2]), BASE + S_COLS, averaging="nodal_average_global"
    )
    import_fe_csv(one, tmp_path / "h1.h5")  # a single material and instance needs no flag


def test_unknown_column_refused_unless_ignored(tmp_path):
    rows = _rows([1, 2], [1, 2], extra=lambda n, s: {"Volume": 3.0})
    path = _write(tmp_path / "a.csv", rows, BASE + S_COLS + ["Volume"])
    with pytest.raises(ValueError, match="Volume"):
        import_fe_csv(path, tmp_path / "h.h5")
    import_fe_csv(path, tmp_path / "h.h5", ignore_columns=["volume"])


def test_temperature_column_roundtrips(tmp_path):
    rows = _rows([1, 2], [1, 2], extra=lambda n, s: {"T": 20.0 * n + s})
    path = _write(tmp_path / "a.csv", rows, BASE + S_COLS + ["T"], temperature_unit="C")
    meta = import_fe_csv(path, tmp_path / "h.h5")
    block, *_ = _read(tmp_path / "h.h5")
    np.testing.assert_allclose(block.T_hist, [[21.0, 22.0], [41.0, 42.0]])
    assert meta.source_units["temperature"] == "C"


def test_node_major_row_limit_gives_advice(tmp_path):
    path = write_tube_csv(
        tmp_path / "n.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=3, strain="none",
        order="node_major",
    )  # fmt: skip
    with pytest.raises(ValueError, match="step-major"):
        import_fe_csv(path, tmp_path / "h.h5", max_node_major_rows=10)
    assert not (tmp_path / "h.h5").exists()


def test_step_blocks_spanning_several_read_chunks(tmp_path):
    path = write_tube_csv(
        tmp_path / "t.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=3, strain="total"
    )
    import_fe_csv(path, tmp_path / "a.h5", chunk_rows=4)  # chunks cut inside and across steps
    import_fe_csv(path, tmp_path / "b.h5")
    a, *_ = _read(tmp_path / "a.h5")
    b, *_ = _read(tmp_path / "b.h5")
    np.testing.assert_array_equal(a.S_hist, b.S_hist)
    np.testing.assert_array_equal(a.E_hist, b.E_hist)


def test_rows_in_a_different_order_within_each_step_are_matched_by_key(tmp_path):
    rows = _rows([1, 2, 3], [1, 2])
    shuffled = rows[:3] + rows[3:][::-1]  # step 2 lists the nodes in reverse
    block, *_ = _load(_write(tmp_path / "a.csv", shuffled, BASE + S_COLS), tmp_path)
    np.testing.assert_allclose(block.S_hist[0, 1], _tensor(lambda k: _val(1, 2, k)))
    np.testing.assert_allclose(block.S_hist[2, 1], _tensor(lambda k: _val(3, 2, k)))


def test_coordinates_that_change_between_steps_refused(tmp_path):
    rows = _rows([1, 2], [1, 2])
    rows[-1] = dict(rows[-1], x=rows[-1]["x"] + 0.5)
    with pytest.raises(ValueError, match="coordinates"):
        import_fe_csv(_write(tmp_path / "a.csv", rows, BASE + S_COLS), tmp_path / "h.h5")


def test_provenance_records_hash_header_and_alias_mapping(tmp_path):
    path = write_tube_csv(
        tmp_path / "t.csv", GEOM, MAT, LOADS, n_radii=2, n_theta=3, flavour="abaqus",
        shear="engineering", strain="split", analysis_type="elastic_plastic",
        header={"notes": "provenance test"},
    )  # fmt: skip
    meta = import_fe_csv(path, tmp_path / "h.h5")
    assert [(s.name, s.sha256) for s in meta.source_files] == [("t.csv", sha256_of(path))]
    assert meta.solver.startswith("analytical thick tube") and meta.solver_version == "1"
    assert meta.source_frame == "global_cartesian" and meta.averaging == "none"
    assert meta.analysis_type == "elastic_plastic" and meta.cycle == "stabilised"
    assert meta.step_role == "time_sample" and meta.notes == "provenance test"
    assert meta.source_units["strain"] == "1"
    conv = meta.source_conventions
    assert conv["strain_convention"] == "engineering"
    assert conv["component_order"] == "xx,yy,zz,xy,yz,xz"
    aliases = json.loads(conv["column_aliases_used"])
    assert aliases["S_xz"] == "S13" and aliases["Eel_yz"] == "EE23" and aliases["Epl_xy"] == "PE12"
    assert meta.adapter == "multifat.fe.csvio"


def _peak(tmp_path, n_steps, n_nodes=1500):
    loads = [TubeLoad(float(k), 100.0 * k, 10.0 * k) for k in range(n_steps)]
    path = write_tube_csv(
        tmp_path / f"m{n_steps}.csv",
        GEOM,
        MAT,
        loads,
        n_radii=n_nodes // 3,
        n_theta=3,
        strain="total",
    )
    tracemalloc.start()
    try:
        import_fe_csv(path, tmp_path / f"m{n_steps}.h5", chunk_rows=500)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def test_streaming_peak_memory_independent_of_the_number_of_steps(tmp_path):
    few, many = _peak(tmp_path, 3), _peak(tmp_path, 30)
    assert many < 1.5 * few + 1_000_000, (few, many)


EXAMPLE = (
    Path(__file__).resolve().parents[1] / "docs" / "examples" / "fe-csv" / "tube_two_steps.csv"
)


def test_committed_example_equals_the_generated_one(tmp_path):
    assert main(["--csv", str(tmp_path)]) == 0
    generated = (tmp_path / "tube_two_steps.csv").read_bytes()
    committed = EXAMPLE.read_bytes().replace(b"\r\n", b"\n")  # a Windows checkout may use CRLF
    assert committed == generated


def test_committed_example_imports_and_matches_the_oracle(tmp_path):
    meta = import_fe_csv(EXAMPLE, tmp_path / "h.h5")
    ids, xyz, r, theta = tube_points(GEOM, 2, 2)
    two_steps = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 12000.0, 8e4)]
    block, _, cases = _read(tmp_path / "h.h5")
    for k, load in enumerate(two_steps):
        S, E = tube_cartesian(r, theta, load, GEOM, MAT)
        _close(block.S_hist[:, k], S, scale=1e3)
        _close(block.E_hist[:, k], E)
    assert list(cases.steps) == [1, 2] and block.n_locations == 4
    assert meta.solver.startswith("analytical thick tube")
