"""Node assessment: node histories, unit-case superposition, statuses and chunking."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from multifat.fe.csvio import import_fe_csv
from multifat.fe.csvwriter import FeCsvWriter
from multifat.fe.model import HistoryMeta, LoadCases
from multifat.fe.store import HistoryWriter
from multifat.fe.synthetic import (
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    tube_cartesian,
    tube_points,
    write_tube_csv,
)
from multifat.fe.tensors import elastic_strain
from multifat.history.accumulate import SNCurve
from multifat.history.assess import assess_history
from multifat.history.eligibility import IneligibleHistoryError
from multifat.history.nodes import (
    MaterialSetup,
    assess_nodes,
    node_histories,
    read_load_functions,
)
from multifat.history.signal import TimeHistory
from multifat.models.registry import _PROBE_PARAMS, declared_param_names, get_models
from multifat.planes.sampling import build_R_list

pytest.importorskip("h5py")

GEOM, MAT = TubeGeometry(a=6.0, b=7.5), TubeMaterial(E=210000.0, nu=0.3)
E_MOD, NU = MAT.E, MAT.nu
R_LIST = build_R_list(n1=3, n2=6)
N_RADII, N_THETA = 2, 3
(MODEL,) = get_models(["FS"])
PARAMS = {
    k: _PROBE_PARAMS.get(k, 1.0) for k in set(declared_param_names()) | set(MODEL.required_params())
}
PARAMS.update(E=E_MOD, nu=NU)
CURVE = SNCurve(A=2000.0, b=-0.12, method="FS", variant="")
SETUP = MaterialSetup(MODEL, PARAMS, CURVE)

LAMBDAS = [0.0, 1.0, -0.5, 0.8, -1.0, 0.3, 0.0]  # a proportional load path with several reversals


def _tube_loads(lams, scale=(50.0, 12000.0, 8e4)):
    return [TubeLoad(scale[0] * x, scale[1] * x, scale[2] * x) for x in lams]


def _tube_file(tmp_path, lams=LAMBDAS, name="t", **kw):
    csv = write_tube_csv(
        tmp_path / f"{name}.csv", GEOM, MAT, _tube_loads(lams), n_radii=N_RADII,
        n_theta=N_THETA, strain="total", **kw,
    )  # fmt: skip
    import_fe_csv(csv, tmp_path / f"{name}.h5")
    return tmp_path / f"{name}.h5"


def _oracle(lams):
    ids, xyz, r, theta = tube_points(GEOM, N_RADII, N_THETA)
    S = np.empty((len(ids), len(lams), 3, 3))
    E = np.empty_like(S)
    for k, load in enumerate(_tube_loads(lams)):
        S[:, k], E[:, k] = tube_cartesian(r, theta, load, GEOM, MAT)
    return S, E


def _direct(S, E, axis, **kw):
    h = TimeHistory.from_arrays(axis, S, E, stress_unit="MPa", axis_kind="cycle", time_unit="step")
    return assess_history(h, MODEL, PARAMS, CURVE, R_list=R_LIST, **kw)


def test_two_steps_give_one_cycle_and_the_constant_amplitude_value(tmp_path):
    path = _tube_file(tmp_path, [-1.0, 1.0])
    S, E = _oracle([-1.0, 1.0])
    rows = assess_nodes(path, SETUP, R_list=R_LIST)
    assert len(rows) == 6
    for n, row in enumerate(rows):
        assert row["status"] == "ok" and row["n_events"] == 1
        a, b = S[n, 0], S[n, 1]
        ea, eb = E[n, 0], E[n, 1]
        harmonics = {
            "S_mean": 0.5 * (a + b), "S_sin": 0.5 * (a - b), "S_cos": np.zeros((3, 3)),
            "E_mean": 0.5 * (ea + eb), "E_sin": 0.5 * (ea - eb), "E_cos": np.zeros((3, 3)),
        }  # fmt: skip
        dp = MODEL.evaluate_case(
            S0=a, E0=ea, S1=b, E1=eb, params=PARAMS, R_list=R_LIST, harmonics=harmonics
        ).values[""]
        n_f = CURVE.cycles_to_failure(np.array([dp]))[0]
        assert row["damage_per_block"] == pytest.approx(1.0 / n_f, rel=1e-9)
        assert row["blocks_to_failure"] == pytest.approx(n_f, rel=1e-9)
        direct = _direct(S[n], E[n], [1.0, 2.0])
        assert row["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-12)


def test_time_sample_history_equals_direct_assess_history_for_every_node(tmp_path):
    path = _tube_file(tmp_path)
    S, E = _oracle(LAMBDAS)
    rows = assess_nodes(path, SETUP, R_list=R_LIST)
    assert [r["location_id"] for r in rows] == [1, 2, 3, 4, 5, 6]
    for n, row in enumerate(rows):
        direct = _direct(S[n], E[n], np.arange(1.0, len(LAMBDAS) + 1))
        assert row["status"] == "ok"
        assert row["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-9)
        assert row["n_events"] == len(direct.events)
        assert row["cycles_to_failure"] == pytest.approx(direct.cycles_to_failure, rel=1e-9)
        assert row["instance"] == "" and row["material"] == ""


def test_axis_is_time_when_present_else_step_index_as_cycle(tmp_path):
    with_time = _tube_file(tmp_path, name="a", with_time=True)
    without = _tube_file(tmp_path, name="b")
    ref, h = next(node_histories(with_time))
    assert h.axis_kind == "time" and np.allclose(h.time, np.linspace(0.0, 1.0, len(LAMBDAS)))
    assert h.units["time"] == "s"
    ref, h = next(node_histories(without))
    assert h.axis_kind == "cycle" and np.array_equal(h.time, np.arange(1.0, len(LAMBDAS) + 1))
    assert h.units["time"] == "step"
    assert h.provenance["source"] == "fe_node" and h.provenance["location_id"] == ref.location_id
    assert h.provenance["source_files"][0]["name"] == "b.csv"


def _unit_case_file(tmp_path, name="u"):
    """Three unit cases (pressure, axial force, torque), one step each."""
    unit = [TubeLoad(1.0, 0.0, 0.0), TubeLoad(0.0, 1.0, 0.0), TubeLoad(0.0, 0.0, 1.0)]
    csv = write_tube_csv(
        tmp_path / f"{name}.csv", GEOM, MAT, unit, n_radii=N_RADII, n_theta=N_THETA,
        strain="total", header={"step_role": "unit_case"},
    )  # fmt: skip
    import_fe_csv(csv, tmp_path / f"{name}.h5")
    return tmp_path / f"{name}.h5"


def _functions_csv(tmp_path, lams=LAMBDAS, name="f.csv", **override):
    t = np.arange(len(lams), dtype=float)
    frame = pd.DataFrame(
        {"time": t, "f_1": 50.0 * np.array(lams), "f_2": 12000.0 * np.array(lams),
         "f_3": 8e4 * np.array(lams)}  # fmt: skip
    )
    for key, value in override.items():
        frame[key] = value
    frame.to_csv(tmp_path / name, index=False)
    return tmp_path / name


def test_unit_case_superposition_equals_the_directly_built_history(tmp_path):
    path = _unit_case_file(tmp_path)
    functions = _functions_csv(tmp_path)
    S, E = _oracle(LAMBDAS)
    got = list(node_histories(path, load_functions=functions))
    assert len(got) == 6
    for n, (ref, history) in enumerate(got):
        np.testing.assert_allclose(history.S_hist, S[n], rtol=1e-11, atol=1e-9)
        np.testing.assert_allclose(history.E_hist, E[n], rtol=1e-11, atol=1e-14)
        assert history.axis_kind == "time" and history.provenance["step_role"] == "unit_case"
        assert history.provenance["load_functions"] == {"file": "f.csv"}
    rows = assess_nodes(path, SETUP, R_list=R_LIST, load_functions=functions)
    for n, row in enumerate(rows):
        direct = _direct(S[n], E[n], np.arange(len(LAMBDAS), dtype=float))
        assert row["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-9)


def test_unit_case_file_needs_functions_and_a_time_sample_file_refuses_them(tmp_path):
    unit = _unit_case_file(tmp_path)
    sample = _tube_file(tmp_path, name="s")
    with pytest.raises(ValueError, match="load_functions"):
        next(node_histories(unit))
    with pytest.raises(ValueError, match="unit cases"):
        next(node_histories(sample, load_functions=_functions_csv(tmp_path)))
    one_step = tmp_path / "one.csv"
    write_tube_csv(one_step, GEOM, MAT, _tube_loads([1.0]), n_radii=2, n_theta=2, strain="total")
    import_fe_csv(one_step, tmp_path / "one.h5")
    with pytest.raises(ValueError, match="two or more"):
        next(node_histories(tmp_path / "one.h5"))


def test_unit_case_on_elastic_plastic_refused(tmp_path):
    csv = write_tube_csv(
        tmp_path / "ep.csv", GEOM, MAT, [TubeLoad(1.0, 0.0, 0.0), TubeLoad(0.0, 1.0, 0.0)],
        n_radii=2, n_theta=2, strain="split", analysis_type="elastic_plastic",
        header={"step_role": "unit_case"},
    )  # fmt: skip
    with pytest.raises(ValueError, match="unit_case"):  # refused when importing
        import_fe_csv(csv, tmp_path / "ep.h5")
    # a file written by another tool: declared elastic-plastic, step_role unit_case
    meta = HistoryMeta(
        material_id="m", solver="x", solver_version="1", analysis_step="s", adapter="t",
        adapter_version="0", source_files=(), source_units={"stress": "MPa", "length": "mm"},
        source_frame="global_cartesian", averaging="none", location_kind="node",
        analysis_type="elastic_plastic", step_role="unit_case",
    )  # fmt: skip
    z = np.zeros((1, 3, 3))
    with HistoryWriter(
        tmp_path / "x.h5", meta, np.array([1]), np.zeros((1, 3)),
        LoadCases(steps=np.array([1, 2])), 2, with_strain=False,
    ) as w:  # fmt: skip
        w.write_step(0, z, None)
        w.write_step(1, z, None)
    with pytest.raises(ValueError, match="linear-elastic"):
        next(
            node_histories(
                tmp_path / "x.h5", load_functions=(np.array([0.0, 1.0]), np.ones((2, 2)), "time")
            )
        )


def test_load_function_columns_matched_by_name_missing_and_extra_refused(tmp_path):
    steps = [1, 2, 3]
    ok = _functions_csv(tmp_path)
    axis, F, kind = read_load_functions(ok, steps)
    assert kind == "time" and F.shape == (len(LAMBDAS), 3)
    shuffled = pd.read_csv(ok)[["f_3", "time", "f_1", "f_2"]]
    shuffled.to_csv(tmp_path / "s.csv", index=False)
    axis2, F2, _ = read_load_functions(tmp_path / "s.csv", steps)
    np.testing.assert_array_equal(F, F2)
    pd.read_csv(ok).drop(columns=["f_2"]).to_csv(tmp_path / "m.csv", index=False)
    with pytest.raises(ValueError, match="f_2"):
        read_load_functions(tmp_path / "m.csv", steps)
    with pytest.raises(ValueError, match="f_9"):
        read_load_functions(_functions_csv(tmp_path, name="e.csv", f_9=0.0), steps)
    with pytest.raises(ValueError, match="Axis|axis"):
        pd.read_csv(ok).drop(columns=["time"]).to_csv(tmp_path / "n.csv", index=False)
        read_load_functions(tmp_path / "n.csv", steps)
    with pytest.raises(ValueError, match="axis"):
        read_load_functions(_functions_csv(tmp_path, name="b.csv", index=0.0), steps)
    index_axis = pd.read_csv(ok).rename(columns={"time": "index"})
    index_axis.to_csv(tmp_path / "i.csv", index=False)
    assert read_load_functions(tmp_path / "i.csv", steps)[2] == "cycle"
    bad = pd.read_csv(ok)
    bad.loc[2, "time"] = bad.loc[1, "time"]
    bad.to_csv(tmp_path / "d.csv", index=False)
    with pytest.raises(ValueError, match="increasing"):
        read_load_functions(tmp_path / "d.csv", steps)


# --- hand-made nodes: proportional and non-proportional histories, materials and instances ------


def _node_history(kind, scale=1.0):
    lam = np.array(LAMBDAS)
    other = np.array([0.0, 0.2, 1.0, -0.4, 0.5, -1.0, 0.0])  # a different shape: not proportional
    S = np.zeros((len(lam), 3, 3))
    S[:, 2, 2] = 30.0 + scale * 100.0 * lam
    shape = lam if kind == "proportional" else other
    S[:, 0, 2] = S[:, 2, 0] = scale * 40.0 * shape
    return S


def _craft(tmp_path, nodes, name="c", constants=None):
    """Import a CSV of hand-made nodes: ``nodes`` are ``(instance, node, material, S)``."""
    columns = ["node_id", "instance", "material", "step", "x", "y", "z"]
    columns += [f"S_{c}" for c in ("xx", "yy", "zz", "xy", "yz", "xz")]
    header = {
        "format": "multifat-fe-csv 1.0", "solver": "unit test", "solver_version": "1",
        "stress_unit": "MPa", "length_unit": "mm", "coordinate_system": "global_cartesian",
        "averaging": "none", "analysis_type": "linear_elastic", "step_role": "time_sample",
    }  # fmt: skip
    n_steps = nodes[0][3].shape[0]
    with FeCsvWriter(tmp_path / f"{name}.csv", header, columns=columns) as w:
        for step in range(n_steps):
            for n, (inst, node, material, S) in enumerate(nodes):
                row = {"node_id": node, "instance": inst, "material": material, "step": step + 1,
                       "x": float(n), "y": 0.0, "z": 0.0}  # fmt: skip
                for c, (i, j) in zip(
                    ("xx", "yy", "zz", "xy", "yz", "xz"),
                    ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)),
                ):
                    row[f"S_{c}"] = float(S[step, i, j])
                w.write_row(row)
    constants = constants or {"": (E_MOD, NU)}
    import_fe_csv(tmp_path / f"{name}.csv", tmp_path / f"{name}.h5", elastic_constants=constants)
    return tmp_path / f"{name}.h5"


def _history_of(S, material="", provenance=None):
    return TimeHistory.from_arrays(
        np.arange(1.0, len(S) + 1), S, elastic_strain(S, E_MOD, NU), stress_unit="MPa",
        axis_kind="cycle", time_unit="step",
    )  # fmt: skip


def test_per_node_status_records_ineligible_and_continues(tmp_path):
    nodes = [("", 1, "", _node_history("proportional")), ("", 2, "", _node_history("other")),
             ("", 3, "", _node_history("proportional", 1.5))]  # fmt: skip
    path = _craft(tmp_path, nodes)
    rows = assess_nodes(path, SETUP, R_list=R_LIST)
    assert [r["status"][:10] for r in rows] == ["ok", "ineligible", "ok"]
    assert "proportional" in rows[1]["status"]
    assert rows[1]["damage_per_block"] is None and rows[1]["n_events"] is None
    assert (
        rows[0]["damage_per_block"] > 0.0
        and rows[2]["damage_per_block"] > rows[0]["damage_per_block"]
    )


def test_on_refusal_raise_aborts_at_the_first(tmp_path):
    nodes = [("", 1, "", _node_history("other")), ("", 2, "", _node_history("proportional"))]
    path = _craft(tmp_path, nodes)
    with pytest.raises(IneligibleHistoryError):
        assess_nodes(path, SETUP, R_list=R_LIST, on_refusal="raise")
    with pytest.raises(ValueError, match="on_refusal"):
        assess_nodes(path, SETUP, R_list=R_LIST, on_refusal="skip")


def test_nonproportional_node_refused_then_allowed_with_the_pass_through_flag(tmp_path):
    path = _craft(tmp_path, [("", 1, "", _node_history("other"))])
    (refused,) = assess_nodes(path, SETUP, R_list=R_LIST)
    assert refused["status"].startswith("ineligible")
    (allowed,) = assess_nodes(path, SETUP, R_list=R_LIST, allow_nonproportional_approximation=True)
    assert allowed["status"] == "ok" and allowed["approximation"]
    direct = assess_history(
        _history_of(_node_history("other")), MODEL, PARAMS, CURVE, R_list=R_LIST,
        allow_nonproportional_approximation=True,
    )  # fmt: skip
    assert allowed["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-9)


def test_material_without_setup_gets_status_and_run_continues(tmp_path):
    nodes = [("", 1, "steel", _node_history("proportional")),
             ("", 2, "alu", _node_history("proportional", 0.5))]  # fmt: skip
    path = _craft(tmp_path, nodes, constants={"steel": (E_MOD, NU), "alu": (70000.0, 0.33)})
    rows = assess_nodes(path, {"steel": SETUP}, R_list=R_LIST)
    assert [r["status"] for r in rows] == ["ok", "no_material_setup"]
    assert rows[1]["material"] == "alu" and rows[1]["damage_per_block"] is None
    with pytest.raises(ValueError, match="alu"):
        assess_nodes(path, {"steel": SETUP}, R_list=R_LIST, on_refusal="raise")


def test_interface_node_assessed_once_per_material_and_repeated_ids_across_instances_distinct(
    tmp_path,
):
    nodes = [
        ("A", 1, "steel", _node_history("proportional", 1.0)),
        ("B", 1, "steel", _node_history("proportional", 1.3)),
        ("A", 2, "steel", _node_history("proportional", 0.9)),
        ("A", 2, "alu", _node_history("proportional", 0.6)),
    ]
    path = _craft(tmp_path, nodes, constants={"steel": (E_MOD, NU), "alu": (E_MOD, NU)})
    setups = {"steel": SETUP, "alu": SETUP}
    rows = assess_nodes(path, setups, R_list=R_LIST)
    assert [(r["instance"], r["location_id"], r["material"]) for r in rows] == [
        ("A", 1, "steel"), ("B", 1, "steel"), ("A", 2, "steel"), ("A", 2, "alu"),
    ]  # fmt: skip
    for row, (inst, node, material, S) in zip(rows, nodes):
        direct = assess_history(_history_of(S), MODEL, PARAMS, CURVE, R_list=R_LIST)
        assert row["status"] == "ok"
        assert row["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-9)
    assert len({r["damage_per_block"] for r in rows}) == 4  # four distinct locations


def test_stress_only_file_gets_hooke_strain_from_the_setup_and_is_flagged(tmp_path):
    S = _node_history("proportional")
    columns = ["node_id", "step", "x", "y", "z"] + [
        f"S_{c}" for c in ("xx", "yy", "zz", "xy", "yz", "xz")
    ]
    header = {
        "format": "multifat-fe-csv 1.0", "solver": "unit test", "solver_version": "1",
        "stress_unit": "MPa", "length_unit": "mm", "coordinate_system": "global_cartesian",
        "averaging": "none", "analysis_type": "linear_elastic", "step_role": "time_sample",
    }  # fmt: skip
    with FeCsvWriter(tmp_path / "s.csv", header, columns=columns) as w:
        for step in range(len(S)):
            row = {"node_id": 1, "step": step + 1, "x": 0.0, "y": 0.0, "z": 0.0}
            for c, (i, j) in zip(
                ("xx", "yy", "zz", "xy", "yz", "xz"),
                ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)),
            ):
                row[f"S_{c}"] = float(S[step, i, j])
            w.write_row(row)
    import_fe_csv(tmp_path / "s.csv", tmp_path / "s.h5")  # no strain in the file
    (without,) = assess_nodes(tmp_path / "s.h5", SETUP, R_list=R_LIST)
    assert without["status"].startswith("ineligible")
    with_hooke = MaterialSetup(MODEL, PARAMS, CURVE, E=E_MOD, nu=NU)
    (row,) = assess_nodes(tmp_path / "s.h5", with_hooke, R_list=R_LIST)
    direct = assess_history(_history_of(S), MODEL, PARAMS, CURVE, R_list=R_LIST)
    assert row["status"] == "ok"
    assert row["damage_per_block"] == pytest.approx(direct.damage_per_block, rel=1e-9)
    with pytest.raises(ValueError, match="both E and nu"):
        MaterialSetup(MODEL, PARAMS, CURVE, E=E_MOD)


def test_positions_subset_and_chunking_do_not_change_results(tmp_path):
    path = _tube_file(tmp_path)
    full = assess_nodes(path, SETUP, R_list=R_LIST)
    chunked = assess_nodes(path, SETUP, R_list=R_LIST, chunk_locations=1)
    assert chunked == full
    subset = assess_nodes(path, SETUP, R_list=R_LIST, positions=[4, 1], chunk_locations=1)
    assert [r["position"] for r in subset] == [4, 1]
    assert subset[0] == full[4] and subset[1] == full[1]
    assert [r.position for r, _ in node_histories(path, positions=[5, 0], chunk_locations=2)] == [
        5,
        0,
    ]
