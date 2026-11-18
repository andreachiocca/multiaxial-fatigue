"""Time-history contract: validation, constructors, reproducible record."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from multifat.fe.model import HistoryBlock
from multifat.history.signal import TimeHistory


def _sym(rng, shape):
    a = rng.normal(size=shape + (3, 3)) * 50.0
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def _arrays(n=8, seed=0, strain=False):
    rng = np.random.default_rng(seed)
    time = np.linspace(0.0, 1.0, n)
    S = _sym(rng, (n,))
    E = _sym(rng, (n,)) * 1e-5 if strain else None
    return time, S, E


def test_rejects_non_increasing_or_duplicate_times():
    time, S, _ = _arrays()
    bad = time.copy()
    bad[3] = bad[2]
    with pytest.raises(ValueError, match="increasing"):
        TimeHistory.from_arrays(bad, S, stress_unit="MPa")
    bad = time.copy()
    bad[4], bad[5] = bad[5], bad[4]
    with pytest.raises(ValueError, match="increasing"):
        TimeHistory.from_arrays(bad, S, stress_unit="MPa")


def test_rejects_nan_and_asymmetric_tensors():
    time, S, _ = _arrays()
    nan = S.copy()
    nan[2, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        TimeHistory.from_arrays(time, nan, stress_unit="MPa")
    asym = S.copy()
    asym[1, 0, 1] += 10.0
    with pytest.raises(ValueError, match="symmetric"):
        TimeHistory.from_arrays(time, asym, stress_unit="MPa")
    with pytest.raises(ValueError, match="finite"):
        TimeHistory.from_arrays(np.where(np.arange(8) == 3, np.nan, time), S, stress_unit="MPa")


def test_rejects_single_sample_and_shape_mismatch():
    time, S, _ = _arrays()
    with pytest.raises(ValueError, match="at least two"):
        TimeHistory.from_arrays(time[:1], S[:1], stress_unit="MPa")
    with pytest.raises(ValueError, match="shape"):
        TimeHistory.from_arrays(time, S[:-1], stress_unit="MPa")
    with pytest.raises(ValueError, match="shape"):
        TimeHistory.from_arrays(time, S, _arrays(strain=True)[2][:-1], stress_unit="MPa")


def test_from_arrays_requires_stress_unit_and_knows_units():
    time, S, _ = _arrays()
    with pytest.raises(TypeError):
        TimeHistory.from_arrays(time, S)  # type: ignore[call-arg]
    with pytest.raises(ValueError, match="unit"):
        TimeHistory.from_arrays(time, S, stress_unit="furlong")


def test_stress_unit_pa_is_converted_to_mpa():
    time, S, _ = _arrays()
    h = TimeHistory.from_arrays(time, S * 1e6, stress_unit="Pa")
    assert np.allclose(h.S_hist, S) and h.units["stress"] == "MPa"
    assert h.units["source_stress_unit"] == "Pa"


def test_from_csv_roundtrip_and_missing_unit_rejected(tmp_path):
    time, S, E = _arrays(strain=True)
    cols = {"t": time}
    for name, (i, j) in zip(
        ("xx", "yy", "zz", "xy", "yz", "xz"), ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))
    ):
        cols["s" + name] = S[:, i, j]
        cols["e" + name] = 2.0 * E[:, i, j] if i != j else E[:, i, j]  # engineering shear
    path = tmp_path / "h.csv"
    pd.DataFrame(cols).to_csv(path, index=False)
    h = TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s", strain_shear="engineering")
    assert np.allclose(h.time, time) and np.allclose(h.S_hist, S) and np.allclose(h.E_hist, E)
    assert h.provenance["source"] == "csv" and len(h.provenance["sha256"]) == 64
    with pytest.raises(TypeError):
        TimeHistory.from_csv(path, stress_unit="MPa")  # type: ignore[call-arg]
    with pytest.raises(TypeError):
        TimeHistory.from_csv(path, time_unit="s")  # type: ignore[call-arg]


def test_from_csv_custom_columns_and_missing_column_rejected(tmp_path):
    time, S, _ = _arrays()
    names = {
        "time": "time_s",
        "xx": "Sx",
        "yy": "Sy",
        "zz": "Sz",
        "xy": "Txy",
        "yz": "Tyz",
        "xz": "Txz",
    }
    frame = pd.DataFrame({"time_s": time, "Sx": S[:, 0, 0], "Sy": S[:, 1, 1], "Sz": S[:, 2, 2],
                          "Txy": S[:, 0, 1], "Tyz": S[:, 1, 2], "Txz": S[:, 0, 2]})  # fmt: skip
    path = tmp_path / "h.csv"
    frame.to_csv(path, index=False, sep=";")
    h = TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s", delimiter=";", columns=names)
    assert np.allclose(h.S_hist, S)
    frame.drop(columns=["Tyz"]).to_csv(path, index=False, sep=";")
    with pytest.raises(ValueError, match="Tyz"):
        TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s", delimiter=";", columns=names)


def _tension_torsion(time, *, amplitude=100.0, shear=40.0, mean=20.0):
    S = np.zeros((len(time), 3, 3))
    lam = np.sin(2.0 * np.pi * time)
    S[:, 2, 2] = mean + amplitude * lam
    S[:, 0, 2] = S[:, 2, 0] = shear * lam
    return S


def test_from_function_records_name_and_params():
    time = np.linspace(0.0, 2.0, 41)
    h = TimeHistory.from_function(_tension_torsion, time, amplitude=80.0, shear=30.0)
    assert np.allclose(h.S_hist, _tension_torsion(time, amplitude=80.0, shear=30.0))
    assert h.provenance["source"] == "function"
    assert h.provenance["function"].endswith("_tension_torsion")
    assert h.provenance["params"] == {"amplitude": 80.0, "shear": 30.0}
    named = TimeHistory.from_function(_tension_torsion, time, name="block A")
    assert named.provenance["function"] == "block A"


def test_from_function_accepts_a_dict_with_strain():
    time = np.linspace(0.0, 1.0, 11)

    def both(t):
        S = _tension_torsion(t)
        return {"S": S, "E": S / 2.1e5}

    h = TimeHistory.from_function(both, time)
    assert h.E_hist is not None and h.E_hist.shape == (11, 3, 3)


def test_from_function_wrong_shape_rejected():
    time = np.linspace(0.0, 1.0, 11)
    with pytest.raises(ValueError, match="shape"):
        TimeHistory.from_function(lambda t: np.zeros((len(t), 3)), time)
    with pytest.raises(ValueError, match="shape"):
        TimeHistory.from_function(lambda t: {"E": np.zeros((len(t), 3, 3))}, time)


def test_from_function_non_json_params_rejected():
    time = np.linspace(0.0, 1.0, 11)
    with pytest.raises(ValueError, match="JSON"):
        TimeHistory.from_function(lambda t, f: _tension_torsion(t), time, f=object())


def test_record_roundtrip_reproduces_arrays_and_hash(tmp_path):
    time, S, E = _arrays(strain=True)
    h = TimeHistory.from_arrays(time, S, E, stress_unit="MPa", provenance={"note": "test"})
    path = h.record(tmp_path / "rec")
    assert path.suffix == ".json" and (tmp_path / "rec.npz").exists()
    back = TimeHistory.from_record(tmp_path / "rec")
    assert back.sha256() == h.sha256() and back.provenance == h.provenance
    assert np.array_equal(back.S_hist, h.S_hist) and np.array_equal(back.E_hist, h.E_hist)
    tampered = np.load(tmp_path / "rec.npz")
    np.savez(tmp_path / "rec.npz", **{k: tampered[k] * (1.0 + 1e-9) for k in tampered.files})
    with pytest.raises(ValueError, match="hash"):
        TimeHistory.from_record(tmp_path / "rec")


def test_sha256_changes_with_any_array():
    time, S, _ = _arrays()
    a = TimeHistory.from_arrays(time, S, stress_unit="MPa")
    S2 = S.copy()
    S2[3, 1, 1] += 1e-6
    assert a.sha256() != TimeHistory.from_arrays(time, S2, stress_unit="MPa").sha256()
    assert a.sha256() == TimeHistory.from_arrays(time, S, stress_unit="MPa").sha256()


def test_to_block_matches_fe_contract_and_meta_is_local():
    time, S, E = _arrays(strain=True)
    h = TimeHistory.from_arrays(time, S, E, stress_unit="MPa")
    block = h.to_block(location_id=5)
    assert isinstance(block, HistoryBlock) and block.S_hist.shape == (1, 8, 3, 3)
    assert block.location_id[0] == 5 and np.array_equal(block.E_hist[0], E)
    meta = h.meta()
    assert meta.local_only and meta.strain_kind == "total_mechanical"
    assert TimeHistory.from_arrays(time, S, stress_unit="MPa").meta().strain_kind is None


def test_from_file_reads_one_location_history(tmp_path):
    pytest.importorskip("h5py")
    from multifat.fe.synthetic import (
        TubeGeometry,
        TubeLoad,
        TubeMaterial,
        write_tube_history,
    )

    loads = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 1000.0, 0.0), TubeLoad(50.0, 0.0, 8e4)]
    path = tmp_path / "one.h5"
    write_tube_history(
        path, TubeGeometry(6.0, 7.5), TubeMaterial(210000.0, 0.3), loads, n_radii=1, n_theta=1
    )
    h = TimeHistory.from_file(path)
    assert h.S_hist.shape == (3, 3, 3) and h.E_hist is not None
    assert np.allclose(h.time, [0.0, 0.5, 1.0]) and h.provenance["source"] == "history_file"
    path2 = tmp_path / "many.h5"
    write_tube_history(
        path2, TubeGeometry(6.0, 7.5), TubeMaterial(210000.0, 0.3), loads, n_radii=2, n_theta=3
    )
    with pytest.raises(ValueError, match="one location"):
        TimeHistory.from_file(path2)
    assert TimeHistory.from_file(path2, location_id=2).S_hist.shape == (3, 3, 3)


def test_axis_kind_and_time_unit_are_accepted_and_validated(tmp_path):
    time, S, _ = _arrays()
    h = TimeHistory.from_arrays(time, S, stress_unit="MPa", time_unit="s", axis_kind="time")
    assert h.units["time"] == "s"
    with pytest.raises(ValueError, match="time_unit"):
        TimeHistory.from_arrays(time, S, stress_unit="MPa", time_unit="")
    cycles = TimeHistory.from_arrays(np.arange(8.0), S, stress_unit="MPa", axis_kind="cycle")
    assert cycles.axis_kind == "cycle"
    cols = {"t": np.arange(8.0), "sxx": S[:, 0, 0], "syy": S[:, 1, 1], "szz": S[:, 2, 2],
            "sxy": S[:, 0, 1], "syz": S[:, 1, 2], "sxz": S[:, 0, 2]}  # fmt: skip
    path = tmp_path / "c.csv"
    pd.DataFrame(cols).to_csv(path, index=False)
    loaded = TimeHistory.from_csv(path, stress_unit="MPa", time_unit="cycles", axis_kind="cycle")
    assert loaded.axis_kind == "cycle"
    with pytest.raises(ValueError, match="time_unit"):
        TimeHistory.from_csv(path, stress_unit="MPa", time_unit="")


def test_record_detects_edited_metadata_and_json_rejects_nonfinite_provenance(tmp_path):
    import json

    time, S, _ = _arrays()
    h = TimeHistory.from_arrays(time, S, stress_unit="MPa", provenance={"note": "x"})
    path = h.record(tmp_path / "rec")
    info = json.loads(path.read_text(encoding="utf-8"))
    info["provenance"]["note"] = "edited"
    path.write_text(json.dumps(info), encoding="utf-8")
    with pytest.raises(ValueError, match="metadata"):
        TimeHistory.from_record(tmp_path / "rec")
    with pytest.raises(ValueError, match="JSON"):
        TimeHistory.from_arrays(time, S, stress_unit="MPa", provenance={"x": float("nan")})


def _node_csv_frame(time, S, E, shear_factor=1.0, *, strain_prefix="E_"):
    cols = {"time": time}
    pairs = zip(
        ("xx", "yy", "zz", "xy", "yz", "xz"), ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))
    )
    for name, (i, j) in pairs:
        cols["S_" + name] = S[:, i, j]
        if E is not None:
            cols[strain_prefix + name] = (shear_factor if i != j else 1.0) * E[:, i, j]
    return pd.DataFrame(cols)


def test_from_csv_accepts_the_node_csv_names(tmp_path):
    time, S, E = _arrays(strain=True)
    path = tmp_path / "h.csv"
    _node_csv_frame(time, S, E, 2.0).to_csv(path, index=False)
    h = TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s", strain_shear="engineering")
    assert np.allclose(h.time, time) and np.allclose(h.S_hist, S) and np.allclose(h.E_hist, E)
    _node_csv_frame(time, S, None).to_csv(path, index=False)
    assert TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s").E_hist is None


def test_from_csv_split_strain_columns_sum_to_the_total_and_eel_alone_is_refused(tmp_path):
    time, S, E = _arrays(strain=True)
    path = tmp_path / "h.csv"
    frame = _node_csv_frame(time, S, 0.75 * E, strain_prefix="Eel_")
    plastic = _node_csv_frame(time, S, 0.25 * E, strain_prefix="Epl_")
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="Eel"):
        TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s")
    both = pd.concat([frame, plastic.filter(like="Epl_")], axis=1)
    both.to_csv(path, index=False)
    h = TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s")
    assert np.allclose(h.E_hist, E)


def test_from_csv_refuses_both_names_of_one_quantity(tmp_path):
    time, S, E = _arrays(strain=True)
    path = tmp_path / "h.csv"
    frame = _node_csv_frame(time, S, E)
    frame["t"] = time
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="both"):
        TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s")
    frame = _node_csv_frame(time, S, E)
    frame["sxx"] = S[:, 0, 0]
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="both"):
        TimeHistory.from_csv(path, stress_unit="MPa", time_unit="s")


def _labelled_history_file(tmp_path, *, step_role="time_sample", count=2):
    pytest.importorskip("h5py")
    from multifat.fe.model import HistoryMeta, LoadCases
    from multifat.fe.store import HistoryWriter

    path = tmp_path / "labelled.h5"
    meta = HistoryMeta(
        material_id="",
        solver="test",
        solver_version="1",
        analysis_step="",
        adapter="test",
        adapter_version="1",
        source_files=(),
        source_units={"stress": "MPa", "length": "mm", "time": "s"},
        source_frame="global_cartesian",
        averaging="none",
        location_kind="node",
        analysis_type="linear_elastic",
        step_role=step_role,
    )
    with HistoryWriter(
        path,
        meta,
        np.ones(count, dtype=int),
        np.zeros((count, 3)),
        LoadCases(time=np.array([0.0, 1.0])),
        2,
        with_strain=False,
        instance=np.array(["A", "B"][:count]),
        material=np.array(["steel", "aluminium"][:count]),
    ) as writer:
        for i in range(2):
            writer.write_step(
                i, np.array([np.eye(3) * (10 * (j + 1) + i) for j in range(count)]), None
            )
    return path


def test_from_file_refuses_ambiguous_node_label(tmp_path):
    path = _labelled_history_file(tmp_path)
    with pytest.raises(ValueError, match="ambiguous|instance|material"):
        TimeHistory.from_file(path, location_id=1)


def test_from_file_refuses_unit_cases_even_with_time_axis(tmp_path):
    path = _labelled_history_file(tmp_path, step_role="unit_case", count=1)
    with pytest.raises(ValueError, match="unit_case|load.function"):
        TimeHistory.from_file(path)


def test_from_file_retains_selected_location_identity(tmp_path):
    path = _labelled_history_file(tmp_path, count=1)
    history = TimeHistory.from_file(path, location_id=1)
    assert history.provenance["instance"] == "A"
    assert history.provenance["material"] == "steel"
