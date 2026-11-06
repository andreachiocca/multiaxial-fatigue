"""History file format 1.1: instance, material, split strain, temperature, steps and unit tables."""

from __future__ import annotations

import json

import numpy as np
import pytest

from multifat.fe import store
from multifat.fe.model import HistoryBlock, HistoryMeta, LoadCases, SourceFile
from multifat.fe.select import select_by_ids
from multifat.fe.tensors import elastic_strain, to_mm, to_mpa, to_strain

h5py = pytest.importorskip("h5py")


def _sym(rng, shape):
    a = rng.normal(size=shape + (3, 3)) * 100.0
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def _meta(**override):
    kw = dict(
        material_id="mat",
        solver="synthetic",
        solver_version="1",
        analysis_step="s1",
        adapter="test",
        adapter_version="0",
        source_files=(SourceFile("a.csv", "ab" * 32),),
        source_units={"stress": "MPa", "length": "mm"},
        source_frame="global_cartesian",
        averaging="none",
        location_kind="node",
    )
    kw.update(override)
    return HistoryMeta(**kw)


def test_version_1_1_roundtrip_of_instance_material_split_strain_temperature_and_steps(tmp_path):
    rng = np.random.default_rng(1)
    p, n = 4, 3
    ids = np.array([1, 2, 3, 4])
    inst = np.array(["A", "A", "B", "B"])
    mat = np.array(["steel", "alu", "steel", "alu"])
    xyz = rng.normal(size=(p, 3))
    S = _sym(rng, (p, n))
    Eel = _sym(rng, (p, n)) * 1e-5
    Epl = _sym(rng, (p, n)) * 1e-6
    T = rng.uniform(20.0, 80.0, size=(p, n))
    lc = LoadCases(time=np.array([0.0, 0.5, 1.0]), steps=np.array([3, 7, 12]))
    meta = _meta(
        analysis_type="elastic_plastic",
        cycle="stabilised",
        step_role="time_sample",
        load_description="max-min",
        strain_kind="total_mechanical",
        averaging="nodal_average_within_material",
    )
    with store.HistoryWriter(
        tmp_path / "h.h5", meta, ids, xyz, lc, n, with_strain=True,
        instance=inst, material=mat, with_strain_split=True, with_temperature=True,
    ) as w:  # fmt: skip
        for k in range(n):
            w.write_step(k, S[:, k], Eel[:, k] + Epl[:, k],
                         E_elastic=Eel[:, k], E_plastic=Epl[:, k], T=T[:, k])  # fmt: skip
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert store.FORMAT_VERSION == "1.1"
        assert f.meta.analysis_type == "elastic_plastic"
        assert f.meta.cycle == "stabilised" and f.meta.load_description == "max-min"
        assert f.meta.averaging == "nodal_average_within_material"
        assert list(f.instance) == ["A", "A", "B", "B"]
        assert list(f.material) == ["steel", "alu", "steel", "alu"]
        assert f.has_strain_split and f.has_temperature
        assert list(f.load_cases.steps) == [3, 7, 12]
        b = f.read(slice(None))
        assert np.allclose(b.S_hist, S)
        assert np.allclose(b.E_elastic_hist, Eel) and np.allclose(b.E_plastic_hist, Epl)
        assert np.allclose(b.E_hist, Eel + Epl)
        assert np.allclose(b.T_hist, T)
        assert list(b.instance) == ["A", "A", "B", "B"]
        assert list(b.material) == ["steel", "alu", "steel", "alu"]


def _write_locations(path, ids, inst, mat, p=None):
    p = len(ids)
    rng = np.random.default_rng(0)
    S = _sym(rng, (p, 2))
    with store.HistoryWriter(
        path, _meta(), np.asarray(ids), np.zeros((p, 3)), LoadCases(), 2, with_strain=False,
        instance=inst, material=mat,
    ) as w:  # fmt: skip
        for k in range(2):
            w.write_step(k, S[:, k], None)


def test_same_node_label_in_two_instances_and_an_interface_node_in_two_materials_are_distinct_locations(
    tmp_path,
):
    _write_locations(
        tmp_path / "h.h5", [5, 5, 9, 9], ["A", "B", "A", "A"], ["m1", "m1", "m1", "m2"]
    )
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert f.n_locations == 4


def test_duplicate_location_triple_rejected(tmp_path):
    with pytest.raises(ValueError, match="duplicate"):
        _write_locations(tmp_path / "h.h5", [5, 5], ["A", "A"], ["m1", "m1"])


def test_node_id_only_duplicates_still_rejected_without_instance_or_material(tmp_path):
    with pytest.raises(ValueError, match="duplicate"):
        _write_locations(tmp_path / "h.h5", [5, 5], None, None)


def test_block_uniqueness_follows_the_triple():
    kw = dict(coordinates=np.zeros((2, 3)), S_hist=np.zeros((2, 1, 3, 3)))
    HistoryBlock(np.array([5, 5]), instance=np.array(["A", "B"]), **kw)
    with pytest.raises(ValueError, match="duplicate"):
        HistoryBlock(np.array([5, 5]), **kw)
    with pytest.raises(ValueError, match="duplicate"):
        HistoryBlock(np.array([5, 5]), instance=np.array(["A", "A"]), **kw)


def test_version_1_0_file_still_opens(tmp_path):
    with h5py.File(tmp_path / "old.h5", "w") as f:
        f.attrs["format"], f.attrs["format_version"], f.attrs["complete"] = (
            "multifat-fe-history", "1.0", True,
        )  # fmt: skip
        old = _meta().to_dict()
        for key in ("analysis_type", "cycle", "step_role", "load_description"):
            old.pop(key, None)
        f.create_group("meta").attrs["json"] = json.dumps(old)
        loc = f.create_group("locations")
        loc.create_dataset("id", data=np.arange(3))
        loc.create_dataset("coordinates", data=np.zeros((3, 3)))
        f.create_group("load_cases").create_dataset("index", data=np.arange(2))
        f.create_group("mesh")
        f.create_group("fields").create_dataset("stress", data=np.zeros((3, 2, 6)))
    with store.HistoryFile(tmp_path / "old.h5") as f:
        assert f.instance is None and f.material is None
        assert not f.has_strain_split and not f.has_temperature
        assert f.load_cases.steps is None
        assert f.meta.analysis_type == "unspecified" and f.meta.step_role == "time_sample"
        assert f.read(slice(None)).n_locations == 3


def test_select_by_ids_with_instance_and_material_and_ambiguous_without(tmp_path):
    _write_locations(
        tmp_path / "h.h5", [5, 5, 9, 9], ["A", "B", "A", "A"], ["m1", "m1", "m1", "m2"]
    )
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert list(select_by_ids(f, [5], instance="B")) == [1]
        assert list(select_by_ids(f, [9], material="m2")) == [3]
        assert list(select_by_ids(f, [9], instance="A", material="m1")) == [2]
        assert list(select_by_ids(f, [5])) == [0, 1]  # omitted filters return every match
        assert list(select_by_ids(f, [9])) == [2, 3]
        with pytest.raises(KeyError):
            select_by_ids(f, [5], instance="C")
        with pytest.raises(KeyError):
            select_by_ids(f, [77])


def test_meta_new_fields_roundtrip_and_defaults_for_older_dicts():
    m = _meta(analysis_type="linear_elastic", step_role="unit_case", cycle="3")
    again = HistoryMeta.from_dict(m.to_dict())
    assert again == m
    older = _meta().to_dict()
    for key in ("analysis_type", "cycle", "step_role", "load_description"):
        older.pop(key)
    got = HistoryMeta.from_dict(older)
    assert (got.analysis_type, got.cycle, got.step_role, got.load_description) == (
        "unspecified", "unspecified", "time_sample", "",
    )  # fmt: skip
    with pytest.raises(ValueError, match="analysis_type"):
        _meta(analysis_type="plastic")
    with pytest.raises(ValueError, match="step_role"):
        _meta(step_role="sample")
    with pytest.raises(ValueError, match="strain_kind"):
        _meta(strain_kind="weird")
    _meta(strain_kind="derived_linear_elastic")
    _meta(strain_kind="derived_elastic_part")


def test_averaging_values_validated():
    for ok in ("none", "nodal_average", "extrapolated", "nodal_average_within_material",
               "nodal_average_global"):  # fmt: skip
        _meta(averaging=ok)
    with pytest.raises(ValueError, match="averaging"):
        _meta(averaging="smooth")


def test_load_cases_steps_validated():
    LoadCases(steps=np.array([1, 4, 9]))
    with pytest.raises(ValueError, match="steps"):
        LoadCases(steps=np.array([[1, 2]]))
    with pytest.raises(ValueError, match="differ"):
        LoadCases(time=np.array([0.0, 1.0]), steps=np.array([1, 2, 3]))


def test_writer_refuses_missing_split_or_temperature_arguments(tmp_path):
    ids = np.arange(2)
    w = store.HistoryWriter(
        tmp_path / "h.h5", _meta(), ids, np.zeros((2, 3)), LoadCases(), 1, with_strain=True,
        with_strain_split=True, with_temperature=True,
    )  # fmt: skip
    z = np.zeros((2, 3, 3))
    with pytest.raises(ValueError, match="split"):
        w.write_step(0, z, z, T=np.zeros(2))
    with pytest.raises(ValueError, match="temperature"):
        w.write_step(0, z, z, E_elastic=z, E_plastic=z)
    w.abort()


def test_elastic_strain_matches_hooke_by_hand():
    E, nu = 200000.0, 0.3
    S = np.zeros((1, 3, 3))
    S[0, 0, 0] = 100.0  # uniaxial: eps_xx = s/E, eps_yy = eps_zz = -nu s/E, no shear
    e = elastic_strain(S, E, nu)
    assert e[0, 0, 0] == pytest.approx(100.0 / E)
    assert e[0, 1, 1] == pytest.approx(-nu * 100.0 / E)
    assert e[0, 2, 2] == pytest.approx(-nu * 100.0 / E)
    assert e[0, 0, 1] == 0.0
    S2 = np.zeros((1, 3, 3))
    S2[0, 0, 1] = S2[0, 1, 0] = 80.0  # shear: tensor shear eps_xy = (1+nu) tau / E
    e2 = elastic_strain(S2, E, nu)
    assert e2[0, 0, 1] == pytest.approx((1 + nu) * 80.0 / E)
    assert e2[0, 1, 0] == e2[0, 0, 1]
    with pytest.raises(ValueError):
        elastic_strain(S, -1.0, nu)
    with pytest.raises(ValueError):
        elastic_strain(S, E, 0.5)


def test_us_customary_and_strain_units_convert():
    assert to_mpa(np.array([1.0]), "psi")[0] == pytest.approx(6.894757293168361e-3)
    assert to_mpa(np.array([1.0]), "ksi")[0] == pytest.approx(6.894757293168361)
    assert to_mpa(np.array([5.0]), "N/mm^2")[0] == 5.0
    assert to_mm(np.array([2.0]), "in")[0] == pytest.approx(50.8)
    assert to_mm(np.array([1.0]), "ft")[0] == pytest.approx(304.8)
    assert to_strain(np.array([3.0]), "percent")[0] == pytest.approx(0.03)
    assert to_strain(np.array([500.0]), "microstrain")[0] == pytest.approx(5e-4)
    for unit in ("1", "mm/mm", "m/m"):
        assert to_strain(np.array([0.002]), unit)[0] == 0.002
    for fn in (to_mpa, to_mm, to_strain):
        with pytest.raises(ValueError, match="unit"):
            fn(np.array([1.0]), "furlong")
