"""HDF5 history store: streaming writer, reader, chunked iteration and bounded memory."""

from __future__ import annotations

import sys
import tracemalloc

import numpy as np
import pytest

from multifat.fe import store
from multifat.fe.model import HistoryMeta, LoadCases, SourceFile


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
        source_files=(SourceFile("a.vtu", "ab" * 32),),
        source_units={"stress": "Pa", "length": "m"},
        source_frame="cartesian",
        averaging="nodal_average",
        location_kind="node",
        notes="hello",
    )
    kw.update(override)
    return HistoryMeta(**kw)


def _write(path, p=5, n=4, strain=True, seed=0, dtype="float64", chunk_locations=4096, ids=None):
    rng = np.random.default_rng(seed)
    ids = np.arange(100, 100 + p) if ids is None else ids
    xyz = rng.normal(size=(p, 3))
    S = _sym(rng, (p, n))
    E = _sym(rng, (p, n)) * 1e-5
    lc = LoadCases(time=np.linspace(0.0, 1.0, n), labels=tuple(f"step{i}" for i in range(n)))
    with store.HistoryWriter(
        path, _meta(), ids, xyz, lc, n, with_strain=strain, dtype=dtype,
        chunk_locations=chunk_locations,
    ) as w:  # fmt: skip
        for k in range(n):
            w.write_step(k, S[:, k], E[:, k] if strain else None)
    return ids, xyz, S, E, lc


def test_without_h5py_message(monkeypatch):
    monkeypatch.setitem(sys.modules, "h5py", None)
    with pytest.raises(ImportError, match=r"multifat\[fe\]"):
        store.require_h5py()


def test_roundtrip_equals_input(tmp_path):
    pytest.importorskip("h5py")
    ids, xyz, S, E, lc = _write(tmp_path / "h.h5")
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert (f.n_locations, f.n_load_cases, f.has_strain) == (5, 4, True)
        assert np.array_equal(f.location_id, ids) and np.array_equal(f.coordinates, xyz)
        assert np.array_equal(f.load_cases.time, lc.time) and f.load_cases.labels == lc.labels
        block = f.read(slice(None))
    assert np.array_equal(block.S_hist, S) and np.array_equal(block.E_hist, E)


def test_float32_storage_within_tolerance(tmp_path):
    pytest.importorskip("h5py")
    _, _, S, _, _ = _write(tmp_path / "h.h5", dtype="float32")
    with store.HistoryFile(tmp_path / "h.h5") as f:
        block = f.read(slice(None))
    assert np.allclose(block.S_hist, S, rtol=1e-6, atol=1e-4)
    assert not np.array_equal(block.S_hist, S)


def test_missing_step_raises_on_close(tmp_path):
    pytest.importorskip("h5py")
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="load case"):
        with store.HistoryWriter(
            tmp_path / "h.h5", _meta(), np.arange(3), rng.normal(size=(3, 3)),
            LoadCases(), 2, with_strain=False,
        ) as w:  # fmt: skip
            w.write_step(0, _sym(rng, (3,)), None)


def test_provenance_and_sha_roundtrip(tmp_path):
    pytest.importorskip("h5py")
    _write(tmp_path / "h.h5")
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert f.meta == _meta()
        assert f.meta.source_files[0].sha256 == "ab" * 32


def test_duplicate_location_ids_rejected_by_writer(tmp_path):
    pytest.importorskip("h5py")
    with pytest.raises(ValueError, match="duplicate"):
        _write(tmp_path / "h.h5", p=3, ids=np.array([1, 1, 2]))


def test_non_integer_location_ids_rejected(tmp_path):
    pytest.importorskip("h5py")
    with pytest.raises(TypeError, match="integer"):
        _write(tmp_path / "h.h5", p=2, ids=np.array(["a", "b"]))


def test_write_step_rejects_asymmetric_and_nan(tmp_path):
    pytest.importorskip("h5py")
    rng = np.random.default_rng(0)
    w = store.HistoryWriter(
        tmp_path / "h.h5", _meta(), np.arange(2), rng.normal(size=(2, 3)), LoadCases(), 1,
        with_strain=False,
    )  # fmt: skip
    bad = _sym(rng, (2,))
    bad[0, 0, 1] += 5.0
    with pytest.raises(ValueError, match="symmetric"):
        w.write_step(0, bad, None)
    nan = _sym(rng, (2,))
    nan[1, 2, 2] = np.nan
    with pytest.raises(ValueError, match="finite"):
        w.write_step(0, nan, None)
    with pytest.raises(ValueError, match="step"):
        w.write_step(3, _sym(rng, (2,)), None)
    w.abort()


def test_iter_chunks_covers_all_locations_once(tmp_path):
    pytest.importorskip("h5py")
    ids, _, S, _, _ = _write(tmp_path / "h.h5", p=23)
    with store.HistoryFile(tmp_path / "h.h5") as f:
        blocks = list(f.iter_chunks(max_locations=5))
    assert [b.n_locations for b in blocks] == [5, 5, 5, 5, 3]
    assert np.array_equal(np.concatenate([b.location_id for b in blocks]), ids)
    assert np.array_equal(np.concatenate([b.S_hist for b in blocks]), S)


def test_iter_chunks_honours_byte_budget(tmp_path):
    pytest.importorskip("h5py")
    _write(tmp_path / "h.h5", p=50, n=4)
    with store.HistoryFile(tmp_path / "h.h5") as f:
        per_location = 4 * 6 * 8 * 2 * 3  # N * six components * float64 * (S and E) * expansion
        blocks = list(f.iter_chunks(max_bytes=per_location * 7))
    assert max(b.n_locations for b in blocks) == 7


def test_iter_chunks_over_selection_and_unsorted_read(tmp_path):
    pytest.importorskip("h5py")
    ids, _, S, _, _ = _write(tmp_path / "h.h5", p=10)
    with store.HistoryFile(tmp_path / "h.h5") as f:
        got = list(f.iter_chunks(max_locations=2, index=np.array([1, 4, 7])))
        block = f.read(np.array([7, 1]))
    assert [b.n_locations for b in got] == [2, 1]
    assert np.array_equal(block.location_id, ids[[7, 1]])
    assert np.array_equal(block.S_hist, S[[7, 1]])


def test_iter_chunks_empty_selection_yields_nothing(tmp_path):
    pytest.importorskip("h5py")
    _write(tmp_path / "h.h5")
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert list(f.iter_chunks(max_locations=2, index=np.array([], dtype=int))) == []


def test_single_load_case_file(tmp_path):
    pytest.importorskip("h5py")
    _, _, S, _, _ = _write(tmp_path / "h.h5", p=3, n=1)
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert np.array_equal(f.read(slice(None)).S_hist, S)


def _peak(path, max_locations):
    with store.HistoryFile(path) as f:
        tracemalloc.start()
        try:
            for block in f.iter_chunks(max_locations=max_locations):
                block.S_hist.sum()
            return tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()


def test_peak_memory_independent_of_P(tmp_path):
    pytest.importorskip("h5py")
    _write(tmp_path / "small.h5", p=2000, n=8, strain=False, chunk_locations=100)
    _write(tmp_path / "large.h5", p=8000, n=8, strain=False, chunk_locations=100)
    small, large = _peak(tmp_path / "small.h5", 100), _peak(tmp_path / "large.h5", 100)
    assert abs(large - small) < 0.25 * small


def _partial(path, raise_inside=False):
    rng = np.random.default_rng(0)
    with store.HistoryWriter(
        path, _meta(), np.arange(3), rng.normal(size=(3, 3)), LoadCases(), 2, with_strain=False
    ) as w:  # fmt: skip
        w.write_step(0, _sym(rng, (3,)), None)
        if raise_inside:
            raise RuntimeError("solver crashed")


def test_file_with_missing_load_cases_is_refused_on_open(tmp_path):
    pytest.importorskip("h5py")
    with pytest.raises(ValueError, match="load case"):
        _partial(tmp_path / "h.h5")
    with pytest.raises(ValueError, match="incomplete"):
        store.HistoryFile(tmp_path / "h.h5")


def test_file_after_exception_in_with_block_is_refused_on_open(tmp_path):
    pytest.importorskip("h5py")
    with pytest.raises(RuntimeError):
        _partial(tmp_path / "h.h5", raise_inside=True)
    with pytest.raises(ValueError, match="incomplete"):
        store.HistoryFile(tmp_path / "h.h5")


def test_complete_file_opens(tmp_path):
    pytest.importorskip("h5py")
    _write(tmp_path / "h.h5")
    with store.HistoryFile(tmp_path / "h.h5") as f:
        assert f.n_load_cases == 4


def _chunk_bytes(path, dataset="fields/stress"):
    import h5py

    with h5py.File(path, "r") as f:
        chunks, itemsize = f[dataset].chunks, f[dataset].dtype.itemsize
    return int(np.prod(chunks)) * itemsize


@pytest.mark.parametrize("n_load_cases,n_locations", [(8, 100), (3000, 64), (30_000, 4096)])
def test_default_chunk_stays_small_for_long_histories(tmp_path, n_load_cases, n_locations):
    pytest.importorskip("h5py")
    rng = np.random.default_rng(0)
    w = store.HistoryWriter(
        tmp_path / "h.h5", _meta(), np.arange(n_locations), rng.normal(size=(n_locations, 3)),
        LoadCases(), n_load_cases, with_strain=True,
    )  # fmt: skip
    w.abort()
    assert _chunk_bytes(tmp_path / "h.h5") <= 4 * 2**20
    assert _chunk_bytes(tmp_path / "h.h5", "fields/strain") <= 4 * 2**20


def _craft(path, ids, chunks):
    """A file in the format written by another tool: contiguous (unchunked) stress, no strain."""
    import json

    import h5py

    with h5py.File(path, "w") as f:
        f.attrs["format"], f.attrs["format_version"], f.attrs["complete"] = (
            "multifat-fe-history", "1.0", True,
        )  # fmt: skip
        f.create_group("meta").attrs["json"] = json.dumps(_meta().to_dict())
        loc = f.create_group("locations")
        loc.create_dataset("id", data=ids)
        loc.create_dataset("coordinates", data=np.zeros((len(ids), 3)))
        f.create_group("load_cases").create_dataset("index", data=np.arange(2))
        f.create_group("mesh")
        f.create_group("fields").create_dataset(
            "stress", data=np.zeros((len(ids), 2, 6)), chunks=chunks
        )


def test_unchunked_file_from_another_tool_opens_and_iterates(tmp_path):
    pytest.importorskip("h5py")
    _craft(tmp_path / "c.h5", np.arange(10), chunks=None)
    with store.HistoryFile(tmp_path / "c.h5") as f:
        assert [b.n_locations for b in f.iter_chunks(max_locations=4)] == [4, 4, 2]
        assert len(list(f.iter_chunks())) == 1


def test_duplicate_ids_rejected_when_opening(tmp_path):
    pytest.importorskip("h5py")
    _craft(tmp_path / "d.h5", np.array([1, 2, 2]), chunks=(3, 2, 6))
    with pytest.raises(ValueError, match="duplicate"):
        store.HistoryFile(tmp_path / "d.h5")
