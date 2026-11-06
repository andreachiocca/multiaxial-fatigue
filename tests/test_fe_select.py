"""Selection helpers and screening of FE locations."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("h5py")

from multifat.fe import store  # noqa: E402
from multifat.fe.model import HistoryMeta, LoadCases  # noqa: E402
from multifat.fe.select import (  # noqa: E402
    screen_by,
    select_by_ids,
    select_by_mask,
    select_in_box,
    von_mises_max,
)


def _file(tmp_path, p=40, n=5):
    rng = np.random.default_rng(3)
    ids = np.arange(1000, 1000 + p)[::-1].copy()  # IDs are not in file order
    xyz = rng.uniform(0.0, 10.0, size=(p, 3))
    a = rng.normal(size=(p, n, 3, 3)) * 50.0
    S = 0.5 * (a + np.swapaxes(a, -1, -2))
    meta = HistoryMeta(
        material_id="m", solver="s", solver_version="1", analysis_step="a", adapter="t",
        adapter_version="0", source_files=(), source_units={}, source_frame="cartesian",
        averaging="none", location_kind="node",
    )  # fmt: skip
    path = tmp_path / "h.h5"
    with store.HistoryWriter(path, meta, ids, xyz, LoadCases(), n, with_strain=False) as w:
        for k in range(n):
            w.write_step(k, S[:, k], None)
    return path, ids, xyz, S


def test_select_by_ids_returns_sorted_positions_and_unknown_raises(tmp_path):
    path, ids, _, _ = _file(tmp_path)
    with store.HistoryFile(path) as f:
        pos = select_by_ids(f, [ids[7], ids[2], ids[30]])
        assert list(pos) == [2, 7, 30]
        with pytest.raises(KeyError, match="999"):
            select_by_ids(f, [ids[0], 999])


def test_select_by_mask(tmp_path):
    path, _, xyz, _ = _file(tmp_path)
    with store.HistoryFile(path) as f:
        mask = xyz[:, 0] > 5.0
        assert np.array_equal(select_by_mask(f, mask), np.flatnonzero(mask))
        with pytest.raises(ValueError, match="mask"):
            select_by_mask(f, mask[:-1])


def test_select_in_box(tmp_path):
    path, _, xyz, _ = _file(tmp_path)
    lo, hi = np.array([2.0, 0.0, 0.0]), np.array([6.0, 10.0, 5.0])
    expected = np.flatnonzero(np.all((xyz >= lo) & (xyz <= hi), axis=1))
    with store.HistoryFile(path) as f:
        assert np.array_equal(select_in_box(f, lo, hi), expected)


def test_screen_by_top_k_matches_bruteforce(tmp_path):
    path, _, _, S = _file(tmp_path)
    dev = S - np.trace(S, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    vm = np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1))).max(axis=1)
    with store.HistoryFile(path) as f:
        top = screen_by(f, von_mises_max, top=6, max_locations=7)
        thr = screen_by(f, von_mises_max, threshold=float(np.sort(vm)[-4]), max_locations=7)
        with pytest.raises(ValueError, match="exactly one"):
            screen_by(f, von_mises_max)
    assert set(top) == set(np.argsort(vm)[-6:]) and list(top) == sorted(top)
    assert set(thr) == set(np.flatnonzero(vm >= np.sort(vm)[-4]))


def test_screen_by_rejects_nonpositive_top(tmp_path):
    path, _, _, _ = _file(tmp_path)
    with store.HistoryFile(path) as f:
        for top in (0, -3):
            with pytest.raises(ValueError, match="top"):
                screen_by(f, von_mises_max, top=top)


def test_screen_by_ignores_nan_scalars(tmp_path):
    path, ids, _, _ = _file(tmp_path)

    def scalar(block):
        values = von_mises_max(block)
        values[block.location_id % 3 == 0] = np.nan
        return values

    with store.HistoryFile(path) as f:
        top = screen_by(f, scalar, top=8, max_locations=7)
        thr = screen_by(f, scalar, threshold=0.0, max_locations=7)
    assert len(top) == 8 and not np.any(ids[top] % 3 == 0)
    assert not np.any(ids[thr] % 3 == 0) and len(thr) == int(np.sum(ids % 3 != 0))
