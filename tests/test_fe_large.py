"""Large FE inputs are streamed with memory set by the chunk budget, not by the location count."""

from __future__ import annotations

import os
import tracemalloc

import numpy as np
import pytest

pytest.importorskip("h5py")

from multifat.fe.store import HistoryFile  # noqa: E402
from multifat.fe.synthetic import (  # noqa: E402
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    write_tube_history,
)

GEOM, MAT = TubeGeometry(6.0, 7.5), TubeMaterial(210000.0, 0.3)


def _loads(n):
    return [TubeLoad(50.0 * np.sin(np.pi * k / n), 1000.0 * k, 5000.0 * (n - k)) for k in range(n)]


def _stream(path, n_locations, budget, dtype):
    seen = 0
    ids_sum = 0
    with HistoryFile(path) as f:
        assert f.n_locations == n_locations
        tracemalloc.start()
        try:
            for block in f.iter_chunks(max_bytes=budget):
                seen += block.n_locations
                ids_sum += int(block.location_id.sum())
                block.S_hist.sum()
            peak = tracemalloc.get_traced_memory()[1]
        finally:
            tracemalloc.stop()
    return seen, ids_sum, peak


def test_stream_many_locations_bounded(tmp_path):
    p, n, budget = 40_000, 8, 4_000_000
    write_tube_history(tmp_path / "big.h5", GEOM, MAT, _loads(n), n_radii=200, n_theta=200)
    seen, ids_sum, peak = _stream(tmp_path / "big.h5", p, budget, "float64")
    assert seen == p and ids_sum == p * (p + 1) // 2  # every location exactly once
    assert peak < 3 * budget


@pytest.mark.skipif(
    os.environ.get("MULTIFAT_SLOW") != "1", reason="set MULTIFAT_SLOW=1 (writes about 1 GB)"
)
def test_large_scale_opt_in(tmp_path):
    p, n, budget = 1_000_000, 16, 64_000_000
    write_tube_history(
        tmp_path / "huge.h5", GEOM, MAT, _loads(n), n_radii=1000, n_theta=1000, dtype="float32"
    )
    seen, ids_sum, peak = _stream(tmp_path / "huge.h5", p, budget, "float32")
    assert seen == p and ids_sum == p * (p + 1) // 2
    assert peak < 4 * budget
