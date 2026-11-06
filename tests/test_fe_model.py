"""Data model of the FE interface: validation of blocks, metadata round trip, file hashing."""

from __future__ import annotations

import hashlib

import numpy as np
import pytest

from multifat.fe.model import HistoryBlock, HistoryMeta, LoadCases, SourceFile, sha256_of


def _sym(rng, shape):
    a = rng.normal(size=shape + (3, 3))
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def _block(m=4, n=3, strain=True, **override):
    rng = np.random.default_rng(0)
    kw = dict(
        location_id=np.arange(m),
        coordinates=rng.normal(size=(m, 3)),
        S_hist=_sym(rng, (m, n)),
        E_hist=_sym(rng, (m, n)) * 1e-3 if strain else None,
    )
    kw.update(override)
    return HistoryBlock(**kw)


def _meta(**override):
    kw = dict(
        material_id="42CrMo4_QT",
        solver="synthetic",
        solver_version="1",
        analysis_step="step-1",
        adapter="multifat.fe.synthetic",
        adapter_version="1.0",
        source_files=(SourceFile("a.vtu", "0" * 64),),
        source_units={"stress": "Pa", "length": "m"},
        source_frame="cylindrical",
        averaging="none",
        location_kind="node",
    )
    kw.update(override)
    return HistoryMeta(**kw)


def test_block_accepts_valid():
    b = _block()
    assert b.S_hist.shape == (4, 3, 3, 3) and b.E_hist.shape == (4, 3, 3, 3)
    assert _block(strain=False).E_hist is None


def test_block_rejects_nan():
    S = _block().S_hist.copy()
    S[1, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        _block(S_hist=S)


def test_block_rejects_asymmetric():
    S = _block().S_hist.copy()
    S[0, 0, 0, 1] += 1.0
    with pytest.raises(ValueError, match="symmetric"):
        _block(S_hist=S)


def test_block_rejects_duplicate_ids():
    with pytest.raises(ValueError, match="duplicate"):
        _block(location_id=np.array([1, 2, 2, 3]))


def test_block_rejects_shape_mismatch_between_S_and_E():
    E = _block().E_hist[:, :2]
    with pytest.raises(ValueError, match="shape"):
        _block(E_hist=E)


def test_block_rejects_empty_history():
    with pytest.raises(ValueError, match="load case"):
        _block(S_hist=np.zeros((4, 0, 3, 3)), E_hist=None)


def test_meta_dict_roundtrip():
    meta = _meta(notes="x", periodic=True)
    assert HistoryMeta.from_dict(meta.to_dict()) == meta


def test_meta_rejects_unknown_averaging():
    with pytest.raises(ValueError, match="averaging"):
        _meta(averaging="smoothed")


def test_load_cases_accepts_optional_axes():
    assert LoadCases(time=None, labels=None).time is None
    assert LoadCases(time=np.array([0.0, 1.0]), labels=("a", "b")).labels == ("a", "b")


def test_sha256_of_known_bytes(tmp_path):
    p = tmp_path / "f.bin"
    p.write_bytes(b"abc")
    assert sha256_of(p) == hashlib.sha256(b"abc").hexdigest()


def test_meta_from_dict_ignores_fields_of_newer_minor_versions():
    data = _meta().to_dict()
    data["future_field"] = "x"
    assert HistoryMeta.from_dict(data) == _meta()


def test_meta_conventions_roundtrip_and_default_for_older_files():
    meta = _meta(source_conventions={"stress_order": "xx,yy,zz,xy,yz,xz", "strain_shear": "tensor"})
    assert HistoryMeta.from_dict(meta.to_dict()) == meta
    data = _meta().to_dict()
    data.pop("source_conventions")
    assert HistoryMeta.from_dict(data).source_conventions == {}
