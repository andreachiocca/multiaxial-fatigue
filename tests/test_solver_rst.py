"""Optional native RST contract tests and opt-in public-result verification."""

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from multifat.fe.solvers import rst
from multifat.fe.csvio import import_fe_csv
from multifat.fe.store import HistoryFile

OPTIONS = dict(
    stress_unit="MPa",
    length_unit="mm",
    analysis_type="linear_elastic",
    step_role="time_sample",
    coordinate_system="global_cartesian",
    time_unit="s",
)


def fake_result():
    mesh = SimpleNamespace(
        nnum=np.array([1, 2]),
        nodes=np.array([[0, 0, 0], [1, 2, 3]]),
        material_type=np.array([1]),
        ekey=np.array([[1, 186]]),
        key_option={},
    )
    return SimpleNamespace(
        mesh=mesh,
        nsets=2,
        time_values=np.array([0.0, 1.0]),
        _is_cyclic=False,
        _resultheader=dict(
            kan=0,
            cpxrst=0,
            nlgeom=0,
            units="User Defined",
            verstring="test",
            ls_table=np.array([[1, 1, 1], [1, 2, 2]]),
        ),
        nodal_stress=lambda i: (
            np.array([2, 1]),
            np.array([[41, 52, 63, 42, 53, 43], [11, 22, 33, 12, 23, 13]]) * (i + 1),
        ),
    )


def test_native_component_order_and_node_join(tmp_path, monkeypatch):
    path = tmp_path / "input.rst"
    path.write_bytes(b"contract placeholder, not solver data")
    monkeypatch.setattr(rst, "_read_binary", lambda path: fake_result())
    header = rst.ansys_rst_to_csv(path, tmp_path / "out.csv", **OPTIONS)
    assert header.averaging == "nodal_average_global"
    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        np.testing.assert_array_equal(
            block.S_hist[0, 0], [[41, 42, 43], [42, 52, 53], [43, 53, 63]]
        )
        np.testing.assert_array_equal(block.coordinates[0], [1, 2, 3])


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda r: r._resultheader.update(kan=2), "static|transient"),
        (lambda r: r._resultheader.update(cpxrst=1), "complex"),
        (lambda r: setattr(r, "_is_cyclic", True), "cyclic"),
        (lambda r: setattr(r.mesh, "material_type", np.array([1, 2])), "material"),
        (lambda r: setattr(r.mesh, "ekey", np.array([[1, 181]])), "solid"),
        (lambda r: setattr(r, "time_values", np.array([1, 0])), "time"),
        (lambda r: r._resultheader.update(nlgeom=1), "geometr"),
        (lambda r: r._resultheader.update(units="SI"), "unit"),
        (
            lambda r: setattr(
                r, "nodal_stress", lambda i: (np.array([1]), np.full((1, 6), np.nan))
            ),
            "finite",
        ),
    ],
)
def test_native_refusals(tmp_path, monkeypatch, change, match):
    path = tmp_path / "input.rst"
    path.write_bytes(b"placeholder")
    result = fake_result()
    change(result)
    monkeypatch.setattr(rst, "_read_binary", lambda path: result)
    with pytest.raises(ValueError, match=match):
        rst.ansys_rst_to_csv(path, tmp_path / "out.csv", **OPTIONS)
    assert not (tmp_path / "out.csv").exists()


def test_inhouse_rst_fatigue_cycle(tmp_path):
    reader = pytest.importorskip("ansys.mapdl.reader")
    path = Path(__file__).parent / "data/solvers/ansys/hole_plate_fatigue.rst"
    if not path.exists():
        pytest.skip("in-house holed plate fatigue fixture not found")

    result = reader.read_binary(path)
    assert result.nsets == 2
    ids, stress_0 = result.nodal_stress(0)
    finite = np.isfinite(stress_0).all(axis=1)
    selected = ids[finite]
    assert len(selected) == 2864

    header = rst.ansys_rst_to_csv(path, tmp_path / "out.csv", nodes=selected, **OPTIONS)
    assert header.n_nodes == 2864
    assert header.n_steps == 2

    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        assert block.n_locations == 2864
        assert block.S_hist.shape[1] == 2

        coords = block.coordinates
        # Find Notch A at (10, 0, 0) and Notch B at (0, 10, 0)
        dist_A = np.hypot(coords[:, 0] - 10.0, coords[:, 1] - 0.0)
        dist_B = np.hypot(coords[:, 0] - 0.0, coords[:, 1] - 10.0)
        idx_A = np.argmin(dist_A)
        idx_B = np.argmin(dist_B)

        s1 = block.S_hist[:, 0]
        s2 = block.S_hist[:, 1]
        # Step 1: Notch A is in tension (SY ~ 304), Notch B in compression (SX ~ -79)
        assert s1[idx_A, 1, 1] > 295.0
        assert s1[idx_B, 0, 0] < -70.0

        # Step 2: Notch B is in tension (SX ~ 304), Notch A in compression (SY ~ -79)
        assert s2[idx_B, 0, 0] > 295.0
        assert s2[idx_A, 1, 1] < -70.0

        # Peak fatigue stress range over the 2-step cycle:
        delta_s_A = s1[idx_A, 1, 1] - s2[idx_A, 1, 1]
        delta_s_B = s2[idx_B, 0, 0] - s1[idx_B, 0, 0]
        assert delta_s_A > 360.0
        assert delta_s_B > 360.0


@pytest.mark.skipif(
    not os.environ.get("MULTIFAT_TEST_RST"),
    reason="set MULTIFAT_TEST_RST to public vol_test.rst; no test downloads",
)
def test_real_public_rst(tmp_path):
    reader = pytest.importorskip("ansys.mapdl.reader")
    path = Path(os.environ["MULTIFAT_TEST_RST"])
    result = reader.read_binary(path)
    ids, stress = result.nodal_stress(0)
    finite = np.isfinite(stress).all(axis=1)
    # Quadratic midside nodes without stress are selected out explicitly, never silently.
    selected = ids[finite]
    header = rst.ansys_rst_to_csv(path, tmp_path / "out.csv", nodes=selected, **OPTIONS)
    assert header.n_nodes == int(finite.sum())
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        expected = stress[finite]
        np.testing.assert_allclose(block.S_hist[:, 0, 0, 2], expected[:, 5])
        np.testing.assert_allclose(block.S_hist[:, 0, 1, 2], expected[:, 4])
        np.testing.assert_allclose(block.S_hist[:, 0, 0, 0], expected[:, 0])
    # Independently assemble global element tensors at one node using the reader's
    # explicitly global element API, rather than only repeating nodal_stress.
    _, tensors, element_nodes = result.element_stress(0, in_element_coord_sys=False)
    node = int(selected[0])
    contributions = [s[np.asarray(n) == node] for s, n in zip(tensors, element_nodes) if node in n]
    np.testing.assert_allclose(
        expected[0], np.concatenate(contributions).mean(axis=0), rtol=1e-10, atol=1e-12
    )
