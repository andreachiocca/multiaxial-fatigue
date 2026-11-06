"""Analytical thick-tube generator: checked against physical resultants, not its own formulas."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.fe.synthetic import (
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    tube_cartesian,
    tube_points,
    tube_strain_cyl,
    tube_stress_cyl,
    write_tube_history,
)

GEOM = TubeGeometry(a=6.0, b=7.5)
MAT = TubeMaterial(E=210000.0, nu=0.3)
LOAD = TubeLoad(p_i=50.0, F=12000.0, T=80000.0)


def _von_mises(T):
    dev = T - np.trace(T, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    return np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1)))


def test_boundary_conditions():
    s = tube_stress_cyl(np.array([GEOM.a, GEOM.b]), LOAD, GEOM)
    assert s[0, 0, 0] == pytest.approx(-LOAD.p_i, abs=1e-9)
    assert s[1, 0, 0] == pytest.approx(0.0, abs=1e-9)


def test_radial_equilibrium():
    r = np.linspace(GEOM.a, GEOM.b, 401)
    s = tube_stress_cyl(r, LOAD, GEOM)
    d_srr = np.gradient(s[:, 0, 0], r, edge_order=2)
    residual = d_srr + (s[:, 0, 0] - s[:, 1, 1]) / r
    assert np.max(np.abs(residual[2:-2])) < 1e-5 * LOAD.p_i


def _integral(values, r):
    return np.trapezoid(values * 2.0 * np.pi * r, r) if hasattr(np, "trapezoid") else np.trapz(
        values * 2.0 * np.pi * r, r
    )  # fmt: skip


def test_torque_resultant():
    r = np.linspace(GEOM.a, GEOM.b, 4001)
    tau = tube_stress_cyl(r, LOAD, GEOM)[:, 1, 2]
    assert _integral(tau * r, r) == pytest.approx(LOAD.T, rel=1e-6)


def test_axial_resultant_closed_and_open_ends():
    r = np.linspace(GEOM.a, GEOM.b, 4001)
    thrust = LOAD.p_i * np.pi * GEOM.a**2  # pressure acting on the closed end
    closed = tube_stress_cyl(r, LOAD, GEOM, closed_ends=True)[:, 2, 2]
    opened = tube_stress_cyl(r, LOAD, GEOM, closed_ends=False)[:, 2, 2]
    assert _integral(closed, r) == pytest.approx(LOAD.F + thrust, rel=1e-6)
    assert _integral(opened, r) == pytest.approx(LOAD.F, rel=1e-6)


def test_hooke_volumetric_relation():
    r = np.linspace(GEOM.a, GEOM.b, 9)
    s = tube_stress_cyl(r, LOAD, GEOM)
    e = tube_strain_cyl(s, MAT)
    tr_e = np.trace(e, axis1=-2, axis2=-1)
    tr_s = np.trace(s, axis1=-2, axis2=-1)
    assert np.allclose(tr_e, (1 - 2 * MAT.nu) / MAT.E * tr_s, rtol=1e-12)
    # tensor shear strain: eps_tz = tau (1 + nu) / E
    assert np.allclose(e[:, 1, 2], s[:, 1, 2] * (1 + MAT.nu) / MAT.E, rtol=1e-12)


def test_points_cover_radii_and_angles():
    ids, xyz, r, theta = tube_points(GEOM, n_radii=3, n_theta=8)
    assert len(ids) == 24 and len(np.unique(ids)) == 24
    assert np.allclose(np.hypot(xyz[:, 0], xyz[:, 1]), r)
    assert np.isclose(r.min(), GEOM.a) and np.isclose(r.max(), GEOM.b)


def test_cartesian_matches_cylindrical_at_known_angles_and_preserves_invariants():
    r = np.array([6.5, 6.5, 7.0])
    theta = np.array([0.0, np.pi / 2, 0.7])
    S, E = tube_cartesian(r, theta, LOAD, GEOM, MAT)
    cyl = tube_stress_cyl(r, LOAD, GEOM)
    # theta = 0: e_r = x, e_theta = y, e_z = z, so the tensors coincide
    assert np.allclose(S[0], cyl[0])
    # theta = pi/2: e_r = y, e_theta = -x
    assert np.isclose(S[1][1, 1], cyl[1][0, 0]) and np.isclose(S[1][0, 0], cyl[1][1, 1])
    assert np.isclose(S[1][0, 2], -cyl[1][1, 2])
    assert np.allclose(_von_mises(S), _von_mises(cyl))
    assert np.allclose(
        np.trace(E, axis1=-2, axis2=-1), np.trace(tube_strain_cyl(cyl, MAT), axis1=-2, axis2=-1)
    )


def test_write_tube_history_roundtrip(tmp_path):
    pytest.importorskip("h5py")
    from multifat.fe.store import HistoryFile

    loads = [
        TubeLoad(p, F, T) for p, F, T in [(0.0, 0.0, 0.0), (50.0, 12000.0, 0.0), (50.0, 0.0, 8e4)]
    ]
    path = tmp_path / "tube.h5"
    write_tube_history(path, GEOM, MAT, loads, n_radii=3, n_theta=6)
    ids, xyz, r, theta = tube_points(GEOM, n_radii=3, n_theta=6)
    with HistoryFile(path) as f:
        block = f.read(slice(None))
        assert f.meta.solver.startswith("analytical thick tube") and f.meta.local_only
    assert np.array_equal(block.location_id, ids) and np.allclose(block.coordinates, xyz)
    for k, load in enumerate(loads):
        S, E = tube_cartesian(r, theta, load, GEOM, MAT)
        assert np.array_equal(block.S_hist[:, k], S) and np.array_equal(block.E_hist[:, k], E)


def test_documented_example_runs(tmp_path):
    pytest.importorskip("h5py")
    from multifat.fe.store import HistoryFile
    from multifat.fe.synthetic import main
    from multifat.fe.tensors import CANONICAL
    from multifat.fe.vtk import TensorSpec, import_pvd

    assert main(["--out", str(tmp_path)]) == 0
    direct, pvd = tmp_path / "tube_history.h5", tmp_path / "vtk" / "tube.pvd"
    assert direct.exists() and pvd.exists()
    import_pvd(
        pvd, tmp_path / "from_vtk.h5", stress=TensorSpec("Stress", CANONICAL, unit="MPa"),
        strain=TensorSpec("Strain", CANONICAL), id_array="NodeID", material_id="synthetic_tube",
        solver="analytical thick tube", solver_version="1", analysis_step="synthetic", length_unit="mm",
    )  # fmt: skip
    with HistoryFile(direct) as a, HistoryFile(tmp_path / "from_vtk.h5") as b:
        x, y = a.read(slice(None)), b.read(slice(None))
    assert np.array_equal(x.location_id, y.location_id)
    assert np.allclose(x.S_hist, y.S_hist, rtol=1e-10, atol=1e-9)
    assert np.allclose(x.E_hist, y.E_hist, rtol=1e-10, atol=1e-15)
