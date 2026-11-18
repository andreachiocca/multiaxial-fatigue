"""Tensor conventions of the FE interface: component order, shear convention, units, frames."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.fe.tensors import (
    CANONICAL,
    VTK9,
    components_to_sym3,
    cylindrical_rotations,
    rotate_to_global,
    six_to_sym3,
    sym3_to_six,
    to_mm,
    to_mpa,
)


def _random_sym(rng, shape=()):
    a = rng.normal(size=shape + (3, 3))
    return 0.5 * (a + np.swapaxes(a, -1, -2))


def von_mises(T):
    dev = T - np.trace(T, axis1=-2, axis2=-1)[..., None, None] / 3.0 * np.eye(3)
    return np.sqrt(1.5 * np.sum(dev * dev, axis=(-2, -1)))


def test_six_roundtrip():
    T = _random_sym(np.random.default_rng(0), (7,))
    assert np.array_equal(six_to_sym3(sym3_to_six(T)), T)


def test_component_order_detects_permutation():
    row = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]])
    T = components_to_sym3(row, CANONICAL, kind="stress")[0]
    assert (T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 2], T[0, 2]) == (1, 2, 3, 4, 5, 6)
    swapped = components_to_sym3(row, ("xx", "yy", "zz", "xz", "yz", "xy"), kind="stress")[0]
    assert swapped[0, 1] == 6 and swapped[0, 2] == 4


def test_engineering_shear_halved_for_strain_only():
    row = np.array([[0.0, 0.0, 0.0, 4.0, 6.0, 8.0]])
    strain = components_to_sym3(row, CANONICAL, kind="strain", shear="engineering")[0]
    stress = components_to_sym3(row, CANONICAL, kind="stress", shear="engineering")[0]
    assert (strain[0, 1], strain[1, 2], strain[0, 2]) == (2.0, 3.0, 4.0)
    assert (stress[0, 1], stress[1, 2], stress[0, 2]) == (4.0, 6.0, 8.0)


def test_vtk9_roundtrip_and_asymmetric_pair_rejected():
    T = _random_sym(np.random.default_rng(1), (3,))
    nine = T.reshape(3, 9)
    assert np.allclose(components_to_sym3(nine, VTK9, kind="stress"), T)
    bad = nine.copy()
    bad[1, 1] += 1.0
    with pytest.raises(ValueError, match="symmetric"):
        components_to_sym3(bad, VTK9, kind="stress")


def test_unit_scaling():
    assert to_mpa(np.array([2.0e6]), "Pa")[0] == pytest.approx(2.0)
    assert to_mpa(np.array([3.0]), "GPa")[0] == pytest.approx(3000.0)
    assert to_mm(np.array([0.5]), "m")[0] == pytest.approx(500.0)
    with pytest.raises(ValueError, match="unit"):
        to_mpa(np.array([1.0]), "furlong")


def test_cylindrical_to_cartesian_analytic():
    pts = np.array([[2.0, 0.0, 5.0], [0.0, 3.0, -1.0]])
    R = cylindrical_rotations(pts)
    cyl = np.tile(np.diag([10.0, 20.0, 30.0]), (2, 1, 1))
    cart = rotate_to_global(cyl, R)
    # at (r, 0): e_r = x, e_theta = y; at (0, r): e_r = y, e_theta = -x
    assert np.allclose(np.diag(cart[0]), [10.0, 20.0, 30.0])
    assert np.allclose(np.diag(cart[1]), [20.0, 10.0, 30.0])


def test_cylindrical_rotation_preserves_invariants_at_random_points():
    rng = np.random.default_rng(2)
    pts = rng.normal(size=(30, 3)) + np.array([1.5, 0.5, 0.0])
    axis = np.array([0.2, -0.3, 1.0])
    R = cylindrical_rotations(pts, axis=axis, origin=(0.1, 0.2, 0.3))
    assert np.allclose(R @ np.swapaxes(R, -1, -2), np.eye(3), atol=1e-12)
    assert np.allclose(np.linalg.det(R), 1.0)
    cyl = _random_sym(rng, (30,))
    cart = rotate_to_global(cyl, R)
    assert np.allclose(np.trace(cart, axis1=-2, axis2=-1), np.trace(cyl, axis1=-2, axis2=-1))
    assert np.allclose(von_mises(cart), von_mises(cyl))


def test_point_on_axis_rejected():
    with pytest.raises(ValueError, match="axis"):
        cylindrical_rotations(np.array([[0.0, 0.0, 1.0]]))


def test_nan_component_reports_finiteness():
    row = np.array([[1.0, np.nan, 3.0, 4.0, 5.0, 6.0]])
    with pytest.raises(ValueError, match="finite"):
        components_to_sym3(row, CANONICAL, kind="stress")
