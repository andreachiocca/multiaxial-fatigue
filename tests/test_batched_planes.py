"""Batched plane evaluation must equal the per-plane scalar path bit for bit."""

import numpy as np
import pytest

from multifat.eval.evaluator import best_cp_for_case, evaluate_planes, prepare_planes
from multifat.models.registry import available_cp_methods
from multifat.models.utils.batch import pypow2
from multifat.planes.sampling import build_R_list

PARAMS = dict(
    E=200000.0,
    nu=0.3,
    Su=900.0,
    Sy=600.0,
    Sigm1=400.0,
    Taum1=250.0,
    k_FS=0.6,
    k_FI=0.3,
    k_MGSE=0.8,
    a_MSWT=0.4,
    L_LI=2.0,
    sigma_f=1100.0,
    tau_f=640.0,
    b_fat=-0.09,
    b0=-0.1,
    epsilon_f=0.4,
    c_fat=-0.6,
    gamma_f=0.7,
    c0=-0.6,
)


class _ScalarOnly:
    """Hide ``evaluate_on_planes`` so the evaluator takes the per-plane fallback."""

    def __init__(self, method):
        self._m, self.name = method, method.name
        self.np_shear_rule = getattr(method, "np_shear_rule", "rss")
        self.ext_variant = getattr(method, "ext_variant", True)

    def evaluate_on_plane(self, **kw):
        return self._m.evaluate_on_plane(**kw)


def _case(rng, harmonic):
    sym = lambda a: a + a.T  # noqa: E731
    mean, sin, cos = (sym(rng.normal(size=(3, 3))) * s for s in (60.0, 120.0, 80.0))
    if not harmonic:
        cos = np.zeros((3, 3))
    strain = lambda s: ((1 + 0.3) * s - 0.3 * np.trace(s) * np.eye(3)) / 200000.0  # noqa: E731
    plastic = sym(rng.normal(size=(3, 3))) * 1e-4
    e_mean, e_sin, e_cos = strain(mean), strain(sin) + plastic, strain(cos)
    return dict(
        S0=mean + sin,
        S1=mean - sin,
        E0=e_mean + e_sin,
        E1=e_mean - e_sin,
        S_sin=sin,
        S_cos=cos,
        E_sin=e_sin,
        E_cos=e_cos,
    )


@pytest.mark.parametrize("name", sorted(available_cp_methods()))
@pytest.mark.parametrize("harmonic", [False, True])
def test_batch_equals_scalar_bitwise(name, harmonic):
    method = available_cp_methods()[name]
    R_list = build_R_list()
    rng = np.random.default_rng(7)
    for _ in range(4):
        case = _case(rng, harmonic)
        planes = prepare_planes(
            **case, R_list=R_list, tau_max=True, strain_ranges=True, E=2e5, nu=0.3
        )
        d_b, m_b = evaluate_planes(method, planes, PARAMS)
        d_s, m_s = evaluate_planes(_ScalarOnly(method), planes, PARAMS)
        assert np.array_equal(d_b, d_s) and np.array_equal(m_b, m_s)
        a = best_cp_for_case(**case, R_list=R_list, method=method, params=PARAMS)
        b = best_cp_for_case(**case, R_list=R_list, method=_ScalarOnly(method), params=PARAMS)
        assert (a.ext, a.by_metric) == (b.ext, b.by_metric)


def test_non_orthonormal_plane_rejected():
    R = build_R_list()[:3].copy()
    R[1] *= 1.01
    z = np.zeros((3, 3))
    with pytest.raises(ValueError, match="orthonormal"):
        prepare_planes(S0=z, E0=z, S1=z, E1=z, R_list=R)


def test_pypow2_matches_cpython_float_pow():
    edge = np.array(
        [0.0, -0.0, 1.0, -1.0, np.inf, -np.inf, np.nan, 1e-160, 1e-310, -1e-310, 1e100, -1e100]
    )
    rng = np.random.default_rng(0)
    scales = [rng.normal(0.0, s, size=2500) for s in (1e-10, 1.0, 1e5, 1e100)]
    samples = np.concatenate([edge, *scales])

    for arr in (
        samples,
        samples.reshape(4, -1),
        samples.reshape(2, 2, -1),
        np.array(-3.25),
        np.array([1, -2, 3], dtype=int),
        np.empty((0,), dtype=float),
        np.empty((2, 0), dtype=float),
    ):
        got = pypow2(arr)
        expected_flat = [float(v) ** 2 for v in np.asarray(arr, dtype=float).ravel().tolist()]
        expected = np.array(expected_flat, dtype=float).reshape(np.asarray(arr).shape)
        assert isinstance(got, np.ndarray)
        assert got.dtype == np.float64
        assert got.shape == expected.shape
        assert np.array_equal(got, expected, equal_nan=True)
        assert np.array_equal(np.signbit(got), np.signbit(expected))


def test_pypow2_overflow_matches_cpython():
    with pytest.raises(OverflowError):
        pypow2(np.array([1e200]))
