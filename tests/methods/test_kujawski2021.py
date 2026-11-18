"""KUJAWSKI2021 (N26): Kujawski (2021), Fatigue Fract Eng Mater Struct 44:2202-2209,
doi:10.1111/ffe.13488, checked against the PDF supplied by the maintainer.

Eqs. (11a)-(13b), pp.2205-2207. The paper validates mSWT only with figures (Figs. 8-11, ductile
iron and 7075-T651 data of other authors) and prints no table, so the checks are closed forms
(maintainer decision): the SWT branch for sigma_m >= 0, the |sigma_m|/3 augmentation for
sigma_m < 0, the stress-based forms of Eqs. (12a)-(12b) and (13a)-(13b) in uniaxial elastic
loading, and the critical-plane extension (identical to SWT without a compressive normal mean).
"""

import numpy as np
import pytest

from multifat.models.cp_methods.kujawski2021 import Kujawski2021, mswt_stress
from multifat.models.cp_methods.swt import SWT
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

E_MOD, NU = 71700.0, 0.33  # 7075-T651, the material of Figs. 8, 10 and 11
R = build_R_list()


def _case(sigma_a=0.0, sigma_m=0.0, tau_a=0.0, delta=0.0):
    """Axial (yy) stress and torsion (xy) with elastic strains, in multifat's harmonic form."""
    d = np.radians(delta)

    def t(x, y):
        return np.array([[0.0, y, 0.0], [y, x, 0.0], [0.0, 0.0, 0.0]])

    def eps(s):
        return ((1 + NU) * s - NU * np.trace(s) * np.eye(3)) / E_MOD

    S_mean, S_sin, S_cos = (
        t(sigma_m, 0.0),
        t(sigma_a, tau_a * np.cos(d)),
        t(0.0, -tau_a * np.sin(d)),
    )
    return dict(
        S0=S_mean + S_sin, S1=S_mean - S_sin, E0=eps(S_mean + S_sin), E1=eps(S_mean - S_sin),
        R_list=R, params={"E": E_MOD, "nu": NU},
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos,
                   "E_mean": eps(S_mean), "E_sin": eps(S_sin), "E_cos": eps(S_cos)},
    )  # fmt: skip


def _dp(name, **load):
    (model,) = get_models([name])
    return model.evaluate_case(**_case(**load)).values


@pytest.mark.parametrize(
    "sigma_max, sigma_m, expected",
    [
        (300.0, 100.0, 300.0),  # Eq. (11a): sigma_m >= 0, original SWT
        (200.0, 0.0, 200.0),  # R = -1: both branches agree
        (110.0, -90.0, 140.0),  # Eq. (11b): 110 + 90/3
        (-50.0, -150.0, 0.0),  # compression-compression, bracket 100 - 2*150/3 = 0
        (-20.0, -120.0, 20.0),  # compression-compression with sigma_max < 0 but damaging (Fig. 10)
        (-100.0, -300.0, 0.0),  # bracket < 0: clipped (repository choice)
    ],
)
def test_eq11_stress_factor(sigma_max, sigma_m, expected):
    assert mswt_stress(sigma_max, sigma_m) == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("sigma_a, sigma_m", [(200.0, 0.0), (150.0, 100.0), (200.0, -90.0)])
def test_eq12_13_uniaxial_equivalent_amplitude(sigma_a, sigma_m):
    """Uniaxial elastic loading: sqrt(E W_D) is sigma_ar of Eqs. (12a)-(12b) and (13a)-(13b)."""
    W_D = _dp("KUJAWSKI2021", sigma_a=sigma_a, sigma_m=sigma_m)[""]
    sigma_max = sigma_m + sigma_a
    bracket = sigma_max + abs(sigma_m) / 3 if sigma_m < 0 else sigma_max
    assert np.sqrt(E_MOD * W_D) == pytest.approx(np.sqrt(bracket * sigma_a), rel=1e-9)  # Eq. 12
    eps_a = sigma_a / E_MOD
    assert np.sqrt(E_MOD * W_D) == pytest.approx(np.sqrt(E_MOD * eps_a * bracket), rel=1e-9)


def test_negative_mean_is_more_damaging_than_swt():
    """The purpose of Eq. (11b): SWT is nonconservative for sigma_m < 0 (R < -1)."""
    swt = _dp("SWT", sigma_a=200.0, sigma_m=-90.0)[""]
    mswt = _dp("KUJAWSKI2021", sigma_a=200.0, sigma_m=-90.0)[""]
    assert mswt / swt == pytest.approx(140.0 / 110.0, rel=1e-9)


@pytest.mark.parametrize(
    "load",
    [
        dict(sigma_a=200.0),
        dict(tau_a=120.0),
        dict(sigma_a=150.0, sigma_m=80.0, tau_a=90.0, delta=0.0),
        dict(sigma_a=150.0, sigma_m=80.0, tau_a=90.0, delta=90.0),
        dict(sigma_a=150.0, tau_a=100.0, delta=45.0),
    ],
)
def test_identical_to_swt_without_compressive_normal_mean(load):
    """Critical-plane extension: sigma_n,m >= 0 on every plane, so both variants equal SWT."""
    assert _dp("KUJAWSKI2021", **load) == pytest.approx(_dp("SWT", **load), rel=1e-12)


def test_plane_level_equals_swt_for_nonnegative_means():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(50, 3, 3)) * 100
    a = 0.5 * (a + a.transpose(0, 2, 1))
    mean = np.zeros((50, 3, 3))
    mean[:, 2, 2] = np.abs(rng.normal(size=50)) * 50
    S0, S1, E0, E1 = mean + a, mean - a, a / E_MOD, -a / E_MOD
    kw = dict(S0r=S0, S1r=S1, E0r=E0, E1r=E1, params={})
    k_dp, k_m = Kujawski2021().evaluate_on_planes(**kw)
    s_dp, s_m = SWT().evaluate_on_planes(**kw)
    np.testing.assert_array_equal(k_dp, s_dp)
    np.testing.assert_array_equal(k_m, s_m)
