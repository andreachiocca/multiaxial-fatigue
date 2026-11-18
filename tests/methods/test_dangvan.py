"""DANGVAN, Dang Van criterion as restated by Boehme, Papuga & Lange (2026) FFEMS,
doi:10.1111/ffe.70244 (Zotero 5T8PJV87): Eq. (3) p.4 with a = kappa, b = 3 - 1.5 kappa,
kappa = s_-1/t_-1 (Eq. 4); Eq. (16) closed forms p.8; Table 4 DV column p.8
(crop _audit-trail/phase2/crops/5T8PJV87/p8.png). Documented variant of Dang Van (1973), the
registered reference; neither it nor the later restatements (Dang Van 1993, 1999; background
only) give closed-form a, b.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.cp_methods.dangvan import dang_van_equivalent
from multifat.models.registry import get_models

Z = np.zeros((3, 3))


@pytest.mark.parametrize("kappa", [1.0, 1.3, np.sqrt(3), 2.0])
def test_eq16_closed_forms(kappa):
    sa, t = 100.0, 250.0
    s = kappa * t
    bi = dang_van_equivalent(Z, np.diag([sa, sa, 0]), Z, s, t)
    tri = dang_van_equivalent(Z, np.diag([sa, sa, sa]), Z, s, t)
    assert bi == pytest.approx(2 * sa * (1 - kappa / 4), rel=1e-9)
    assert tri == pytest.approx(3 * sa * (1 - kappa / 2), abs=1e-9 * sa)


@pytest.mark.parametrize("s,t", [(433.0, 344.0), (241.0, 144.0), (500.0, 250.0)])
def test_calibration_limits(s, t):
    assert dang_van_equivalent(Z, np.diag([s, 0, 0]), Z, s, t) == pytest.approx(s, rel=1e-9)
    shear = np.array([[0, t, 0], [t, 0, 0], [0, 0, 0.0]])
    assert dang_van_equivalent(Z, shear, Z, s, t) == pytest.approx(s, rel=1e-9)


def test_rotation_and_phase_invariance():
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    A = np.array([[200.0, 80, 0], [80, -50, 30], [0, 30, 10]])
    B = np.array([[0.0, 120, 0], [120, 40, 0], [0, 0, -20]])
    M = np.diag([100.0, 20, 0])
    ref = dang_van_equivalent(M, A, B, 400.0, 260.0)
    rot = lambda X: Q @ X @ Q.T  # noqa: E731
    assert dang_van_equivalent(rot(M), rot(A), rot(B), 400.0, 260.0) == pytest.approx(ref, rel=1e-6)
    c, s = np.cos(0.7), np.sin(0.7)
    shifted = dang_van_equivalent(M, A * c - B * s, A * s + B * c, 400.0, 260.0)
    assert shifted == pytest.approx(ref, rel=1e-6)


def test_table4_dv_statistics():
    (dv,) = get_models(["DANGVAN"])
    err = []
    for r in load_table("bohme2026_table2").itertuples():
        mean, sin = Z.copy(), Z.copy()
        mean[0, 0] = r.sigma_xm
        sin[0, 1] = sin[1, 0] = r.tau_xya
        h = dict(S_mean=mean, S_sin=sin, S_cos=Z)
        v = dv.evaluate_case(S0=Z, E0=Z, S1=Z, E1=Z, params=dict(Sigm1=r.s_m1, Taum1=r.t_m1),
                             harmonics=h).values[""]  # fmt: skip
        err.append(100 * (v - 1))
    e = np.array(err)
    assert len(e) == 18
    # printed to 0.1: Max 20.06 and SD 11.84 print as 20.0 and 11.9 (truncation/rounding)
    assert e.max() == pytest.approx(20.0, abs=0.07)
    assert e.min() == pytest.approx(-18.6, abs=0.05)
    assert e.max() - e.min() == pytest.approx(38.7, abs=0.1)
    assert e.mean() == pytest.approx(3.6, abs=0.05)
    assert e.std() == pytest.approx(11.9, abs=0.07)
    assert 100 * np.mean(abs(e) <= 10) == pytest.approx(44.4, abs=0.05)
