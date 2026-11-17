"""Ellyin & Kujawski (1993), ASTM STP 1191, doi:10.1520/STP24795S (Zotero K3WD859Q), roadmap G16.

Checked on crops _audit-trail/phase2/crops/K3WD859Q/p4-p8: Eqs. (5b), (9)-(18), Fig. 2. The paper
prints no per-test values; the checks are its own limits: p_bar = 1 (uniaxial), 1 + nu (reversed
torsion), 1 - nu (equibiaxial) (p.60), Eq. (13) -> nu_e for elastic strains, and Eq. (18)
sigma_a / sigma_-1 = (1 + sigma_m / (2 sigma_a))^(-1/2) for alpha = 1/2.
"""

import numpy as np
import pytest

from multifat.models.energy_based_methods.ellyin_kujawski import (
    EllyinKujawski,
    _elastic_strain,
    constraint_factor,
    ek_damage,
)

E, NU = 206000.0, 0.3
Z = np.zeros((3, 3))


def _t(xx=0.0, yy=0.0, xy=0.0):
    return np.array([[xx, xy, 0.0], [xy, yy, 0.0], [0.0, 0.0, 0.0]])


def _psi(mean, amp, phase_cos=None, **kw):
    cos = Z if phase_cos is None else phase_cos
    return ek_damage(
        mean, amp, cos, _elastic_strain(amp, E, NU), _elastic_strain(cos, E, NU), E, NU, **kw
    )


@pytest.mark.parametrize(
    "S, p", [(_t(yy=300), 1.0), (_t(xy=200), 1.0 + NU), (_t(xx=250, yy=250), 1.0 - NU)]
)
def test_constraint_factor_limits(S, p):
    assert constraint_factor(_elastic_strain(S, E, NU), Z, NU) == pytest.approx(p)


@pytest.mark.parametrize("ratio", [-0.8, 0.0, 0.5, 1.0, 2.0, 4.0])
def test_eq18_uniaxial_mean_stress(ratio):
    f_m1 = 400.0
    sa = f_m1 / np.sqrt(1.0 + 0.5 * ratio)  # Eq. (18)
    psi = _psi(_t(yy=ratio * sa), _t(yy=sa))
    assert psi == pytest.approx(_psi(Z, _t(yy=f_m1)))


@pytest.mark.parametrize("n, ratio", [(1.0, np.sqrt((1 + NU) / 3)), (0.0, 1 / np.sqrt(3))])
def test_torsion_to_tension_strength_ratio(n, ratio):
    s = 400.0
    t = ratio * s  # equal psi in reversed torsion and reversed tension
    assert _psi(Z, _t(xy=t), n=n) == pytest.approx(_psi(Z, _t(yy=s), n=n))


def test_mean_shear_has_no_effect_and_phase_origin_invariance():
    a = _psi(Z, _t(yy=300, xy=150))
    assert _psi(_t(xy=100), _t(yy=300, xy=150)) == pytest.approx(a)
    ph = np.deg2rad(37.0)
    S_sin, S_cos = _t(yy=300 * np.cos(ph), xy=150 * np.cos(ph)), _t(
        yy=300 * np.sin(ph), xy=150 * np.sin(ph)
    )
    assert _psi(Z, S_sin, S_cos) == pytest.approx(a, rel=1e-4)


def test_case_model():
    H = {"S_mean": Z, "S_sin": _t(yy=300), "S_cos": Z}
    H.update(E_sin=_elastic_strain(H["S_sin"], E, NU), E_cos=Z)
    res = EllyinKujawski().evaluate_case(
        S0=Z, E0=Z, S1=Z, E1=Z, params={"E": E, "nu": NU}, harmonics=H
    )
    assert res.values[""] == pytest.approx((600.0**2) * (1 + NU) / (3 * E))


def test_plastic_term_vanishes_for_stiff_cyclic_curve():
    a = _psi(Z, _t(yy=300, xy=120))
    b = _psi(Z, _t(yy=300, xy=120), K_prime=1e12, n_prime=0.15)
    assert b == pytest.approx(a, rel=1e-9)


def test_plastic_uniaxial_energy_and_constraint():
    from multifat.models.energy_based_methods.ellyin_kujawski import effective_poisson

    Kp, npr, s = 1100.0, 0.12, 420.0
    ep = (s / Kp) ** (1 / npr)
    E_star = 3 * E / (2 * (1 + NU))
    psi = _psi(Z, _t(yy=s), K_prime=Kp, n_prime=npr)
    expected = (2 * s) ** 2 / (2 * E_star) + (1 - npr) / (2 * (1 + npr)) * (2 * s) * (2 * ep)
    assert psi == pytest.approx(expected)  # p_bar stays 1 in uniaxial loading
    ee = s / E
    eps_a, eps_t = ee + ep, -NU * ee - 0.5 * ep  # uniaxial elastic + plastic strains
    eq13 = (0.5 * (1 - NU) * (eps_a + eps_t) + (NU - 0.5) * (ee - NU * ee)) / (
        (1 - NU) * (eps_a + eps_t) + (NU - 0.5) * (ee - NU * ee)
    )
    assert effective_poisson(ee, ep, NU) == pytest.approx(eq13)


def test_eq_18_mean_stress_range_is_the_loading_scope():
    """Eq. (18) holds for -1.06 <= sigma_m/sigma_a <= 4.4 (p.61); outside, the case is out of
    scope (maintainer decision, 17 Nov 2026), the bounds themselves are inside."""
    import numpy as np
    import pytest

    from multifat.models.energy_based_methods.ellyin_kujawski import ek_damage
    from multifat.models.errors import LoadingScopeError

    E, nu, s_a = 200000.0, 0.3, 100.0
    S_sin = np.diag([0.0, s_a, 0.0])
    E_sin = np.diag([-nu, 1.0, -nu]) * s_a / E
    Z = np.zeros((3, 3))
    for ratio in (-1.06, 0.0, 4.4):
        S_mean = np.diag([0.0, ratio * s_a, 0.0])
        assert ek_damage(S_mean, S_sin, Z, E_sin, Z, E, nu) > 0.0
    for ratio in (-1.07, 4.41):
        with pytest.raises(LoadingScopeError):
            ek_damage(np.diag([0.0, ratio * s_a, 0.0]), S_sin, Z, E_sin, Z, E, nu)
