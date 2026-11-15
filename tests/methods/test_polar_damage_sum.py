"""PDS_SWT / PDS_FS (N59): Albinmousa & Al Hussain (2021), Forces in Mechanics 4:100025,
doi:10.1016/j.finmec.2021.100025, checked against the open-access PDF supplied by the maintainer.

Eqs. (7)-(14), pp.3-5. The paper's examples are strain-controlled low-cycle tests with measured
hysteresis loops (Figs. 4 and 7) and fitted constants (Tables 3 and 5); no per-test value is printed,
and the elastic repository version cannot reproduce the measured loops. The evidence is therefore
the definitions: the in-plane transformations (Eqs. 7a-7d), the polar area (Eqs. 10-12) and the
closed form of a uniaxial elastic state.
"""

import numpy as np
import pytest

from multifat.models.errors import LoadingScopeError
from multifat.models.integral_methods.polar_damage_sum import (
    PolarDamageSumFS,
    PolarDamageSumSWT,
    polar_damage,
    polar_damage_sum,
)

E, NU = 200000.0, 0.3
ZERO = np.zeros((3, 3))


def strain(S):
    """Linear-elastic strain of a stress tensor."""
    return ((1 + NU) * S - NU * np.trace(S) * np.eye(3)) / E


def plane(sx, sy, txy):
    return np.array([[sx, txy, 0.0], [txy, sy, 0.0], [0.0, 0.0, 0.0]])


def test_eq7_transformations():
    """sigma_phi, eps_phi, gamma_phi/2 of Eqs. (7a), (7c), (7d) for sigma_x = 0 (phi from y)."""
    sy, txy = 300.0, 120.0
    S = plane(0.0, sy, txy)
    Ee = strain(S)
    phi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    sigma_7a = sy / 2 + sy / 2 * np.cos(2 * phi) - txy * np.sin(2 * phi)
    ex, ey, gxy = Ee[0, 0], Ee[1, 1], 2 * Ee[0, 1]
    eps_7c = (ex + ey) / 2 - (ex - ey) / 2 * np.cos(2 * phi) - gxy / 2 * np.sin(2 * phi)
    half_gam_7d = gxy / 2 * np.cos(2 * phi) - (ex - ey) / 2 * np.sin(2 * phi)
    # the repository plane at angle psi has normal (cos psi, sin psi): phi of Fig. 1 is psi - 90 deg
    psi = phi + np.pi / 2
    n = np.column_stack([np.cos(psi), np.sin(psi), 0 * psi])
    t = np.column_stack([-np.sin(psi), np.cos(psi), 0 * psi])
    assert np.einsum("ni,ij,nj->n", n, S, n) == pytest.approx(sigma_7a)
    assert np.einsum("ni,ij,nj->n", n, Ee, n) == pytest.approx(eps_7c)
    assert np.abs(np.einsum("ni,ij,nj->n", t, Ee, n)) == pytest.approx(np.abs(half_gam_7d))


def test_uniaxial_closed_form_and_critical_plane_value():
    """R = -1 tension: D_SWT(phi) = sigma_a eps_a cos^2 phi |cos^2 phi - nu sin^2 phi| (normal along y)."""
    sa = 250.0
    S = plane(0.0, sa, 0.0)
    D = polar_damage("SWT", ZERO, S, ZERO, strain(S), ZERO)
    psi = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    c2, s2 = np.sin(psi) ** 2, np.cos(psi) ** 2  # angle from the y axis
    expected = sa * (sa / E) * c2 * np.abs(c2 - NU * s2)
    assert D == pytest.approx(expected, rel=1e-12, abs=1e-18)
    assert D.max() == pytest.approx(
        sa * sa / E
    )  # the critical-plane SWT value is the peak (Fig. 2b)
    x = np.linspace(0, 2 * np.pi, 400000, endpoint=False)
    c2x, s2x = np.sin(x) ** 2, np.cos(x) ** 2
    area = 0.5 * np.mean((sa * sa / E * c2x * np.abs(c2x - NU * s2x)) ** 2) * 2 * np.pi
    assert polar_damage_sum("SWT", ZERO, S, ZERO, strain(S), ZERO) == pytest.approx(area, rel=1e-6)


def test_riemann_sum_converges_for_out_of_phase():
    """Eq. (11) -> Eq. (12): 720 planes agree with a dense sum for a 90-degree tension-torsion path."""
    S_sin, S_cos = plane(0.0, 250.0, 0.0), plane(0.0, 0.0, 140.0)
    args = ("FS", plane(0.0, 80.0, 0.0), S_sin, S_cos, strain(S_sin), strain(S_cos), 1.2, 600.0)
    coarse = polar_damage_sum(*args)
    psi = np.linspace(0, 2 * np.pi, 200000, endpoint=False)
    n = np.column_stack([np.cos(psi), np.sin(psi), 0 * psi])
    t = np.column_stack([-np.sin(psi), np.cos(psi), 0 * psi])

    def amp(A, B, u, v):
        return np.hypot(np.einsum("ni,ij,nj->n", u, A, v), np.einsum("ni,ij,nj->n", u, B, v))

    smax = np.einsum("ni,ij,nj->n", n, plane(0.0, 80.0, 0.0), n) + amp(S_sin, S_cos, n, n)
    D = 2 * amp(strain(S_sin), strain(S_cos), t, n) * (1 + 1.2 * smax / 600.0)
    assert coarse == pytest.approx(0.5 * np.mean(np.maximum(D, 0) ** 2) * 2 * np.pi, rel=1e-5)


def test_fs_constant_raises_the_tension_side():
    S = plane(0.0, 200.0, 0.0)
    low = polar_damage_sum("FS", ZERO, S, ZERO, strain(S), ZERO, 0.0, 600.0)
    high = polar_damage_sum("FS", plane(0.0, 100.0, 0.0), S, ZERO, strain(S), ZERO, 1.0, 600.0)
    assert high > low > 0.0


def test_models_and_scope():
    S = plane(0.0, 200.0, 90.0)
    harm = {
        "S_mean": ZERO,
        "S_sin": S,
        "S_cos": ZERO,
        "E_mean": ZERO,
        "E_sin": strain(S),
        "E_cos": ZERO,
    }
    swt = PolarDamageSumSWT().evaluate_case(
        S0=S, E0=None, S1=-S, E1=None, params={}, harmonics=harm
    )
    assert swt.values[""] == pytest.approx(polar_damage_sum("SWT", ZERO, S, ZERO, strain(S), ZERO))
    fs = PolarDamageSumFS().evaluate_case(
        S0=S, E0=None, S1=-S, E1=None, params={"k_PDS": 0.5, "Sy": 500.0}, harmonics=harm
    )
    assert fs.values[""] > 0
    S3 = S.copy()
    S3[2, 2] = 50.0
    with pytest.raises(LoadingScopeError, match="plane stress"):
        PolarDamageSumSWT().evaluate_case(
            S0=S3, E0=None, S1=-S3, E1=None, params={}, harmonics={**harm, "S_sin": S3}
        )
