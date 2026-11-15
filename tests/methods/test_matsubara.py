"""MATSUBARA (N09): Matsubara, Hayashida & Kano (2018), Int J Fatigue, accepted manuscript,
doi:10.1016/j.ijfatigue.2017.12.001, checked against the PDF supplied by the maintainer.

Eqs. (7)-(9), (13)-(15), pp.7-11. The paper reports its verification only as histograms of the
error index (Figs. 6-7) and life scatter plots (Figs. 8-12); no per-test value is printed, so the
evidence here is the closed forms that tie Eq. (15) to Eq. (9) and to the uniaxial S-N curve.
"""

import numpy as np
import pytest

from multifat.models.errors import LoadingScopeError, MaterialScopeError, MissingParamsError
from multifat.models.invariant_methods.matsubara import (
    Matsubara,
    equivalent_shear_amplitude,
    max_range_normal_stress,
    true_fracture_strength,
    unified_equivalent_shear,
)

ZERO = np.zeros((3, 3))
SQRT3 = np.sqrt(3.0)


def ax_tor(sigma, tau):
    return np.array([[sigma, tau, 0.0], [tau, 0.0, 0.0], [0.0, 0.0, 0.0]])


@pytest.mark.parametrize("sigma_a", [100.0, 250.0, 400.0])
def test_uniaxial_reversed_maps_onto_itself(sigma_a):
    """On sqrt(J_2,amp)/S_max = 1/sqrt3 (R = -1 uniaxial), Eq. (15) gives sigma_a/sqrt3 (Fig. 4)."""
    for sigma_T in (700.0, 1300.0):
        value = unified_equivalent_shear(ZERO, ax_tor(sigma_a, 0), ZERO, sigma_T)
        assert value == pytest.approx(sigma_a / SQRT3, rel=1e-12)


@pytest.mark.parametrize("sigma_m,tau_m", [(0, 0), (120, 0), (0, 80), (-60, 40)])
def test_eq9_is_eq15_at_the_fatigue_limit(sigma_m, tau_m):
    """States on the Eq. (9) line have sqrt(J_2,amp,eq) = sigma_w/sqrt3."""
    sigma_w, sigma_T = 314.0, 1.84 * 680.0
    S_mean, S_dir = ax_tor(sigma_m, tau_m), ax_tor(1.0, 0.6)
    j2_unit = equivalent_shear_amplitude(S_dir, ZERO)
    # solve Eq. (9) for the amplitude scale c of S_dir (S_max is linear in c for this path)
    lo, hi = 1e-6, 2000.0
    for _ in range(200):
        c = 0.5 * (lo + hi)
        lhs = c * j2_unit + sigma_w / SQRT3 / (sigma_T - sigma_w) * max_range_normal_stress(
            S_mean, c * S_dir, ZERO
        )
        lo, hi = (c, hi) if lhs < sigma_w / SQRT3 * sigma_T / (sigma_T - sigma_w) else (lo, c)
    value = unified_equivalent_shear(S_mean, c * S_dir, ZERO, sigma_T)
    assert value == pytest.approx(sigma_w / SQRT3, rel=1e-9)


def test_mamiya_amplitude_of_one_harmonic():
    """Eq. (13) with the basis of Eq. (8): half-ranges of S_1..S_5 equal the component amplitudes."""
    S_sin, S_cos = ax_tor(200, 50), ax_tor(0, 90)
    t = np.linspace(0, 2 * np.pi, 20001)
    S = S_sin[None] * np.sin(t)[:, None, None] + S_cos[None] * np.cos(t)[:, None, None]
    tr = np.trace(S, axis1=1, axis2=2)
    s = S - tr[:, None, None] / 3 * np.eye(3)
    comps = [
        np.sqrt(3) / 2 * s[:, 0, 0],
        (s[:, 1, 1] - s[:, 2, 2]) / 2,
        s[:, 0, 1],
        s[:, 1, 2],
        s[:, 2, 0],
    ]
    mamiya = np.sqrt(sum(((c.max() - c.min()) / 2) ** 2 for c in comps))
    assert equivalent_shear_amplitude(S_sin, S_cos) == pytest.approx(mamiya, rel=1e-6)


def test_smax_brute_force_and_ties():
    rng = np.random.default_rng(3)
    normals = rng.normal(size=(200000, 3))
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    S_mean, S_sin, S_cos = ax_tor(80, 30), ax_tor(200, 60), ax_tor(0, 110)
    a = np.hypot(*(np.einsum("ni,ij,nj->n", normals, T, normals) for T in (S_sin, S_cos)))
    i = np.argmax(a)
    brute = np.einsum("i,ij,j", normals[i], S_mean, normals[i]) + a[i]
    assert max_range_normal_stress(S_mean, S_sin, S_cos) == pytest.approx(brute, rel=1e-3)
    # pure torsion with a mean shear: the two max-range planes tie; the larger maximum is taken
    assert max_range_normal_stress(ax_tor(0, 40), ax_tor(0, 150), ZERO) == pytest.approx(190.0)
    assert max_range_normal_stress(ax_tor(100, 0), ax_tor(150, 0), ZERO) == pytest.approx(250.0)


def test_sigma_T_rules():
    assert true_fracture_strength({"sigma_T": 1500.0, "Su": 2000.0}) == 1500.0
    assert true_fracture_strength({"Su": 680.0}) == pytest.approx(1.84 * 680.0)
    with pytest.raises(MaterialScopeError, match="400 <= Su <= 1000"):
        true_fracture_strength({"Su": 343.0})
    with pytest.raises(MissingParamsError):
        true_fracture_strength({})


def test_eq15_undefined_is_out_of_scope():
    with pytest.raises(LoadingScopeError, match="Eq. \\(15\\)"):
        unified_equivalent_shear(ax_tor(900, 0), ax_tor(100, 0), ZERO, 900.0)


def test_model_case():
    out = Matsubara().evaluate_case(
        S0=None, E0=None, S1=None, E1=None, params={"Su": 680.0},
        harmonics={"S_mean": ZERO, "S_sin": ax_tor(300, 0), "S_cos": ZERO},
    )  # fmt: skip
    assert out.values[""] == pytest.approx(300 / SQRT3)
