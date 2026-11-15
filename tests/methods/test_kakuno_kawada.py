"""KAKUNO (G04): Kakuno & Kawada (1979), Fatigue Engng Mater Struct 2:229-236,
doi:10.1111/j.1460-2695.1979.tb01358.x, checked against the PDF supplied by the maintainer.

Eqs. (4)-(12), pp.230-234. The paper validates the criterion only on figures (Gough's data,
Figs. 2-6); no table of values is printed.

Findings:
- Eqs. (7)-(9) follow from Eq. (6) at the three reference limits, and Eq. (10) is Eq. (6) squared
  and normalised: both are checked here to machine precision;
- Fig. 5 (m_1 = 0.46, k_1 = 0.83, t_w/s_wb = 0.617): the Eq. (10) curve meets the tau axis at
  0.557 (read from the figure), against 0.559 computed; it meets the sigma axis at about 0.93,
  against 0.906 computed (strict xfail, unexplained: the caption values give one curve whose
  two intercepts cannot both match the drawing);
- Eq. (12): B = C would require s_wb/s_up = 1 - s_wb/(2 sqrt3 t_w), i.e. k_1 = 1 for
  s_wb/t_w = sqrt3, as the paper states.
"""

import numpy as np
import pytest

from multifat.models.errors import LoadingScopeError
from multifat.models.invariant_methods.kakuno_kawada import (
    KakunoKawada,
    kakuno_constants,
    kakuno_index,
)

ZERO = np.zeros((3, 3))


def bt(sigma, tau):
    """Bending sigma (xx) and torsion tau (xy) of a round bar surface."""
    return np.array([[sigma, tau, 0.0], [tau, 0.0, 0.0], [0.0, 0.0, 0.0]])


def eq10(sigma_a, tau_a, sigma_m, s_wb, t_w, k_1):
    """Left-hand side of Eq. (10)."""
    p = 1 - sigma_m / s_wb * (1 - k_1) / k_1
    g = 2 / np.sqrt(3) * s_wb / t_w
    x = sigma_a / s_wb
    return (tau_a / t_w) ** 2 / p**2 + x**2 * (g - 1) / p**2 + x * (2 - g) / p


@pytest.mark.parametrize(
    "s_wb,t_w,k_1", [(400.0, 250.0, 0.83), (313.9, 196.2, 0.9), (200, 160, 0.7)]
)
def test_reference_limits_eqs_7_9(s_wb, t_w, k_1):
    f_0 = 2 * k_1 * s_wb
    assert kakuno_index(ZERO, bt(0, t_w), ZERO, s_wb, t_w, f_0) == pytest.approx(1.0, abs=1e-12)
    assert kakuno_index(ZERO, bt(s_wb, 0), ZERO, s_wb, t_w, f_0) == pytest.approx(1.0, abs=1e-12)
    half = k_1 * s_wb  # pulsating bending: sigma_a = sigma_m = s_up/2
    assert kakuno_index(bt(half, 0), bt(half, 0), ZERO, s_wb, t_w, f_0) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "sigma_a,sigma_m", [(0.0, 0.0), (120.0, 0.0), (200.0, 150.0), (80.0, -60.0)]
)
def test_eq10_is_eq4_for_bending_torsion(sigma_a, sigma_m):
    s_wb, t_w, k_1 = 420.0, 260.0, 0.83
    p = 1 - sigma_m / s_wb * (1 - k_1) / k_1
    g = 2 / np.sqrt(3) * s_wb / t_w
    x = sigma_a / s_wb
    tau_a = t_w * p * np.sqrt(1 - x**2 * (g - 1) / p**2 - x * (2 - g) / p)  # on the Eq. (10) curve
    assert eq10(sigma_a, tau_a, sigma_m, s_wb, t_w, k_1) == pytest.approx(1.0)
    index = kakuno_index(bt(sigma_m, 0), bt(sigma_a, tau_a), ZERO, s_wb, t_w, 2 * k_1 * s_wb)
    assert index == pytest.approx(1.0, abs=1e-12)


def test_torsional_mean_does_not_enter():
    s_wb, t_w = 420.0, 260.0
    plain = kakuno_index(ZERO, bt(150, 100), ZERO, s_wb, t_w)  # no Sig0 needed: tr S_m = 0
    assert kakuno_index(bt(0, 120), bt(150, 100), ZERO, s_wb, t_w) == pytest.approx(plain)


FIG5 = dict(s_wb=1.0, t_w=0.617, k_1=0.83, sigma_m=0.46)


def _fig5_curve(x):
    """tau_a/s_wb on the Eq. (10) curve of Fig. 5 at sigma_a/s_wb = x."""
    c = FIG5
    p = 1 - c["sigma_m"] * (1 - c["k_1"]) / c["k_1"]
    g = 2 / np.sqrt(3) / c["t_w"]
    return c["t_w"] * p * np.sqrt(max(0.0, 1 - x**2 * (g - 1) / p**2 - x * (2 - g) / p))


def test_fig5_tau_intercept():
    assert _fig5_curve(0.0) == pytest.approx(0.557, abs=0.005)  # read from Fig. 5


@pytest.mark.xfail(strict=True, reason="Fig. 5 curve meets the sigma axis at about 0.93; Eq. (10) "
                   "with the caption values gives 0.906 (unexplained)")  # fmt: skip
def test_fig5_sigma_intercept():
    xs = np.linspace(0.8, 1.0, 20001)
    intercept = xs[np.argmax([_fig5_curve(x) == 0.0 for x in xs])]
    assert intercept == pytest.approx(0.93, abs=0.01)


def test_eq12_b_equals_c():
    """B = C (Findley-like maximum octahedral normal stress) forces s_up = 2 s_wb for s_wb/t_w = sqrt3."""
    s_wb = 400.0
    t_w = s_wb / np.sqrt(3)
    k_1 = 1 / (2 * (1 - s_wb / (2 * np.sqrt(3) * t_w)))  # Eq. (12) solved for s_up/(2 s_wb)
    _, B, C = kakuno_constants(s_wb, t_w, 2 * k_1 * s_wb)
    assert k_1 == pytest.approx(1.0) and B == pytest.approx(C, abs=1e-12)


def test_scope():
    model, params = KakunoKawada(), {"Sigm1": 420.0, "Taum1": 260.0}
    out_of_phase = {"S_mean": ZERO, "S_sin": bt(150, 0), "S_cos": bt(0, 100)}
    with pytest.raises(LoadingScopeError, match="proportional"):
        model.evaluate_case(
            S0=None, E0=None, S1=None, E1=None, params=params, harmonics=out_of_phase
        )
    tension_mean = {"S_mean": bt(100, 0), "S_sin": bt(150, 0), "S_cos": ZERO}
    with pytest.raises(LoadingScopeError, match="Sig0"):
        model.evaluate_case(
            S0=None, E0=None, S1=None, E1=None, params=params, harmonics=tension_mean
        )
    out = model.evaluate_case(
        S0=None, E0=None, S1=None, E1=None, params={**params, "Sig0": 700.0}, harmonics=tension_mean
    )
    assert out.values[""] > kakuno_index(ZERO, bt(150, 0), ZERO, 420.0, 260.0)
    in_phase = {"S_mean": ZERO, "S_sin": bt(150, 100) * 0.6, "S_cos": bt(150, 100) * 0.8}
    out = model.evaluate_case(S0=None, E0=None, S1=None, E1=None, params=params, harmonics=in_phase)
    assert out.values[""] == pytest.approx(kakuno_index(ZERO, bt(150, 100), ZERO, 420.0, 260.0))
