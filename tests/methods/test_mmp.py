"""MMP (G05): Papuga & Fojtik (2017), Int J Fatigue 104:27-42, doi:10.1016/j.ijfatigue.2017.07.001,
checked against the PDF supplied by the maintainer.

Eqs. (19)-(25); Table 10 (material parameters, p.40) and Table 11 (Delta FI of MMP_WP, p.40) from
the local stresses of Tables 7-8 (outer and inner surface, pp.36-37). Table 11 prints the value of
the critical surface; plane-bending cases use b_-1 and b_0 (Section 5.3).

Findings: all 24 load cases are reproduced within 0.1 point of Delta FI (inputs printed to 0.1 MPa);
FF061 is critical on the inner surface, as the paper marks.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.invariant_methods.mmp import (
    MMP,
    mmp_beta,
    mmp_constants,
    mmp_index,
    signed_amplitudes,
)

TABLE = load_table("papuga_fojtik2017_tables7_11")
P = dict(f=244.7, f0=385.8, t=169.9, t0=330.8, b=315.9, b0=498.0)  # Table 10


def _tensor(a, t, r, s):
    return np.array([[a, s, 0.0], [s, t, 0.0], [0.0, 0.0, r]])


def _dFI(r):
    d = np.radians(r.phi_at)
    S_mean = _tensor(r.sigma_am, r.sigma_tm, r.sigma_rm, r.tau_m)
    S_sin = _tensor(r.sigma_aa, r.sigma_ta, r.sigma_ra, r.tau_a * np.cos(d))
    S_cos = _tensor(0.0, 0.0, 0.0, r.tau_a * np.sin(d))
    f, f0 = (P["b"], P["b0"]) if r.bending else (P["f"], P["f0"])
    return 100.0 * (mmp_index(S_mean, S_sin, S_cos, f, P["t"], f0, P["t0"]) - 1.0)


def test_table10_constants():
    """kappa and w to the printed digits; X_m within 0.002, its sensitivity to the 0.1 MPa rounding
    of t_-1 and t_0 (d X_m / d t_0 = -0.027 per MPa)."""
    kappa, w, X_m = mmp_constants(P["f"], P["t"], P["f0"], P["t0"])
    assert (kappa, w) == pytest.approx((1.440, 0.657), abs=6e-4)
    assert X_m == pytest.approx(0.233, abs=2e-3)
    kappa_b, w_b, X_mb = mmp_constants(P["b"], P["t"], P["b0"], P["t0"])
    assert (kappa_b, w_b) == pytest.approx((1.859, 0.657), abs=6e-4)
    assert X_mb == pytest.approx(0.301, abs=2e-3)


@pytest.mark.parametrize("mark", list(dict.fromkeys(TABLE["mark"])))
def test_table11_mmp_wp(mark):
    rows = TABLE[TABLE["mark"] == mark]
    printed = rows["dFI_MMP_WP"].iloc[0]
    assert max(_dFI(r) for r in rows.itertuples()) == pytest.approx(printed, abs=0.1)


@pytest.mark.parametrize(
    "S_mean, S_amp, expected",
    [
        (np.zeros((3, 3)), np.diag([244.7, 0, 0]), 1.0),  # fully reversed tension
        (np.zeros((3, 3)), _tensor(0, 0, 0, 169.9), 1.0),  # fully reversed torsion
        (np.diag([192.9, 0, 0]), np.diag([192.9, 0, 0]), 1.0),  # repeated tension, Eq. (21)
        (_tensor(0, 0, 0, 165.4), _tensor(0, 0, 0, 165.4), 1.0),  # repeated torsion, Eq. (25)
    ],
)
def test_calibration_limits(S_mean, S_amp, expected):
    index = mmp_index(S_mean, S_amp, np.zeros((3, 3)), P["f"], P["t"], P["f0"], P["t0"])
    assert index == pytest.approx(expected, rel=1e-12)


def test_eq22_beta():
    assert mmp_beta(1.0, 0.0) == 1.0  # repeated tension
    assert mmp_beta(2.0, -2.0) == 0.5  # repeated torsion (Table 9)
    assert mmp_beta(0.0, -1.0) == -1.0  # repeated compression: beneficial mean
    assert mmp_beta(0.0, 0.0) == 0.0


def test_phase_blind_and_signed_amplitudes():
    """The criterion ignores the phase (paper, Section 7); signs come from the in-phase part."""
    S_sin = _tensor(150.0, -40.0, 0.0, 0.0)
    assert signed_amplitudes(S_sin, np.zeros((3, 3)))[1, 1] == -40.0
    in_phase = mmp_index(
        np.zeros((3, 3)), _tensor(150, 0, 0, 90), np.zeros((3, 3)), 244.7, 169.9, 385.8, 330.8
    )
    out_of_phase = mmp_index(
        np.zeros((3, 3)), _tensor(150, 0, 0, 0), _tensor(0, 0, 0, 90), 244.7, 169.9, 385.8, 330.8
    )
    assert out_of_phase == pytest.approx(in_phase, rel=1e-12)


def test_missing_repeated_limits_are_estimated():
    """Sig0 from Su (Boehme et al. 2026 Eq. 1) and Tau0 from Eq. (23), as for PCR and MMMH."""
    S = np.diag([200.0, 0, 0])
    params = {"Sigm1": 244.7, "Taum1": 169.9, "Su": 560.0}
    f0 = 4 * 244.7 * 560.0 / (244.7 + 2 * 560.0)
    t0 = 4 * 169.9 / (1 + 2 * 244.7 / f0)
    expected = mmp_index(np.zeros((3, 3)), S, np.zeros((3, 3)), 244.7, 169.9, f0, t0)
    got = MMP().evaluate_case(S0=S, S1=-S, E0=S * 0, E1=S * 0, params=params).values[""]
    assert got == pytest.approx(expected, rel=1e-12)
