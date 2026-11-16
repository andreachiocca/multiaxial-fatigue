"""Wang, Shang, Sun & Chen (2018), doi:10.1016/j.tafmec.2017.10.007: Eqs. (12)-(20), Section 5.1."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from multifat.models.cp_methods.wang2018 import (
    Wang2018,
    tension_torsion_parts,
    wang2018_constants,
    wang2018_tau_eq,
    weighted_plane_angle,
)
from multifat.models.errors import LoadingScopeError, MaterialScopeError

DATA = Path(__file__).resolve().parents[1] / "data" / "papers" / "wang2018_tables3_5.csv"
# Table 2 S-N curves at N_0 = 1e6 (Section 5.2)
F_M1 = 888.90 * 1e6**-0.10042
T_M1 = 688.85 * 1e6**-0.12233


def _rows():
    with DATA.open(encoding="utf-8") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def _parts(row):
    """(sigma, tau) as (mean, sin, cos) for tau = tau_a sin(wt - delta) + tau_m (Eq. 8)."""
    d = np.radians(float(row["delta"]))
    sigma = (float(row["sigma_m"]), float(row["sigma_a"]), 0.0)
    ta = float(row["tau_a"])
    return sigma, (float(row["tau_m"]), ta * np.cos(d), -ta * np.sin(d))


def test_table_data_complete():
    rows = _rows()
    assert [int(r["specimen"]) for r in rows] == list(range(1, 53))


@pytest.mark.parametrize("f_0, t_0", [(300.0, 200.0), (None, None)])
def test_eqs_12_to_14_limits_give_t_m1(f_0, t_0):
    """Eq. (15) makes the three calibration limits of Eqs. (12)-(14) satisfy Eq. (16) exactly."""
    f_m1, t_m1 = 222.0, 127.0
    f0 = np.sqrt(2.0) * f_m1 if f_0 is None else f_0
    t0 = 4.0 * t_m1 / (2.0 * f_m1 / f0 + 1.0) if t_0 is None else t_0
    cases = {
        "reversed torsion": ((0, 0, 0), (0, t_m1, 0)),
        "pulsating tension": ((f0 / 2, f0 / 2, 0), (0, 0, 0)),
        "pulsating torsion": ((0, 0, 0), (t0 / 2, t0 / 2, 0)),
    }
    for name, (sigma, tau) in cases.items():
        value, _ = wang2018_tau_eq(sigma, tau, f_m1, t_m1, f_0, t_0)
        assert value == pytest.approx(t_m1, rel=1e-12), name


def test_eq_20_is_eq_15_with_eq_19():
    k_tau, k_sigma = wang2018_constants(F_M1, T_M1)
    assert k_tau == pytest.approx((np.sqrt(2.0) - 1.0) / 2.0, rel=1e-14)
    expected = np.sqrt(24 * T_M1**2 / F_M1**2 - (9 + 6 * np.sqrt(2.0)) / 4) - 1
    assert k_sigma == pytest.approx(expected, rel=1e-14)


def test_section_5_1_estimated_crack_plane_70_deg():
    """Constant amplitude, 30 and 90 deg out of phase: estimated crack orientation 70 deg.

    The plane normal lies at |phi_c| from the axis, so the crack trace is at 90 - |phi_c|.
    The ratios tau_a/sigma_a and tau_m/sigma_m are the same in every Table 5 test, so the
    estimate does not depend on the specimen.
    """
    for row in _rows():
        if row["table"] == "5" and row["delta"] in ("30", "90"):
            _, phi_c = wang2018_tau_eq(*_parts(row), F_M1, T_M1)
            assert abs(90.0 - abs(np.degrees(phi_c)) - 70.0) < 0.5, row["specimen"]


def test_eqs_8_9_constant_tau_max_uses_largest_normal_stress_plane():
    sigma, tau = (0.0, 100.0, 0.0), (0.0, 0.0, -50.0)  # 90 deg, tau_a / sigma_a = 0.5
    assert weighted_plane_angle(sigma, tau) is None
    value, phi_c = wang2018_tau_eq(sigma, tau, F_M1, T_M1)
    assert phi_c == pytest.approx(0.0, abs=1e-12)
    assert value == pytest.approx(np.sqrt(50.0**2 + 100.0**2 / 3.0), rel=1e-12)


def test_reversed_tension_and_torsion_closed_forms():
    assert wang2018_tau_eq((0, 100, 0), (0, 0, 0), F_M1, T_M1)[0] == pytest.approx(100 / 3**0.5)
    assert wang2018_tau_eq((0, 0, 0), (0, 100, 0), F_M1, T_M1)[0] == pytest.approx(100.0)


def test_eq_17_violation_is_out_of_material_scope():
    with pytest.raises(MaterialScopeError):
        wang2018_constants(222.0, 127.0, 300.0, 300.0)  # k_tau = 2 * 127 / 300 - 1 < 0


def test_only_tension_torsion_is_in_scope():
    S = np.zeros((3, 3))
    S[0, 0] = 50.0  # hoop stress
    with pytest.raises(LoadingScopeError):
        tension_torsion_parts(S, np.zeros((3, 3)), np.zeros((3, 3)))


def test_model_uses_axial_y_and_shear_xy():
    S_sin = np.zeros((3, 3))
    S_sin[1, 1] = 144.57
    S_cos = np.zeros((3, 3))
    S_cos[0, 1] = S_cos[1, 0] = -83.47
    S_mean = np.zeros((3, 3))
    S_mean[1, 1], S_mean[0, 1], S_mean[1, 0] = 163.02, 94.12, 94.12
    harmonics = {"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos}
    model = Wang2018()
    assert model.design_loading == "torsion"
    params = {"Sigm1": F_M1, "Taum1": T_M1}
    value = model.evaluate_case(
        S0=None, E0=None, S1=None, E1=None, params=params, harmonics=harmonics
    ).values[""]
    expected = wang2018_tau_eq((163.02, 144.57, 0.0), (94.12, 0.0, -83.47), F_M1, T_M1)[0]
    assert value == pytest.approx(expected, rel=1e-12)
