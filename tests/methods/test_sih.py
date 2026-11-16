"""SIH, Zenner, Simbuerger & Liu (2000), doi:10.1016/S0142-1123(99)00107-3: Eqs. (5)-(7), Table 1."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.models.errors import LoadingScopeError, MaterialScopeError
from multifat.models.integral_methods.bohme_papuga import _plane_components, _sphere_grid
from multifat.models.integral_methods.sih import SIH, SPHERE, sih_constants, sih_index

ZERO = np.zeros((3, 3))
F_W, T_W, F_SCH, T_SCH = 300.0, 190.0, 480.0, 330.0  # a steel-like set, tau_W/sigma_W = 0.633


def axial(s):
    S = np.zeros((3, 3))
    S[1, 1] = s
    return S


def shear(t):
    S = np.zeros((3, 3))
    S[0, 1] = S[1, 0] = t
    return S


@pytest.mark.parametrize(
    "S, moments",
    [
        (axial(1.0), {"t2": 1.0, "s2": 1.5, "t4": 4 / 21, "s3": 15 / 14}),
        (shear(1.0), {"t2": 3.0, "s2": 2.0, "t4": 12 / 7, "s3": 0.0}),
    ],
    ids=["tension", "torsion"],
)
def test_table_1_factors_are_sphere_integrals(S, moments):
    """15/(8 pi) times the sphere integrals behind a, b and the 12/7, 4/21, 15/14 of Table 1."""
    normals, weights = _sphere_grid(48, 96)
    sigma, tau_vec = _plane_components(S, normals)
    tau = np.linalg.norm(tau_vec, axis=1)
    values = {"t2": tau**2, "s2": sigma**2, "t4": tau**4, "s3": sigma**3}
    for key, expected in moments.items():
        assert SPHERE * np.sum(weights * values[key]) == pytest.approx(expected, abs=1e-12), key


@pytest.mark.parametrize("rule", ["longest_chord", "rss"])
@pytest.mark.parametrize("t_0", [T_SCH, None])
def test_four_calibration_limits_are_exact(rule, t_0):
    tau_sch = T_SCH if t_0 is not None else 4 * T_W / (2 * F_W / F_SCH + 1)  # Eq. (7)
    cases = {
        "reversed tension": (ZERO, axial(F_W)),
        "reversed torsion": (ZERO, shear(T_W)),
        "pulsating tension": (axial(F_SCH / 2), axial(F_SCH / 2)),
        "pulsating torsion": (shear(tau_sch / 2), shear(tau_sch / 2)),
    }
    for name, (S_mean, S_sin) in cases.items():
        value = sih_index(S_mean, S_sin, ZERO, F_W, T_W, F_SCH, t_0, rule)
        assert value == pytest.approx(1.0, rel=1e-10), name


def test_eq_5_is_the_case_a_1_b_0_with_ratio_1_over_sqrt3():
    """Eq. (5) yields tau_W/sigma_W = 1/sqrt(3): Table 1 then gives a = 1, b = 0."""
    k = sih_constants(1.0, 1.0 / np.sqrt(3.0))
    assert k["a"] == pytest.approx(1.0) and k["b"] == pytest.approx(0.0, abs=1e-15)


def test_eq_5b_constant_principal_directions_give_octahedral_shear():
    """Eq. (5b): sigma_eq,a = 3/sqrt(2) tau_oct,a for proportional loading (a = 1, b = 0)."""
    S = np.array([[40.0, 25.0, 0.0], [25.0, 120.0, 10.0], [0.0, 10.0, -30.0]])
    dev = S - np.trace(S) / 3 * np.eye(3)
    tau_oct = np.sqrt(np.sum(dev * dev) / 3.0)
    value = sih_index(ZERO, S, ZERO, 1.0, 1.0 / np.sqrt(3.0))  # returns sigma_eq,a / 1
    assert value == pytest.approx(3.0 / np.sqrt(2.0) * tau_oct, rel=1e-10)


def test_section_3_phase_independence_holds_for_the_quadratic_amplitude():
    """Section 3: for tau_xy,a = sigma_x,a / 2 the integral of Eq. (5) does not depend on the
    phase. This holds exactly with the quadratic (rss) amplitude; the longest chord (the
    maintainer's confirmed default) lowers it at 90 deg."""
    out = {}
    for rule in ("rss", "longest_chord"):
        for d in (0, 45, 90):
            r = np.radians(d)
            S_sin, S_cos = axial(1.0) + shear(0.5 * np.cos(r)), shear(0.5 * np.sin(r))
            out[rule, d] = sih_index(ZERO, S_sin, S_cos, 1.0, 1.0 / np.sqrt(3.0), shear_rule=rule)
    assert out["rss", 45] == pytest.approx(out["rss", 0], rel=1e-12)
    assert out["rss", 90] == pytest.approx(out["rss", 0], rel=1e-12)
    assert out["longest_chord", 0] == pytest.approx(out["rss", 0], rel=1e-12)
    assert out["longest_chord", 90] / out["longest_chord", 0] == pytest.approx(0.9294, abs=1e-4)


def test_scope_errors():
    with pytest.raises(MaterialScopeError):
        sih_constants(300.0, 280.0)  # tau_W / sigma_W = 0.93 > 0.866
    with pytest.raises(LoadingScopeError):
        sih_index(axial(50.0), axial(100.0), ZERO, F_W, T_W)  # mean without Sig0
    assert sih_index(ZERO, axial(F_W), ZERO, F_W, T_W) == pytest.approx(1.0)


def test_model_reads_its_shear_rule():
    model = SIH()
    assert model.np_shear_rule == "longest_chord"
    harmonics = {"S_mean": ZERO, "S_sin": axial(200.0), "S_cos": shear(100.0)}
    params = {"Sigm1": F_W, "Taum1": T_W}
    values = {}
    for rule in ("longest_chord", "rss"):
        model.np_shear_rule = rule
        values[rule] = model.evaluate_case(
            S0=None, E0=None, S1=None, E1=None, params=params, harmonics=harmonics
        ).values[""]
        assert values[rule] == pytest.approx(
            sih_index(ZERO, axial(200.0), shear(100.0), F_W, T_W, shear_rule=rule)
        )
    assert values["longest_chord"] < values["rss"]
