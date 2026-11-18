"""LIUMAHADEVAN (N16): Liu & Mahadevan (2005), Int J Fatigue 27:790-800,
doi:10.1016/j.ijfatigue.2005.01.003, checked against the PDF supplied by the maintainer.

Eqs. (6)-(24), Table 1, Figs. 1, 2-5 and 9. Figs. 2-5 plot the mean absolute error index (Eq. 25)
per phase on the Papadopoulos et al. (1997) Table 1 data (tests/data/papers/papadopoulos1997_table1
.csv, phases 120 and 180 deg left out, as in the figures: their Matake and McDiarmid curves are the
phase means of that table's indices). The "present" values below are read from the figures (+-0.2).

Findings:
- Eq. (19) has removable 0/0 points at s = 1/2 and s = 1 (rationalised form used);
- with the worst plane on the alpha-cone and the rss shear amplitude, Fig. 2 (hard steel) is
  reproduced at 0, 60 and 90 deg within 0.5 point (the longest chord gives 4.6 at 90 deg);
- hard steel 30 deg, 42CrMo4, 34Cr4 and 30NCD16 are not reproduced (balanced 90 deg out-of-phase
  rows come out about 10 % non-conservative here, small in the paper); unresolved, strict xfails.
  The life validation (Figs. 10-13) needs data not given in the paper.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.cp_methods.liu_mahadevan import (
    LiuMahadevan,
    fracture_normals,
    liu_mahadevan_index,
    liu_mahadevan_life,
    material_angle,
    mean_stress_factor,
)

Z = np.zeros((3, 3))
TABLE = load_table("papadopoulos1997_table1")
TABLE = TABLE[~TABLE["delta"].isin([120, 180])]
FIG = {  # mean absolute error (%) of the "present" model, read from Figs. 2-5
    "hard steel": {0: 1.9, 30: 3.6, 60: 1.7, 90: 3.0},
    "42CrMo4": {0: 8.3, 90: 8.9},
    "34Cr4": {0: 5.3, 60: 4.3, 90: 5.3},
    "30NCD16": {0: 9.0, 45: 7.1, 60: 6.5, 90: 3.2},
}
LM_LIMITS = {  # Table 2 (f_-1, t_-1), identical to Papadopoulos et al. (1997)
    "hard steel": (313.9, 196.2), "42CrMo4": (398.0, 260.0),
    "34Cr4": (410.0, 256.0), "30NCD16": (660.0, 410.0),
}  # fmt: skip
REPRODUCED = {("hard steel", 0), ("hard steel", 60), ("hard steel", 90), ("30NCD16", 60)}


def _tensors(sigma_a, sigma_m, tau_a, tau_m, delta):
    d = np.radians(delta)
    S_mean = np.array([[sigma_m, tau_m, 0], [tau_m, 0, 0], [0, 0, 0]], float)
    S_sin = np.array([[sigma_a, tau_a * np.cos(d), 0], [tau_a * np.cos(d), 0, 0], [0, 0, 0]])
    S_cos = np.array([[0, -tau_a * np.sin(d), 0], [-tau_a * np.sin(d), 0, 0], [0, 0, 0]])
    return S_mean, S_sin, S_cos


@pytest.mark.parametrize(
    "s, alpha, beta, k",
    [
        (0.5, np.pi / 3, 0.901, 0.0),  # removable 0/0 of Eq. (19)
        (1.0, 0.0, 1.0, 0.0),  # removable 0/0; alpha = 0, Table 1
        (1.28, 0.0, 1.28, 9 * (1.28**2 - 1)),  # Table 1, s > 1
    ],
)
def test_eq19_table1(s, alpha, beta, k):
    assert material_angle(s) == pytest.approx((alpha, beta, k), abs=2e-3)


@pytest.mark.parametrize("s", [0.5, 0.55, 0.58, 0.625, 0.8, 0.95, 1.0, 1.28])
def test_both_reversed_limits_on_the_alpha_cone(s):
    """Eq. (19) / Table 1: bending f_-1 and torsion t_-1 give exactly beta (index 1)."""
    f, t = 400.0, 400.0 * s
    assert liu_mahadevan_index(*_tensors(f, 0, 0, 0, 0), f, t) == pytest.approx(1.0, rel=1e-6)
    assert liu_mahadevan_index(*_tensors(0, 0, t, 0, 0), f, t) == pytest.approx(1.0, rel=1e-6)


def test_below_tresca_ratio_is_rejected():
    """The worst plane of the alpha-cone exceeds the torsion limit for s < 0.5 (Fig. 9 starts at
    s = 0.35 in the paper's 2D construction)."""
    with pytest.raises(ValueError, match="Tresca"):
        material_angle(0.4)


def test_eq21_eta():
    assert mean_stress_factor(1.0) == 1.0 and mean_stress_factor(1.3) == 1.0
    assert mean_stress_factor(0.6) == pytest.approx(0.772, abs=1e-3)  # Fig. 1 curve ~0.75


def test_fracture_plane():
    """Maximum normal stress amplitude: the axis in tension, +-45 deg in torsion."""
    (n,) = fracture_normals(np.diag([1.0, 0, 0]), Z)
    assert abs(n[0]) == pytest.approx(1.0, abs=1e-6)
    normals = fracture_normals(_tensors(0, 0, 1.0, 0, 0)[1], Z)
    assert len(normals) == 2
    assert all(abs(n[0]) == pytest.approx(np.sqrt(0.5), abs=1e-5) for n in normals)


def _groups():
    out = []
    for (mat, delta), rows in TABLE.groupby(["material", "delta"]):
        key = (mat, int(delta))
        mark = [] if key in REPRODUCED else [
            pytest.mark.xfail(strict=True, reason="Figs. 2-5 not reproduced (unresolved)")]  # fmt: skip
        out.append(pytest.param(mat, int(delta), rows, marks=mark, id=f"{mat}-{int(delta)}"))
    return out


@pytest.mark.parametrize("mat, delta, rows", _groups())
def test_figs2_5_mean_absolute_error(mat, delta, rows):
    f, t = LM_LIMITS[mat]
    errors = []
    for r in rows.itertuples():
        index = liu_mahadevan_index(*_tensors(r.sigma_a, r.sigma_m, r.tau_a, r.tau_m, delta), f, t)
        errors.append(100.0 * (1.0 - 1.0 / index))  # Eq. (25): (LHS - RHS) / LHS
    assert np.mean(np.abs(errors)) == pytest.approx(FIG[mat][delta], abs=0.5)


SN = dict(sigma_f=1200.0, b_fat=-0.08, tau_f=700.0, b0=-0.06)  # Basquin, (2N) convention


@pytest.mark.parametrize("N", [1e4, 1e5, 1e6])
def test_eq23_life_reproduces_both_curves(N):
    """Eq. (23)-(24): sigma_a = f_N in tension and tau_a = t_N in torsion give N_f = N."""
    f_n = SN["sigma_f"] * (2 * N) ** SN["b_fat"]
    t_n = SN["tau_f"] * (2 * N) ** SN["b0"]
    assert liu_mahadevan_life(*_tensors(f_n, 0, 0, 0, 0), **SN) == pytest.approx(N, rel=1e-6)
    assert liu_mahadevan_life(*_tensors(0, 0, t_n, 0, 0), **SN) == pytest.approx(N, rel=1e-6)


def test_case_level_model_and_native_life():
    model = LiuMahadevan()
    params = {"Sigm1": 400.0, "Taum1": 260.0, **SN}
    S_mean, S_sin, S_cos = _tensors(300.0, 50.0, 150.0, 0.0, 90.0)
    H = {"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos, "E_sin": Z, "E_cos": Z}
    res = model.evaluate_case(S0=S_mean + S_sin, S1=S_mean - S_sin, E0=Z, E1=Z, params=params,
                              harmonics=H)  # fmt: skip
    assert res.aux["Nf"] == model.native_life({"dp": res.values[""], **res.aux}, params)
    tensile = liu_mahadevan_life(S_mean, S_sin, S_cos, **SN)
    compressive = liu_mahadevan_life(-S_mean, S_sin, S_cos, **SN)
    assert compressive > res.aux["Nf"] == pytest.approx(tensile)  # tensile mean is damaging
