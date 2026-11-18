"""ZHANGYAO (N06): Zhang & Yao (2011), Fatigue Fract Engng Mater Struct 34:337-344,
doi:10.1111/j.1460-2695.2010.01523.x, checked against the PDF supplied by the maintainer.

Eqs. (3), (6), (7), (A1), pp.338-344; Tables 2-8, pp.340-343: 66 fatigue limits (five steels and
LY12CZ) with the error index I = 100 (LHS - t_-1)/t_-1, Eq. (10), of the proposed criterion.

Findings:
- Eq. (7) as printed (Goodman on the plane amplitude) reproduces the five 42CrMo4 mean-stress rows;
  N_max = N_a + N_m and other readings do not;
- the shear amplitude is the major semi-axis (longest chord) of Eq. (A1): 18/66 rows within 0.06
  and 50/66 within 0.3 points, against 6/66 with rss;
- the plane rule (maximum C_a) is ill-conditioned for this criterion: C_a is stationary at its
  maximum while N* varies to first order, so the index moves by up to 0.3 points for a 1 deg
  plane grid (in-phase rows), and over planes within 1e-4 of the maximum C_a it spans several
  points, or tens of points in 90 deg out-of-phase loading where the shear amplitude is equal
  on a continuum of planes. The repository resolves ties towards the larger index; the printed
  values mostly sit at the top of the tied range. Hence the 0.3-point tolerance;
- strict xfails: 90 deg rows whose maximum C_a is shared by a continuum of planes (printed value
  inside the tied range), and mean-stress rows (34Cr4, 30NCD16) not reproduced on any
  near-maximum plane with Eq. (7) and the sigma_u of Table 2 (unresolved; no single sigma_u fits).
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.zhang_yao import ZhangYao, zhang_yao_lhs
from multifat.models.errors import LoadingScopeError, MissingParamsError
from multifat.planes.sampling import build_R_list

TABLE = load_table("zhang_yao2011_tables3_8")
R = build_R_list()
Z = np.zeros((3, 3))
TIED = "maximum C_a shared by a continuum of planes; the printed value lies in the tied range"
UNRESOLVED = "mean-stress row not reproduced with Eq. (7) and Table 2 sigma_u (unresolved)"
XFAIL = {
    **dict.fromkeys(["3-8", "5-3", "5-11", "6-6", "7-6", "8-14", "8-15"], TIED),
    **dict.fromkeys(["5-7", "5-8", "5-10", "6-3", "6-4", "6-7", "6-8", "6-9", "6-10"], UNRESOLVED),
}


def _tension_torsion(r):
    def t(x, y):
        return np.array([[x, y, 0.0], [y, 0.0, 0.0], [0.0, 0.0, 0.0]])

    d = np.radians(r.delta)
    return t(r.sigma_m, r.tau_m), t(r.sigma_a, r.tau_a * np.cos(d)), t(0.0, -r.tau_a * np.sin(d))


def _index(r):
    S_mean, S_sin, S_cos = _tension_torsion(r)
    params = {"Sigm1": r.f_m1, "Taum1": r.t_m1}
    if r.su == r.su:  # Table 8 (LY12CZ) gives no sigma_u and has no mean stresses
        params["Su"] = r.su
    best = best_cp_for_case(
        S0=S_mean + S_sin, S1=S_mean - S_sin, E0=Z, E1=Z, R_list=R, method=ZhangYao(),
        params=params, S_sin=S_sin, S_cos=S_cos, E_sin=Z, E_cos=Z,
    )  # fmt: skip
    return 100 * (best.by_metric - 1)


def _rows():
    return [
        (
            pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=XFAIL[r.row]))
            if r.row in XFAIL
            else r
        )
        for r in TABLE.itertuples()
    ]


@pytest.mark.parametrize("r", _rows(), ids=lambda r: r.row)
def test_tables_3_8_error_index(r):
    assert _index(r) == pytest.approx(r.I_current, abs=0.3)


@pytest.mark.parametrize("f_m1,t_m1", [(313.9, 196.2), (398, 260), (235.4, 137.3)])
def test_eq6_meets_both_reversed_limits(f_m1, t_m1):
    """Torsion: C_a = t_-1, N = 0. Bending: C_a = N = f_-1/2 on the 45 deg plane."""
    assert zhang_yao_lhs(t_m1, 0.0, f_m1, t_m1) == pytest.approx(t_m1, rel=1e-12)
    assert zhang_yao_lhs(f_m1 / 2, f_m1 / 2, f_m1, t_m1) == pytest.approx(t_m1, rel=1e-12)


def test_eq7_mean_stress_and_scope():
    S = np.zeros((3, 3))
    S[0, 2] = S[2, 0] = 100.0
    S0, S1 = S + np.diag([0.0, 0.0, 150.0]), -S + np.diag([0.0, 0.0, 50.0])  # N_a 50, N_m 100
    params = {"Sigm1": 300.0, "Taum1": 200.0, "Su": 800.0}
    r = ZhangYao().evaluate_on_plane(S0r=S0, S1r=S1, E0r=Z, E1r=Z, params=params)
    N_star = 50.0 / (1 - 100.0 / 800.0)
    assert r.damage == pytest.approx(zhang_yao_lhs(100.0, N_star, 300.0, 200.0) / 200.0)
    with pytest.raises(MissingParamsError):
        ZhangYao().evaluate_on_plane(S0r=S0, S1r=S1, E0r=Z, E1r=Z, params={**params, "Su": None})
    with pytest.raises(LoadingScopeError):
        ZhangYao().evaluate_on_plane(S0r=S0, S1r=S1, E0r=Z, E1r=Z, params={**params, "Su": 90.0})
