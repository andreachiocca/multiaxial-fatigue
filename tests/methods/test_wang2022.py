"""WANG2022 (N05): Wang, Sun, Liu, Liu & Shang (2022), Int J Damage Mech 31:1165-1186,
doi:10.1177/10567895221091308, checked against the scanned PDF supplied by the maintainer.

Eqs. (10), (15), (18)-(24); Tables 3-7, pp.1179-1181: lives predicted by the proposed model (column
"Pro."), from t_eq,a on the maximum-C_a plane set equal to the torsional S-N curve (Eq. 24):
Basquin tau'_f (2N)^b1 for LY12CZ, 30CrMnSiA and 6082-T6, the three-parameter form
t_-1 / (1 - kappa N^b2) for SM45C. The fatigue-limit validation (Fig. 3) cannot be reproduced:
Table 1 gives no yield strength for its nine materials (and lists "mild steel" with the limits of
42CrMo4).

Findings:
- every predicted life of Tables 4 (LY12CZ, 15 tests) and 5 (30CrMnSiA, 7) is reproduced within
  1 %, with the longest-chord C_a (rss misses the 90 deg rows by factors 1.4-3.2);
- with Eq. (18) as printed a tensile mean normal stress lowers t_eq,a; a '+' sign reproduces none of
  the tables. Mean-stress loading is out of scope (maintainer decision), so the SM45C rows 10-18
  raise LoadingScopeError;
- strict xfails (unresolved): 6082-T6 rows 6-11 (phases 94 and 125 deg) and SM45C rows 3, 5-9 at
  90 deg; neither the surface-plane restriction nor the maximum t_eq,a over planes reproduces all
  of them.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.cp_methods.wang2022 import Wang2022, wang_k
from multifat.models.errors import LoadingScopeError
from multifat.planes.sampling import build_R_list

TABLE = load_table("wang2022_tables4_7")
R = build_R_list()
Z = np.zeros((3, 3))
MATERIALS = {  # Table 3: f_-1, t_-1, sigma_y and the torsional S-N constants
    4: dict(f=168.73, t=119.62, sy=400.0, tf=602.8, b1=-0.1115),
    5: dict(f=508.0, t=293.283, sy=1105.0, tf=1766.713, b1=-0.12377),
    6: dict(f=126.0, t=74.0, sy=365.0, tf=359.9, b1=-0.0959),
    7: dict(f=442.0, t=311.0, sy=638.0, kappa=62.3, b2=-0.53),
}
UNRESOLVED = {
    **{(6, i): "6082-T6 at 94/125 deg: predicted life not reproduced" for i in range(6, 12)},
    **{(7, i): "SM45C at 90 deg: predicted life not reproduced" for i in (3, 5, 6, 7, 8, 9)},
}


def _case(r):
    d = np.radians(r.delta)

    def t(x, y):
        return np.array([[x, y, 0.0], [y, 0.0, 0.0], [0.0, 0.0, 0.0]])

    S_mean, S_sin, S_cos = (
        t(r.sigma_m, 0.0),
        t(r.sigma_a, r.tau_a * np.cos(d)),
        t(0.0, -r.tau_a * np.sin(d)),
    )
    return dict(
        S0=S_mean + S_sin, S1=S_mean - S_sin, E0=Z, E1=Z, R_list=R,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos, "E_sin": Z, "E_cos": Z},
    )  # fmt: skip


def _life(r):
    m = MATERIALS[r.table]
    params = {"Sigm1": m["f"], "Taum1": m["t"], "Sy": m["sy"]}
    t_eq = Wang2022().evaluate_case(params=params, **_case(r)).values[""]
    if r.table == 7:  # Eq. (23): t_eq = t_-1 / (1 - kappa N^b2)
        return ((1.0 - m["t"] / t_eq) / m["kappa"]) ** (1.0 / m["b2"])
    return 0.5 * (t_eq / m["tf"]) ** (1.0 / m["b1"])  # Eq. (22) with N = 2N_f / 2


def _rows():
    out = []
    for r in TABLE.itertuples():
        if r.sigma_m != 0:
            continue
        key = (r.table, r.row)
        mark = [pytest.mark.xfail(strict=True, reason=UNRESOLVED[key])] if key in UNRESOLVED else []
        out.append(pytest.param(r, marks=mark, id=f"T{r.table}-{r.row}"))
    return out


@pytest.mark.parametrize("r", _rows())
def test_tables_4_7_predicted_lives(r):
    assert _life(r) == pytest.approx(r.N_pro, rel=0.01)


@pytest.mark.parametrize(
    "r", [r for r in TABLE.itertuples() if r.sigma_m != 0], ids=lambda r: f"T7-{r.row}"
)
def test_mean_stress_rows_are_out_of_scope(r):
    m = MATERIALS[7]
    with pytest.raises(LoadingScopeError, match="mean normal stress"):
        Wang2022().evaluate_case(
            params={"Sigm1": m["f"], "Taum1": m["t"], "Sy": m["sy"]}, **_case(r)
        )


@pytest.mark.parametrize("table", [4, 5, 6, 7])
def test_eq19_meets_both_reversed_limits(table):
    """Torsion: C_a = t_-1, N_a = 0. Bending: C_a = N_a = f_-1/2 on the 45 deg plane."""
    m = MATERIALS[table]
    k = wang_k(m["f"], m["t"], m["sy"])
    bending = np.sqrt(
        (m["f"] / 2) ** 2 + k * (1 - (m["f"] / 2) / (2 * m["sy"])) ** 2 * (m["f"] / 2) ** 2
    )
    assert bending == pytest.approx(m["t"], rel=1e-12)
    for sigma_a, tau_a in ((m["f"], 0.0), (0.0, m["t"])):
        r = type("R", (), dict(table=table, sigma_a=sigma_a, sigma_m=0.0, tau_a=tau_a, delta=0.0))
        params = {"Sigm1": m["f"], "Taum1": m["t"], "Sy": m["sy"]}
        assert Wang2022().evaluate_case(params=params, **_case(r)).values[""] == pytest.approx(
            m["t"], rel=1e-9
        )


def test_printed_eq18_makes_a_tensile_mean_beneficial():
    """The finding behind the scope rule, on the plane-level formula."""
    from multifat.models.cp_methods.wang2022 import _WangPlane

    S = np.zeros((3, 3))
    S[0, 2] = S[2, 0] = 100.0
    S[2, 2] = 150.0
    params = {"Sigm1": 442.0, "Taum1": 311.0, "Sy": 638.0}
    reversed_ = _WangPlane().evaluate_on_plane(S0r=S, S1r=-S, E0r=Z, E1r=Z, params=params).damage
    shifted = np.diag([0.0, 0.0, 100.0])
    with_mean = (
        _WangPlane()
        .evaluate_on_plane(S0r=S + shifted, S1r=-S + shifted, E0r=Z, E1r=Z, params=params)
        .damage
    )
    assert with_mean < reversed_
