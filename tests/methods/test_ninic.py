"""NINIC (N11): Ninic (2006), Int J Fatigue 28:103-113, doi:10.1016/j.ijfatigue.2005.04.014,
checked against the PDF supplied by the maintainer.

Eqs. (1)-(4), p.104; Table 1, p.106 (T_e, S_e, k); Tables 2-5, pp.106-111: bending
sigma_y = sigma_a sin(wt) and torsion tau_xy = tau_a sin(wt - phi) (Eqs. 15-16), with the error
index I = 100 (D - 1), Eq. (14), of the Damage Function on its critical plane and the plane
stresses tau_a and sigma_n,max.

Findings:
- Eq. (4) reproduces k of Table 1 for the four materials;
- 29 of the 32 full-search rows of Tables 2-4 are reproduced within 0.08 points by the
  three-dimensional search with the rss shear amplitude (the longest chord misses rows 2-9 and 2-10 by 1.7-1.9 points);
- strict xfails: rows 2-7 and 3-12, where the printed plane is not the maximum of D (the printed
  plane stresses give the printed I; the search finds a larger D, by 0.24 and 1.8 points); row
  3-4, whose printed tau_xy,a = 21.9 MPa contradicts its printed ratio 0.2 (the paper's Table 1 of
  Goncalves et al. 2005 gives 62.8 MPa for this test); row 5-6, whose printed sigma_n,max = 117.7
  MPa exceeds any principal stress of the printed load (ratio 0.31 also contradicts 81.6/67.6);
- Table 5 (aluminium 6082-T6, rows 1-5) agrees within 0.6 points only: the printed plane stresses
  imply slightly different loads (row 1 prints sigma_n,max = 16 for sigma_a = 15 MPa);
- rows 2-11 and 3-13 restrict the plane to the surface normal; 2-11 is reproduced within 0.04
  points of the text value I = 0.0735 % (p.108; the table prints 0.735), 3-13 is 1.8 points off,
  as row 3-12 of the same load;
- the cast-iron uniaxial limit (row 4-5) gives D = sqrt(k): Eq. (4) calibrates the uniaxial limit
  only for S_e/T_e >= sqrt(2).
"""

import numpy as np
import pytest
from _paper import frame_with_normal, load_table

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.ninic import Ninic, ninic_k
from multifat.planes.sampling import build_R_list

TABLE = load_table("ninic2006_tables2_5")
R_ALL = build_R_list()
# planes normal to the surface: normal in the x-y plane (z = surface normal, y = bar axis, Fig. 1)
R_SURFACE = np.array(
    [frame_with_normal([np.cos(t), np.sin(t), 0.0]) for t in np.radians(np.arange(0, 180, 0.05))]
)
Z = np.zeros((3, 3))
XFAIL = {
    "2-7": "printed critical plane is not the maximum of D (search finds +0.24 points)",
    "3-12": "printed critical plane is not the maximum of D (search finds +1.8 points)",
    "3-4": "printed tau_xy,a = 21.9 MPa contradicts the printed ratio 0.2",
    "5-6": "printed sigma_n,max exceeds every principal stress of the printed load",
}
SURFACE_ROWS = {"2-11", "3-13"}


def _load(sigma_a, tau_a, phi):
    p = np.radians(phi)
    bending_torsion = lambda s, x: np.array(
        [[0.0, x, 0.0], [x, s, 0.0], [0.0, 0.0, 0.0]]
    )  # noqa: E731
    return bending_torsion(sigma_a, tau_a * np.cos(p)), bending_torsion(0.0, -tau_a * np.sin(p))


def _index(r, R_list):
    S_sin, S_cos = _load(r.sigma_a, r.tau_a, r.phi)
    best = best_cp_for_case(
        S0=S_sin,
        S1=-S_sin,
        E0=Z,
        E1=Z,
        R_list=R_list,
        refine=R_list is R_ALL,
        method=Ninic(),
        params={"Sigm1": r.Se, "Taum1": r.Te, "k_NIN": r.k},
        S_sin=S_sin,
        S_cos=S_cos,
        E_sin=Z,
        E_cos=Z,
    )
    return 100 * (best.by_metric - 1)


def _rows():
    out = []
    for r in TABLE.itertuples():
        if r.row in SURFACE_ROWS:
            continue
        if r.row in XFAIL:
            out.append(pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=XFAIL[r.row])))
        else:
            out.append(r)
    return out


@pytest.mark.parametrize(
    "material,Te,Se,k",
    [
        ("mild steel", 137.3, 235.4, 0.7793),
        ("hard steel", 196.2, 313.9, 0.9217),
        ("cast iron", 91.20, 96.10, 0.8021),
        ("Al 6082-T6", 87.0, 145.0, 0.8488),
    ],
)
def test_table1_eq4_k(material, Te, Se, k):
    assert ninic_k(Se, Te) == pytest.approx(k, abs=5e-5)


@pytest.mark.parametrize("r", _rows(), ids=lambda r: r.row)
def test_tables_2_5_damage_function(r):
    tol = 0.6 if r.row.startswith("5-") else 0.08  # Table 5: rounded loads, see module docstring
    assert _index(r, R_ALL) == pytest.approx(r.I, abs=tol)


def test_surface_plane_reevaluation_row_2_11():
    """Row 2-11 (row 2-10 on planes normal to the surface): text value I = 0.0735 %, p.108."""
    r = TABLE.set_index("row").loc["2-11"]
    assert _index(r, R_SURFACE) == pytest.approx(0.0735, abs=0.04)


@pytest.mark.xfail(strict=True, reason="row 3-13 is 1.8 points off, as row 3-12 of the same load")
def test_surface_plane_reevaluation_row_3_13():
    r = TABLE.set_index("row").loc["3-13"]
    assert _index(r, R_SURFACE) == pytest.approx(r.I, abs=0.08)


@pytest.mark.parametrize("r", list(TABLE.itertuples()), ids=lambda r: r.row)
def test_printed_plane_stresses_give_the_printed_index(r):
    """Eq. (13) on the printed (tau_a, sigma_n,max); rows 2-11 and 4-8 carry table misprints."""
    if r.row in ("2-11", "4-8"):
        pytest.skip("2-11 prints 0.735 for 0.0735 (text p.108); 4-8 prints tau_a 91.6 for 81.6")
    D = np.sqrt((r.tau_a_plane / r.Te) ** 2 + r.k * (r.sigma_n_max / r.Se) ** 2)
    assert 100 * (D - 1) == pytest.approx(r.I, abs=0.06)


def test_uniaxial_limit_needs_se_over_te_at_least_sqrt2():
    """Eq. (4) gives D = 1 at the uniaxial limit for S_e/T_e >= sqrt(2), D = sqrt(k) below."""
    for Te, Se in ((137.3, 235.4), (91.2, 96.1)):
        S_sin, S_cos = _load(Se, 0.0, 0.0)
        best = best_cp_for_case(
            S0=S_sin, S1=-S_sin, E0=Z, E1=Z, R_list=R_ALL, method=Ninic(),
            params={"Sigm1": Se, "Taum1": Te}, S_sin=S_sin, S_cos=S_cos, E_sin=Z, E_cos=Z,
        )  # fmt: skip
        expected = 1.0 if Se / Te >= np.sqrt(2) else np.sqrt(ninic_k(Se, Te))
        assert best.by_metric == pytest.approx(expected, abs=1e-6)
    with pytest.warns(RuntimeWarning, match="sqrt\\(2\\)"):
        Ninic().evaluate_on_plane(S0r=Z, S1r=Z, E0r=Z, E1r=Z, params={"Sigm1": 96.1, "Taum1": 91.2})


def test_torsion_limit_gives_unit_damage():
    for Te, Se in ((137.3, 235.4), (91.2, 96.1)):
        S_sin, S_cos = _load(0.0, Te, 0.0)
        best = best_cp_for_case(
            S0=S_sin, S1=-S_sin, E0=Z, E1=Z, R_list=R_ALL, method=Ninic(),
            params={"Sigm1": Se, "Taum1": Te, "k_NIN": ninic_k(Se, Te)},
            S_sin=S_sin, S_cos=S_cos, E_sin=Z, E_cos=Z,
        )  # fmt: skip
        assert best.by_metric == pytest.approx(1.0, abs=1e-9)


def test_compressive_normal_stress_adds_no_damage():
    S = np.zeros((3, 3))
    S[0, 2] = S[2, 0] = 50.0
    S[2, 2] = -100.0
    S1 = S.copy()
    S1[0, 2] = S1[2, 0] = -50.0  # shear amplitude 50 MPa, constant compression on the plane
    params = {"Sigm1": 235.4, "Taum1": 137.3}
    r = Ninic().evaluate_on_plane(S0r=S, S1r=S1, E0r=Z, E1r=Z, params=params)
    assert r.damage == pytest.approx(50.0 / 137.3, rel=1e-12)
