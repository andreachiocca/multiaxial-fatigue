"""MATAKE and MCD against the Matake and McDiarmid error-index columns of Papadopoulos et al.
(1997) Int J Fatigue 19:219-235, doi:10.1016/S0142-1123(96)00064-3 (Zotero BW24SJ5K); Eqs. (49)-(51)
and (53)-(54) p.7-8, Table 1 p.13 (crop _audit-trail/phase2/crops/BW24SJ5K/p13.png).

Papadopoulos takes the shear amplitude C_a as the radius of the minimum circumscribed circle of the
shear path. For a single-harmonic path that is half the longest chord, the default of both models
(``np_shear_rule = "longest_chord"``, maintainer decision 17 Nov 2026); the RSS rule used before
reproduces the table only where the shear path on the critical plane is a line. Ties between planes of
equal C_a are resolved to the larger normal stress, as in the table (finding F10).

Rows not reproduced: row 9 Matake (closed form 162.20 + 0.2501 * 149.55 gives +1.7 %, the table prints
-1.4 %, McDiarmid on the same row matches), rows 17 and 28 (mean torsion with a phase shift, both
criteria off by 1 and 6-9 points, unresolved). The table prints values to 0.1 %, and up to 0.25 points
of scatter on in-phase rows whose closed form is exact suggests a coarse plane search in the source.
"""

import numpy as np
import pytest
from _paper import bending_torsion, load_table

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.matake import Matake
from multifat.models.cp_methods.mcdiarmid import McDiarmid
from multifat.planes.sampling import build_R_list

TABLE = load_table("papadopoulos1997_table1")
Z = np.zeros((3, 3))
R_LIST = build_R_list()


def _index(method, params, r):
    d = np.radians(r.delta)
    Sm = bending_torsion(r.sigma_m, r.tau_m)
    amp = bending_torsion(r.sigma_a, r.tau_a)
    best = best_cp_for_case(
        S0=Sm + amp,
        S1=Sm - amp,
        E0=Z,
        E1=Z,
        R_list=R_LIST,
        method=method,
        params=params,
        S_sin=bending_torsion(r.sigma_a, r.tau_a * np.cos(d)),
        S_cos=bending_torsion(0.0, -r.tau_a * np.sin(d)),
        E_sin=Z,
        E_cos=Z,
    )
    return 100 * (best.by_metric / r.t_m1 - 1)


def _rows(column, bad):
    out = []
    for r in TABLE.itertuples():
        if np.isnan(getattr(r, column)):
            continue
        if r.row in bad:
            out.append(pytest.param(r, marks=pytest.mark.xfail(strict=True, reason=bad[r.row])))
        else:
            out.append(r)
    return out


UNRESOLVED = "mean torsion with phase shift, unresolved"


@pytest.mark.parametrize(
    "r", _rows("I_matake", {9: "table calculation error", 17: UNRESOLVED, 28: UNRESOLVED})
)
def test_table1_matake(r):
    I = _index(Matake(), {"Sigm1": r.f_m1, "Taum1": r.t_m1}, r)
    assert I == pytest.approx(r.I_matake, abs=0.25)


@pytest.mark.parametrize("r", _rows("I_mcdiarmid", {17: UNRESOLVED, 28: UNRESOLVED}))
def test_table1_mcdiarmid(r):
    I = _index(McDiarmid(), {"Su": r.sigma_f, "Taum1": r.t_m1}, r)
    assert I == pytest.approx(r.I_mcdiarmid, abs=0.25)


def test_row9_matake_closed_form():
    """In phase: C_a = sqrt((sigma/2)^2 + tau^2), N_max = sigma/2, K = 2 t/f - 1 (Eq. 51)."""
    r = TABLE[TABLE.row == 9].iloc[0]
    K = 2 * r.t_m1 / r.f_m1 - 1
    I = 100 * ((np.hypot(r.sigma_a / 2, r.tau_a) + K * r.sigma_a / 2) / r.t_m1 - 1)
    assert I == pytest.approx(_index(Matake(), {"Sigm1": r.f_m1, "Taum1": r.t_m1}, r), abs=1e-6)
    assert I == pytest.approx(1.73, abs=0.01)


class MatakeRSS(Matake):
    np_shear_rule = "rss"


@pytest.mark.parametrize("row", [7, 21])
def test_rss_rule_differs_out_of_phase(row):
    """The RSS rule (the default before 17 Nov 2026) exceeds the circumscribed-circle amplitude
    of the source restatement out of phase, so it misses these rows; the default now matches."""
    r = TABLE[TABLE.row == row].iloc[0]
    p = {"Sigm1": r.f_m1, "Taum1": r.t_m1}
    assert _index(MatakeRSS(), p, r) != pytest.approx(_index(Matake(), p, r), abs=0.25)
    assert _index(Matake(), p, r) == pytest.approx(r.I_matake, abs=0.25)
