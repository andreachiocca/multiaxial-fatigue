"""MDC (roadmap G06): Margetin, Ďurka & Chmelko (2016) Frattura ed Integrità Strutturale 37:146-152,
doi:10.3221/IGF-ESIS.37.20 (Zotero VWAAQW5Q); Eqs. (4)-(8) p.148 (crop
_audit-trail/phase2/crops/VWAAQW5Q/p3_*.png), Tables 2 and 4 pp.149-150.

Eqs. (7)-(8) write sigma_max for the sigma_a,cr of Eq. (4). Reading it as the normal-stress amplitude
on the maximum-shear plane reproduces every Nf_com of Table 4 within 0.3 % in tau_r (life ratio
0.96-1.005 at the slope 1/b_tau = -13.6; the paper's own plane search is presumably coarser). Other
readings (sigma_max = applied sigma_a, or tau_a,cr = applied tau_a) miss by factors of 2-70, so the
Phase-1 note that Table 4 is not reproducible is withdrawn.

Side checks on the same table: repository MCD with t_A = 244 MPa, sigma_u = 680 MPa (Eq. 10, Table 3)
reproduces Nf_mcd within 4 % in life; repository FIN with k = 0.131 gives lives about 18 % shorter
than Nf_find (1.5 % in stress, unresolved).
"""

import numpy as np
import pytest
from _paper import bending_torsion, load_table

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.mcdiarmid import McDiarmid
from multifat.models.cp_methods.mdc import MDC, mdc_equivalent_shear, mdc_life
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

TABLE = load_table("margetin2016_table4")
P50 = {"sigma_f": 636.0, "b_fat": -0.0531, "tau_f": 550.0, "b0": -0.0736}
CONST = (636.0, -0.0531, 550.0, -0.0736)
Z = np.zeros((3, 3))


def _best(method, params, sigma, tau):
    S = bending_torsion(sigma, tau)
    return best_cp_for_case(
        S0=S, S1=-S, E0=Z, E1=Z, R_list=build_R_list(), method=method, params=params
    )


@pytest.mark.parametrize("N", [1e4, 1e5, 1e6])
def test_boundary_conditions(N):
    """Pure axial and pure torsion return the Basquin curves of Eqs. (5), (6)."""
    sf, bs, tf, bt = CONST
    s = sf * (2 * N) ** bs
    assert mdc_equivalent_shear(s / 2, s / 2, *CONST) == pytest.approx(
        tf * (2 * N) ** bt, rel=1e-12
    )
    t = tf * (2 * N) ** bt
    assert mdc_equivalent_shear(0.0, t, *CONST) == pytest.approx(t, rel=1e-12)
    (model,) = get_models(["MDC"])
    assert model.name == "MDC"
    assert mdc_life(_best(MDC(), P50, s, 0.0).by_metric, tf, bt) == pytest.approx(N, rel=1e-6)


@pytest.mark.parametrize("r", list(TABLE.itertuples()), ids=lambda r: f"n{r.n}")
def test_table4_nf_com(r):
    (model,) = get_models(["MDC"])
    tau_r = _best(MDC(), P50, r.sigma_a, r.tau_a).by_metric
    closed = mdc_equivalent_shear(r.sigma_a / 2, np.hypot(r.sigma_a / 2, r.tau_a), *CONST)
    assert tau_r == pytest.approx(closed, rel=1e-7)
    assert tau_r == pytest.approx(550.0 * (2 * r.Nf_com) ** -0.0736, rel=3e-3)


@pytest.mark.parametrize("r", list(TABLE.itertuples()), ids=lambda r: f"n{r.n}")
def test_table4_mcdiarmid_column(r):
    damage = _best(McDiarmid(), {"Su": 680.0, "Taum1": 244.0}, r.sigma_a, r.tau_a).by_metric
    assert mdc_life(damage, 550.0, -0.0736) == pytest.approx(r.Nf_mcd, rel=0.04)
