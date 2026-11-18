"""MWHS (N28): Wu, Li, Fang & Song (2017), Proc IMechE Part C, doi:10.1177/0954406217708214,
checked against the PDF supplied by the maintainer.

Eq. (19), p.10. The paper's results are figures only; Fig. 14 (p.11) plots, for the GH4169 test
eps_a = 0.612 %, gamma_a = 0.88 %, sigma_a = 788.6 MPa, tau_a = 404 MPa, the maximum shear strain
amplitude (a) and the second term of Eq. (19) (b) against the phase angle. Digitised values
(+-2 % reading accuracy) are checked with the measured strain amplitudes, nu = 0.48 (Table 3,
v_e) for the transverse strains, E = 198.5 GPa and sigma_y = 1083.1 MPa:
- (a) is reproduced at every phase angle with the longest-chord shear range, not with rss;
- (b) is reproduced from 0 to 75 deg; at 90 deg the evaluation gives 0.0039 against 0.0049 read
  from the figure (unresolved, strict xfail).
"""

import numpy as np
import pytest

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.mwhs import MWHS
from multifat.planes.sampling import build_R_list

E, SY, NU = 198500.0, 1083.1, 0.48
EPS_A, GAM_A, SIG_A, TAU_A = 0.00612, 0.0088, 788.6, 404.0
R = build_R_list()
# Fig. 14 digitised: phase angle -> (maximum shear strain amplitude, second term of Eq. 19)
FIG14 = {
    0: (0.0127, 0.0018),
    15: (0.0126, 0.0018),
    30: (0.0122, 0.0024),
    45: (0.0116, 0.0030),
    60: (0.0109, 0.0037),
    75: (0.0100, 0.0044),
    90: (0.0091, 0.0049),
}


class _Part(MWHS):
    """MWHS with one term of Eq. (19) as the plane damage (same plane rule)."""

    def __init__(self, part):
        self.part = part

    def evaluate_on_planes(self, **kw):
        total, d_gam = super().evaluate_on_planes(**kw)
        return (0.5 * d_gam if self.part == "shear" else total - 0.5 * d_gam), d_gam


def _tube(phi_deg):
    p = np.radians(phi_deg)

    def t(x, y):
        return np.array([[x, y, 0.0], [y, 0.0, 0.0], [0.0, 0.0, 0.0]])

    S_sin, S_cos = t(SIG_A, TAU_A * np.cos(p)), t(0.0, -TAU_A * np.sin(p))
    E_sin = np.diag([EPS_A, -NU * EPS_A, -NU * EPS_A]) + t(0.0, GAM_A / 2 * np.cos(p))
    E_cos = t(0.0, -GAM_A / 2 * np.sin(p))
    return S_sin, S_cos, E_sin, E_cos


def _term(phi, part):
    S_sin, S_cos, E_sin, E_cos = _tube(phi)
    best = best_cp_for_case(
        S0=S_sin, S1=-S_sin, E0=E_sin, E1=-E_sin, R_list=R, method=_Part(part),
        params={"Sy": SY, "E": E}, S_sin=S_sin, S_cos=S_cos, E_sin=E_sin, E_cos=E_cos,
    )  # fmt: skip
    return best.by_metric


@pytest.mark.parametrize("phi", sorted(FIG14))
def test_fig14a_maximum_shear_strain(phi):
    assert _term(phi, "shear") == pytest.approx(FIG14[phi][0], rel=0.02)


@pytest.mark.parametrize(
    "phi",
    [
        *sorted(set(FIG14) - {90}),
        pytest.param(
            90,
            marks=pytest.mark.xfail(strict=True, reason="0.0039 against 0.0049 read in Fig. 14b"),
        ),
    ],
)
def test_fig14b_normal_term(phi):
    assert _term(phi, "normal") == pytest.approx(FIG14[phi][1], rel=0.05, abs=1e-4)


def test_uniaxial_and_torsion_closed_forms():
    """Elastic strains: 45 deg plane in tension, the torsion plane carries no normal term."""
    model, s, t = MWHS(), 300.0, 200.0
    G = E / (2 * (1 + 0.3))

    def hooke(S):
        return ((1 + 0.3) * S - 0.3 * np.trace(S) * np.eye(3)) / E

    def run(S):
        return best_cp_for_case(
            S0=S, S1=-S, E0=hooke(S), E1=hooke(-S), R_list=R, method=model,
            params={"Sy": SY, "E": E},
        ).by_metric  # fmt: skip

    tension = (1.3 * s / E) + 0.5 * (1 + s / 2 / SY) * (s / E) * np.sqrt(0.7 / 2)
    assert run(np.diag([s, 0.0, 0.0])) == pytest.approx(tension, rel=1e-9)
    torsion = np.zeros((3, 3))
    torsion[0, 1] = torsion[1, 0] = t
    assert run(torsion) == pytest.approx(t / G, rel=1e-9)


def test_compressive_plane_has_no_normal_term():
    S = np.zeros((3, 3))
    S[2, 2] = -100.0
    e0 = np.zeros((3, 3))
    e0[0, 2] = e0[2, 0] = 0.001
    e0[2, 2] = 0.0005
    r = MWHS().evaluate_on_plane(S0r=S, S1r=S, E0r=e0, E1r=-e0, params={"Sy": SY, "E": E})
    assert r.damage == pytest.approx(0.002, rel=1e-12)
