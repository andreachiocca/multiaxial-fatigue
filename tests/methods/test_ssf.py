"""SSF (N24): Anes, Reis, Li, Fonte & de Freitas (2014), Int J Fatigue 62:21-33,
doi:10.1016/j.ijfatigue.2013.05.004, checked against the PDF supplied by the maintainer.

Eqs. (29)-(31); Table 4 (S-N trend lines, p.27), Tables 5-9 (ssf calibration, pp.28-29), Section
4.5 and Figs. 13-14 (non-proportional case 6, p.31).

Findings:
- Tables 5-8 follow exactly from the Table 4 trend lines (ssf = (AC - AB)/AD), except two printed
  cells: Table 5 prints BC = 0 for case 1 (it is AC - AB = AC; its ssf column uses AC) and Table 8
  prints AC = 508 at 1e3 cycles (it is 567; its BC = 59 uses 567). Strict xfails;
- Eq. (31) with the printed three-digit coefficients misses Table 9 by up to 0.06;
- case 6 (sigma_a = 510, tau_a = 294 MPa, 90 deg): the peak of Eq. (29) is 415 MPa against the
  434 MPa of the text (strict xfail); the ssf history spans 0.37-0.79 as in Fig. 13.
"""

import numpy as np
import pytest
from _paper import load_table

from multifat.models.errors import LoadingScopeError, MissingParamsError
from multifat.models.invariant_methods.ssf import (
    PUBLISHED,
    SSF,
    ssf_equivalent_shear,
    stress_scale_factor,
)
from multifat.models.registry import get_models

TABLE = load_table("anes2014_tables5_9")
TREND = {1: (1204.3, -0.071, 0.0, 0.0), 3: (938.81, -0.052, 312.94, -0.052),
         4: (1025.7, -0.066, 592.3, -0.066), 5: (1132.5, -0.116, 1132.5, -0.116)}  # fmt: skip
MISPRINTS = {(1, 1e3, "BC"), (1, 1e4, "BC"), (1, 1e5, "BC"), (1, 1e6, "BC"), (1, 1e7, "BC"),
             (5, 1e3, "AC")}  # fmt: skip
Z = np.zeros((3, 3))


def _shear_ref(N):
    return 864.78 * N**-0.061  # Table 4, case 2 (pure shear)


def _cells():
    out = []
    for r in TABLE.itertuples():
        A, b, B, bb = TREND[r.case]
        AB = B * r.N**bb if B else 0.0
        computed = dict(AC=_shear_ref(r.N), AB=AB, AD=A * r.N**b)
        computed["BC"] = computed["AC"] - AB
        computed["ssf"] = computed["BC"] / computed["AD"]
        for col, value in computed.items():
            key = (r.case, r.N, col)
            mark = [pytest.mark.xfail(strict=True, reason="misprint")] if key in MISPRINTS else []
            out.append(pytest.param(getattr(r, col), value, col, marks=mark,
                                    id=f"case{r.case}-N{r.N:.0e}-{col}"))  # fmt: skip
    return out


@pytest.mark.parametrize("printed, computed, col", _cells())
def test_tables5_8_from_table4_trend_lines(printed, computed, col):
    """The calibration procedure: ssf = (tau_ref(N) - tau_a(N)) / sigma_a(N)."""
    assert computed == pytest.approx(printed, abs=0.006 if col == "ssf" else 0.6)


def test_table9_eq31_printed_coefficients():
    """Eq. (31) as printed reproduces Table 9 only to 0.06 (coefficient rounding)."""
    fitted = stress_scale_factor(TABLE["AD"], TABLE["lambda"], PUBLISHED)
    gap = np.abs(fitted - TABLE["ssf"])
    assert gap.max() == pytest.approx(0.06, abs=0.01)
    assert np.sqrt(np.mean(gap**2)) < 0.035


def _case(sigma_a, tau_a, delta=0.0):
    d = np.radians(delta)
    S_sin = np.array([[0, tau_a * np.cos(d), 0], [tau_a * np.cos(d), sigma_a, 0], [0, 0, 0]])
    S_cos = np.array([[0, -tau_a * np.sin(d), 0], [-tau_a * np.sin(d), 0, 0], [0, 0, 0]])
    return Z, S_sin, S_cos


def test_case6_ssf_history_fig13():
    """Fig. 13: along the 90 deg cycle the ssf oscillates between about 0.38 and 0.8, with sharp
    drops to zero where lambda reaches pi/2 (sigma = 0); the drops are left out here."""
    phase = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    s, t = np.abs(510 * np.sin(phase)), np.abs(294 * np.cos(phase))
    lam = np.arctan2(t, s)
    ssf = stress_scale_factor(s, lam, PUBLISHED)[lam < 1.2]
    assert ssf.min() == pytest.approx(0.37, abs=0.02)
    assert ssf.max() == pytest.approx(0.79, abs=0.02)


def test_case6_equivalent_stress_within_coefficient_gap():
    """Text p.31: 434 MPa; the printed coefficients give 415 MPa (Table 9 gap: 0.06 x 510 MPa)."""
    tau_eqv = ssf_equivalent_shear(*_case(510.0, 294.0, 90.0), PUBLISHED)
    assert tau_eqv == pytest.approx(434.0, abs=0.06 * 510.0)


@pytest.mark.xfail(strict=True, reason="printed Eq. (31) coefficients give 415 MPa, text 434 MPa")
def test_case6_equivalent_stress_as_printed():
    assert ssf_equivalent_shear(*_case(510.0, 294.0, 90.0), PUBLISHED) == pytest.approx(
        434.0, rel=0.01
    )


@pytest.mark.parametrize("sigma_a, tau_a", [(500.0, 0.0), (480.0, 277.0), (298.0, 298.0)])
def test_proportional_peak(sigma_a, tau_a):
    """Proportional loading: tau_eqv = tau_a + ssf(sigma_a, lambda) sigma_a at the peak."""
    lam = np.arctan2(tau_a, sigma_a)
    expected = tau_a + stress_scale_factor(sigma_a, lam, PUBLISHED) * sigma_a
    assert ssf_equivalent_shear(*_case(sigma_a, tau_a), PUBLISHED) == pytest.approx(
        expected, rel=1e-6
    )


def test_pure_torsion_is_the_reference():
    assert ssf_equivalent_shear(*_case(0.0, 300.0), PUBLISHED) == pytest.approx(300.0, rel=1e-9)


def test_scope_and_missing_constants():
    (model,) = get_models(["SSF"])
    with pytest.raises(LoadingScopeError, match="mean stress"):
        ssf_equivalent_shear(np.diag([0, 50.0, 0]), *_case(300.0, 100.0)[1:], PUBLISHED)
    with pytest.raises(LoadingScopeError, match="tension-torsion"):
        ssf_equivalent_shear(Z, np.diag([300.0, 200.0, 0]), Z, PUBLISHED)
    assert SSF().required_params() == set(PUBLISHED)
    S = _case(300.0, 100.0)[1]
    with pytest.raises((MissingParamsError, KeyError)):
        model.evaluate_case(S0=S, S1=-S, E0=Z, E1=Z, params={})


def test_published_constants_record():
    """Calibrations/42CrMo4_Reis/SSF_published.json loads through --constants as printed."""
    from pathlib import Path

    from multifat.cli.pipeline import load_constants

    path = Path(__file__).resolve().parents[2] / "Calibrations/42CrMo4_Reis/SSF_published.json"
    constants, origin = load_constants(path, "42CrMo4_Reis")
    assert constants == PUBLISHED
    assert origin["kind"] == "published constants"
    with pytest.raises(ValueError, match="not '42CrMo4_Santus'"):
        load_constants(path, "42CrMo4_Santus")
