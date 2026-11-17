"""Non-proportional shear range on a plane: per-source rule (finding F3).

FS (1988), MKBM (Li et al. 2011, Eq. 19), LI (Li et al. 2021, Eqs. 17-19) and LIU2021 (Eq. 6)
define the plane shear range over the cycle, i.e. the longest chord of the shear path (also MWHS,
whose Fig. 14 needs it, and LIU2025, the half range of the plane shear); the other
critical-plane methods keep the repository's RSS of the ellipse semi-axes.
"""

import numpy as np
import pytest

from multifat.eval.evaluator import prepare_planes, shear_rule_of
from multifat.models.registry import available_cp_methods

CHORD = {"FS", "MKBM", "LI", "LIU2021", "MWCM", "MMCCM", "VSE", "MWHS", "LIU2025", "ZHANGYAO"}


def test_rule_assignment():
    rules = {n: shear_rule_of(m) for n, m in available_cp_methods().items()}
    assert {n for n, r in rules.items() if r == "longest_chord"} == CHORD
    assert all(r == "rss" for n, r in rules.items() if n not in CHORD)


@pytest.mark.parametrize("rule,expected", [("rss", np.hypot(2.0, 1.0)), ("longest_chord", 2.0)])
def test_elliptical_path_ranges(rule, expected):
    """Shear ellipse with semi-axes a=2 (sin) and b=1 (cos, 90 deg out of phase) on plane z."""
    z = np.zeros((3, 3))
    s_sin, s_cos = z.copy(), z.copy()
    s_sin[0, 2] = s_sin[2, 0] = 2.0
    s_cos[1, 2] = s_cos[2, 1] = 1.0
    p = prepare_planes(
        S0=s_sin,
        S1=-s_sin,
        E0=z,
        E1=z,
        R_list=np.eye(3)[None],
        S_sin=s_sin,
        S_cos=s_cos,
        shear_rule=rule,
    )
    d = p.S0r[0] - p.S1r[0]
    assert np.hypot(d[0, 2], d[1, 2]) == pytest.approx(2 * expected, rel=1e-14)


def test_rules_agree_for_proportional_paths():
    z = np.zeros((3, 3))
    s_sin = z.copy()
    s_sin[0, 2] = s_sin[2, 0] = 1.5
    s_sin[1, 2] = s_sin[2, 1] = -0.5
    a, b = (
        prepare_planes(
            S0=s_sin,
            S1=-s_sin,
            E0=z,
            E1=z,
            R_list=np.eye(3)[None],
            S_sin=s_sin,
            S_cos=z,
            shear_rule=r,
        )
        for r in ("rss", "longest_chord")
    )
    np.testing.assert_allclose(a.S0r, b.S0r, rtol=1e-14, atol=1e-15)
