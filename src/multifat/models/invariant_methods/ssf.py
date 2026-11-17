"""Stress scale factor (SSF) equivalent shear stress of Anes, Reis, Li, Fonte & de Freitas (2014) (N24).

V. Anes, L. Reis, B. Li, M. Fonte, M. de Freitas, Int. J. Fatigue 62 (2014) 21-33,
doi:10.1016/j.ijfatigue.2013.05.004, Eqs. (29)-(31)::

    tau_eqv = max_t ( tau(t) + ssf(lambda, sigma) sigma(t) )                               (29)
    lambda = arctan(tau / sigma)                                                           (30)
    ssf(sigma, lambda) = a + b sigma + c sigma^2 + d sigma^3
                         + f lambda^2 + g lambda^3 + h lambda^4 + i lambda^5                (31)

with, for the 42CrMo4 of the paper (Reis batch, Table 2), a = 2.69, b = -9.90e-3, c = 1.69e-5,
d = -9.52e-9, f = -5.99, g = 11.72, h = -8.04, i = 1.63 (sigma in MPa, lambda in rad). The ssf
surface is regressed on proportional tension-torsion S-N trend lines at five amplitude ratios
(Tables 4-9): ssf = (tau_ref(N) - tau_a) / sigma_a, tau_ref the pure-shear curve (case 2). Under
non-proportional loading the instantaneous stresses enter Eq. (31) (Section 4.5, Figs. 12-14), and
tau_eqv is compared with the pure-shear S-N curve.

Finding: with the printed three-digit coefficients Eq. (31) misses the paper's own Table 9 ssf by
up to 0.06, and the case-6 peak of Fig. 14 gives 415 MPa against the 434 MPa of the text.

Repository choices:
  - maintainer decision: only the published surface; the eight coefficients are required
    constants (a_SSF..i_SSF, no e term), so without them the model is skipped. The values of the
    paper are shipped in Calibrations/42CrMo4_Reis/SSF_published.json for the 42CrMo4_Reis workbook,
    whose tests are the ones the surface was fitted on (in-sample);
  - maintainer decision: Eq. (31) as printed (no refit);
  - loading scope: fully reversed (no mean stress) tension-torsion, one normal and one shear
    component; any other stress state raises LoadingScopeError;
  - Eqs. (29)-(30) on magnitudes: lambda = arctan(|tau|/|sigma|) in [0, pi/2] and
    tau_eqv = max_t (|tau| + ssf |sigma|), so the result does not depend on the sign conventions
    of the axes (it equals the signed form of Fig. 14 whenever tau and sigma share the sign);
  - the cycle is sampled at SAMPLES points; a warning is issued when |sigma| exceeds the largest
    axial amplitude of the fit (737 MPa, Table 5);
  - life: master curve fitted on the pure torsion R = -1 failures (design_loading = "torsion"), the
    pure-shear reference curve of the paper.

User guide: docs/method-guides/SSF.md.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError
from .stress_invariants import harmonic_parts

SAMPLES = 4096  # points per cycle for the maximum of Eq. (29)
SIGMA_FIT_MAX = 737.0  # largest axial amplitude behind the ssf regression (Table 5), MPa
COEFFICIENTS = ("a_SSF", "b_SSF", "c_SSF", "d_SSF", "f_SSF", "g_SSF", "h_SSF", "i_SSF")
PUBLISHED = dict(zip(COEFFICIENTS, (2.69, -9.90e-3, 1.69e-5, -9.52e-9, -5.99, 11.72, -8.04, 1.63)))
_SCALE_TOL = 1e-9


def stress_scale_factor(sigma, lam, coefficients: Mapping[str, float]):
    """ssf of Eq. (31) for the axial stress magnitude sigma (MPa) and lam = arctan(tau/sigma) (rad)."""
    a, b, c, d, f, g, h, i = (float(coefficients[k]) for k in COEFFICIENTS)
    sigma, lam = np.asarray(sigma, float), np.asarray(lam, float)
    return (
        a
        + b * sigma
        + c * sigma**2
        + d * sigma**3
        + f * lam**2
        + g * lam**3
        + h * lam**4
        + i * lam**5
    )


def equivalent_shear_history(sigma, tau, coefficients: Mapping[str, float]) -> np.ndarray:
    """|tau| + ssf |sigma| at each instant of a sampled history (Eqs. 29-31)."""
    s, t = np.abs(np.asarray(sigma, float)), np.abs(np.asarray(tau, float))
    lam = np.arctan2(t, s)
    return t + stress_scale_factor(s, lam, coefficients) * s


def _tension_torsion_channels(S_mean, S_sin, S_cos):
    """(sigma_sin, sigma_cos, tau_sin, tau_cos) of a fully reversed tension-torsion harmonic."""
    scale = max(float(np.max(np.abs(S_sin))), float(np.max(np.abs(S_cos))), 1e-300)
    if np.max(np.abs(S_mean)) > _SCALE_TOL * scale:
        raise LoadingScopeError("SSF: mean stress; the ssf surface is fitted at R = -1 only")
    active = np.abs(S_sin) + np.abs(S_cos) > _SCALE_TOL * scale
    normals = [k for k in range(3) if active[k, k]]
    shears = [(i, j) for i, j in ((0, 1), (1, 2), (0, 2)) if active[i, j]]
    if len(normals) > 1 or len(shears) > 1 or (normals and shears and normals[0] not in shears[0]):
        raise LoadingScopeError(
            "SSF: only one normal and one shear stress acting in the same plane (tension-torsion)"
        )
    n = normals[0] if normals else 0
    i, j = shears[0] if shears else (0, 1)
    return S_sin[n, n], S_cos[n, n], S_sin[i, j], S_cos[i, j]


def ssf_equivalent_shear(S_mean, S_sin, S_cos, coefficients: Mapping[str, float]) -> float:
    """tau_eqv of Eq. (29) for one single-harmonic tension-torsion cycle (MPa)."""
    S_mean, S_sin, S_cos = (np.asarray(x, float) for x in (S_mean, S_sin, S_cos))
    s_sin, s_cos, t_sin, t_cos = _tension_torsion_channels(S_mean, S_sin, S_cos)
    phase = np.linspace(0.0, 2.0 * np.pi, SAMPLES, endpoint=False)
    sigma = s_sin * np.sin(phase) + s_cos * np.cos(phase)
    tau = t_sin * np.sin(phase) + t_cos * np.cos(phase)
    if np.max(np.abs(sigma)) > SIGMA_FIT_MAX:
        warnings.warn(
            f"SSF: axial stress {np.max(np.abs(sigma)):.0f} MPa beyond the ssf fit "
            f"(<= {SIGMA_FIT_MAX:.0f} MPa); Eq. (31) is extrapolated",
            RuntimeWarning,
            stacklevel=2,
        )
    return float(np.max(equivalent_shear_history(sigma, tau, coefficients)))


class SSF:
    """Anes et al. (2014) stress scale factor equivalent shear stress (see the module docstring)."""

    name = "SSF"
    # Eq. (29) is compared with the pure-shear S-N curve: master curve on pure torsion failures.
    design_loading = "torsion"

    def required_params(self) -> set[str]:
        return set(COEFFICIENTS)

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": ssf_equivalent_shear(S_mean, S_sin, S_cos, params)})
