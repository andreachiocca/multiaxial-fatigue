"""Robert, Fogue and Bahuaud (1994) integral fatigue criterion.

Published in *ASTM STP 1231* (1994), pp. 369-387, doi:10.1520/STP13959S (roadmap G41),
Eqs. (8), (15)-(25), Appendix::

    E_h = a * max_psi tau_{h,psi,a} + b * sigma_{hh,a} + d * sigma_{hh,m}
    E   = (1 / s_-1) * sqrt( (1 / (4 * pi)) * int_S E_h^2 dS ) = 1
    b^2 = (15 - sqrt(D)) / 2,   D = 225 * (1 - 8 * (k - 3)^2 / 25),   k = s_-1^2 / t_-1^2
    a^2 = (12 * k - 21 + b^2) / 2
    d   = (-(3 * b + 2 * a) + sqrt((3 * b + 2 * a)^2 + 3 * (60 * s_-1^2 / s_0^2 - 15))) / 3

Valid for ``0.577 < t_-1 / s_-1 < 0.866``, where ``max_psi tau_{h,psi,a}`` is the largest
half-projection of the shear stress path onto a line in plane ``h`` (for a single harmonic,
the radius of the minimum circumscribed circle) and ``c = 0`` (Eq. 25).

Repository choices:
- Source Eq. (16) prints ``1 / pi`` instead of ``1 / (4 pi)`` required by Eqs. (8) and
  (18)-(24), and the Appendix omits ``1 / 2`` under the root of ``a`` required by Eqs. (20)
  and (24); the implementation restores both corrected factors (and reads Eq. 18 as ``E^2``).
- Material constants outside 0.577 < t_-1 / s_-1 < 0.866 raise ``MaterialScopeError``: the
  model is skipped on that material (eligibility ``out_of_scope``), not failed.
- Spherical integration uses a 48 x 192 Gauss-Legendre grid (``GRID = (48, 192)``).
- Endurance limits (``Sigm1``, ``Taum1``, and ``Sig0`` or the Boehme et al. 2026 Eq. 1
  estimate of ``s_0`` from ``Su``) replace fixed-life strengths ``N``, returning the scalar
  validity index ``E`` (``I`` in Eq. 27) for master-curve life prediction.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import MaterialScopeError
from ..invariant_methods.stress_invariants import harmonic_parts
from .bohme_papuga import _plane_amplitudes, _sphere_grid, _strengths

GRID = (48, 192)


def rfb_parameters(f_m1: float, t_m1: float, f_0: float) -> tuple[float, float, float]:
    k = (f_m1 / t_m1) ** 2
    D = 225.0 * (1.0 - 8.0 * (k - 3.0) ** 2 / 25.0)
    if not (4.0 / 3.0 < k < 3.0) or D < 0.0:
        raise MaterialScopeError(
            f"RFB requires 0.577 < t_-1/s_-1 < 0.866 (received {t_m1 / f_m1:.3f}); "
            "see the Appendix of Robert et al."
        )
    b = np.sqrt((15.0 - np.sqrt(D)) / 2.0)
    a = np.sqrt((12.0 * k - 21.0 + b**2) / 2.0)
    combined = 3.0 * b + 2.0 * a
    d = (-combined + np.sqrt(combined**2 + 3.0 * (60.0 * f_m1**2 / f_0**2 - 15.0))) / 3.0
    return float(a), float(b), float(d)


def rfb_index(
    S_mean,
    S_sin,
    S_cos,
    f_m1: float,
    t_m1: float,
    f_0: float,
) -> float:
    a, b, d = rfb_parameters(f_m1, t_m1, f_0)
    normals, weights = _sphere_grid(*GRID)
    sigma_a, sigma_m, tau_a, _ = _plane_amplitudes(S_mean, S_sin, S_cos, normals)
    E_h = a * tau_a + b * sigma_a + d * sigma_m
    return float(np.sqrt(np.sum(weights * E_h**2)) / f_m1)


class RobertFogueBahuaud:
    """Eq. (16) with the 1/(4 pi) normalization taken from Eq. (8), returned as E (= I, Eq. 27)."""

    name = "RFB"
    optional_params = {
        "Sig0": "estimated from Su and Sigm1 by Eq. (1) of Boehme et al. 2026",
        "Su": "used only to estimate Sig0",
    }

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def evaluate_case(
        self,
        *,
        S0,
        E0,
        S1,
        E1,
        params: Mapping[str, float],
        R_list=None,
        harmonics=None,
    ) -> CaseResult:
        f_m1, t_m1, f_0, _ = _strengths(params, name="RFB")
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        return CaseResult(values={"": rfb_index(S_mean, S_sin, S_cos, f_m1, t_m1, f_0)})
