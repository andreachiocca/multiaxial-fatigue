"""Wang, Sun, Liu, Liu & Shang (2022) multiaxial HCF damage parameter and life model (N05).

Y. Wang, G.-Q. Sun, J. Liu, X. Liu, D. Shang, Int. J. Damage Mech. 31 (2022) 1165-1186,
doi:10.1177/10567895221091308, Eqs. (7)-(11), (13), (15), (18)-(24), on the plane of maximum shear
stress amplitude C_a::

    t_eq,a = sqrt( C_a^2 + k (1 - N*_max / (2 sigma_y))^2 N_a^2 ) <= lambda = t_-1         (18)
    N*_max = sqrt(N_max N_a),   N_max = N_a + N_m        (SWT stress function, Eq. 15)
    k = (4 t_-1^2 - f_-1^2) / (f_-1^2 (1 - f_-1 / (4 sigma_y))^2)                          (19)
    t_eq,a = f(N_f): torsional S-N curve, Basquin tau'_f (2N_f)^b1 or t_-1/(1 - kappa N_f^b2)  (24)

C_a is the major semi-axis of the shear path on the plane (Eq. 10, the minimum circumscribed
ellipse method); N_a, N_m the normal stress amplitude and mean on the plane; sigma_y the yield
strength. k follows from the two fully reversed limits; k > 0 requires t_-1/f_-1 > 0.5 (Eq. 21).

Finding: with Eq. (18) as printed, a tensile mean normal stress raises N*_max and so LOWERS
t_eq,a (a beneficial tensile mean), contrary to the paper's text; the printed form nevertheless
reproduces exactly the predicted lives of Tables 4-5 (no mean stress), a '+' sign reproduces none,
and the mean-stress lives of Table 7 are reproduced by neither.

Repository choices:
  - maintainer decision: Eq. (18) exactly as printed, with loading that puts a mean normal stress
    on the critical plane declared out of scope (LoadingScopeError, model skipped). The critical
    plane is found first with the same rule, and |N_m| there is checked;
  - C_a by the longest-chord rule (Eq. 10); ties of the maximum C_a towards the larger t_eq,a;
  - master curve fitted on the pure torsion R = -1 failures (``design_loading = "torsion"``),
    as the torsional S-N curve of Eq. (24); a single power law replaces the three-parameter form;
  - registered as a case-level model (the scope check needs the critical plane), so it has no
    _ext variant.

User guide: docs/method-guides/WANG2022.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CaseResult, CPPlaneResult
from ..cp_adapter import CPModelAdapter
from ..errors import LoadingScopeError

MEAN_TOL = 1e-9  # relative tolerance on |N_m| at the critical plane (fully reversed loading)


def wang_k(f_m1: float, t_m1: float, sy: float) -> float:
    """Material constant k of Eq. (19)."""
    return (4.0 * t_m1**2 - f_m1**2) / (f_m1**2 * (1.0 - f_m1 / (4.0 * sy)) ** 2)


def _plane_parts(S0r, S1r):
    C_a = 0.5 * np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
    N_a = 0.5 * np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
    N_m = 0.5 * (S0r[:, 2, 2] + S1r[:, 2, 2])
    return C_a, N_a, N_m


class _WangPlane:
    """Plane-level Eq. (18): damage t_eq,a, metric C_a."""

    name = "WANG2022"
    np_shear_rule = "longest_chord"

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1", "Sy"}

    def evaluate_on_plane(self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float]):
        damage, metric = self.evaluate_on_planes(
            S0r=np.asarray(S0r, float)[None],
            E0r=np.asarray(E0r, float)[None],
            S1r=np.asarray(S1r, float)[None],
            E1r=np.asarray(E1r, float)[None],
            params=params,
        )
        return CPPlaneResult(damage=float(damage[0]), metric=float(metric[0]))

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        check_planes(S0r, E0r, S1r, E1r)
        f_m1, t_m1, sy = (float(params[k]) for k in ("Sigm1", "Taum1", "Sy"))
        k = wang_k(f_m1, t_m1, sy)
        if not k > 0.0:
            raise ValueError(f"WANG2022: k = {k:.4g} <= 0 needs t_-1/f_-1 > 0.5 (Eq. 21)")
        C_a, N_a, N_m = _plane_parts(S0r, S1r)
        N_star = np.sqrt(np.maximum(N_a + N_m, 0.0) * N_a)
        return np.sqrt(C_a**2 + k * (1.0 - N_star / (2.0 * sy)) ** 2 * N_a**2), C_a


class _CriticalMean(_WangPlane):
    """Same plane rule; damage |N_m|, to read the mean normal stress on the critical plane."""

    def evaluate_on_planes(
        self, *, S0r, E0r, S1r, E1r, params: Mapping[str, float], tau_max=None, strain_ranges=None
    ):
        check_planes(S0r, E0r, S1r, E1r)
        C_a, _, N_m = _plane_parts(S0r, S1r)
        return np.abs(N_m), C_a


class Wang2022:
    """Case-level WANG2022: scope check on the critical plane, then Eq. (18)."""

    name = "WANG2022"
    # Eq. (24) is the torsional S-N curve: the master curve is fitted on pure torsion failures.
    design_loading = "torsion"

    def __init__(self):
        self._plane = CPModelAdapter(_WangPlane())
        self._mean = CPModelAdapter(_CriticalMean())

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1", "Sy"}

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        kwargs = dict(S0=S0, E0=E0, S1=S1, E1=E1, params=params, R_list=R_list, harmonics=harmonics)
        mean_on_plane = self._mean.evaluate_case(**kwargs).values[""]
        scale = max(float(np.max(np.abs(S0))), float(np.max(np.abs(S1))), 1e-300)
        if mean_on_plane > MEAN_TOL * scale:
            raise LoadingScopeError(
                "WANG2022: mean normal stress on the critical plane; Eq. (18) as printed makes a "
                "tensile mean beneficial, so such loading is out of scope"
            )
        return CaseResult(values={"": self._plane.evaluate_case(**kwargs).values[""]})
