"""Zhang & Yao (2011) critical-plane multiaxial high-cycle fatigue criterion (N06).

C.-C. Zhang, W.-X. Yao, Fatigue Fract. Engng Mater. Struct. 34 (2011) 337-344,
doi:10.1111/j.1460-2695.2010.01523.x, Eqs. (3), (6), (7) and (A1), on the plane of maximum shear
stress amplitude C_a (as McDiarmid's criterion)::

    sqrt( C_a^2 - C_a N* + (t_-1^2 / f_-1^2) (f_-1 + 2 C_a) N* ) <= t_-1               (6)
    N* = N_a / (1 - N_m / sigma_u)                                                     (7)

C_a is the shear stress amplitude on the plane (Eq. A1: the major semi-axis of the shear ellipse),
N_a and N_m the normal stress amplitude and mean on that plane, sigma_u the ultimate tensile
strength. Eq. (6) meets both fully reversed limits (bending f_-1 and torsion t_-1); Eq. (9) and
Papadopoulos' bound restrict it to 1.25 < f_-1/t_-1 < 2. The returned fatigue index is
LHS / t_-1 (1 on the limit); the paper's error index is I = 100 (index - 1), Eq. (10).

Repository choices:
  - C_a by the longest-chord rule (the major semi-axis of Eq. A1); tables are not reproduced
    with rss;
  - ties of the maximum C_a are resolved towards the larger index (repository rule). The paper
    gives no tie rule, and the index varies by many percent over planes within 1e-4 of the
    maximum C_a, most of all in 90 deg out-of-phase loading, where the shear amplitude can be
    equal on a continuum of planes; the printed values mostly sit at the top of that range;
  - sigma_u (``Su``) is needed only where a plane carries a mean normal stress; without it such a
    case raises MissingParamsError. N_m >= sigma_u, or a negative radicand in Eq. (6), raises
    LoadingScopeError;
  - a warning outside 1.25 < f_-1/t_-1 < 2;
  - the fatigue-limit index is used as the damage parameter of the shared master curve
    DP = A N^b for finite-life comparisons (not part of the criterion).

User guide: docs/method-guides/ZHANGYAO.md.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult
from ..errors import LoadingScopeError, MissingParamsError


def zhang_yao_lhs(C_a, N_star, f_m1: float, t_m1: float):
    """Left-hand side of Eq. (6) (MPa) for shear amplitude C_a and corrected normal stress N*."""
    C_a, N_star = np.asarray(C_a, float), np.asarray(N_star, float)
    radicand = C_a**2 - C_a * N_star + (t_m1 / f_m1) ** 2 * (f_m1 + 2.0 * C_a) * N_star
    if np.any(radicand < 0.0):
        raise LoadingScopeError("ZHANGYAO: negative radicand in Eq. (6) (outside its validity)")
    return np.sqrt(radicand)


class ZhangYao:
    """Zhang-Yao criterion on the plane of maximum shear stress amplitude."""

    name = "ZHANGYAO"
    np_shear_rule = "longest_chord"
    optional_params = {"Su": "needed only when a plane carries a mean normal stress"}

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

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
        self,
        *,
        S0r,
        E0r,
        S1r,
        E1r,
        params: Mapping[str, float],
        tau_max=None,
        strain_ranges=None,
    ):
        """All planes at once: index LHS/t_-1 of Eq. (6) and the shear amplitude C_a."""
        check_planes(S0r, E0r, S1r, E1r)
        f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
        if not (f_m1 > 0.0 and t_m1 > 0.0):
            raise ValueError(f"ZHANGYAO needs positive limits, got {f_m1}, {t_m1}")
        kappa = f_m1 / t_m1
        if not (1.25 < kappa < 2.0):
            warnings.warn(
                f"ZHANGYAO: f_-1/t_-1 = {kappa:.3f} outside the validity range (1.25, 2)",
                RuntimeWarning,
                stacklevel=3,
            )
        C_a = 0.5 * np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        N_a = 0.5 * np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
        N_m = 0.5 * (S0r[:, 2, 2] + S1r[:, 2, 2])
        su = params.get("Su")
        if su is None:
            if np.any(np.abs(N_m) > 1e-9 * max(1.0, float(np.max(np.abs(S0r))))):
                raise MissingParamsError(
                    self.name, ["Su"], "sigma_u corrects the mean normal stress, Eq. (7)"
                )
            N_star = N_a
        else:
            su = float(su)
            if np.any(N_m >= su):
                raise LoadingScopeError(
                    "ZHANGYAO: mean normal stress >= sigma_u, Eq. (7) undefined"
                )
            N_star = N_a / (1.0 - N_m / su)
        return zhang_yao_lhs(C_a, N_star, f_m1, t_m1) / t_m1, C_a
