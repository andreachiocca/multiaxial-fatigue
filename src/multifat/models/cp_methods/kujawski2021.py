"""Kujawski (2021) negative-mean-stress modified SWT function mSWT, critical-plane extension (N26).

D. Kujawski, Fatigue Fract. Eng. Mater. Struct. 44 (2021) 2202-2209, doi:10.1111/ffe.13488,
Eqs. (9)-(13)::

    W_D = sigma_max eps_a                       sigma_m >= 0   (SWT)               (11a)
    W_D = (sigma_max + |sigma_m| / 3) eps_a     sigma_m < 0    (mSWT)              (11b)
    sigma_ar = sqrt(sigma_max sigma_a)  or  sqrt((sigma_max + |sigma_m|/3) sigma_a)   (12a, 12b)

W_D is the total damaging energy density in the first quadrant of the damaging axes (Eq. 9);
|sigma_m|/3 is the stabilised part of a negative mean stress (Fig. 7). With sigma_m >= 0 the
function is the original SWT; with sigma_m < 0 the SWT maximum stress is augmented by |sigma_m|/3.

The paper is uniaxial and validated on figures only (ductile iron, 7075-T651).

Repository choices:
  - maintainer decision: the critical-plane extension of the repository's SWT, with sigma_n,max,
    sigma_n,m and eps_n,a of the plane normal; DP = W_D on the plane of largest normal strain
    range, KUJAWSKI2021_ext the largest DP over all planes. Without a compressive normal mean the
    result is identical to SWT;
  - the bracket of Eq. (11b), sigma_a - 2|sigma_m|/3, is clipped at zero (like the SWT clip of a
    fully compressive cycle), so DP = 0 when |sigma_m| >= 3 sigma_a / 2;
  - linear-elastic strains; the applied mean stress is used, not a measured stabilised one;
  - DP = W_D (MPa) rather than sigma_ar; life from the shared master curve DP = A N^b.

User guide: docs/method-guides/KUJAWSKI2021.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult


def mswt_stress(sigma_max, sigma_m):
    """Stress factor of Eqs. (11a)-(11b): sigma_max, augmented by |sigma_m|/3 when sigma_m < 0.

    Parameters
    ----------
    sigma_max, sigma_m : float or ndarray
        Maximum and mean normal stress of the cycle (MPa).

    Returns
    -------
    float or ndarray
        The bracket of Eq. (11b) (or sigma_max), clipped at zero (MPa).
    """
    sigma_max = np.asarray(sigma_max, float)
    sigma_m = np.asarray(sigma_m, float)
    s = np.where(sigma_m < 0.0, sigma_max + np.abs(sigma_m) / 3.0, sigma_max)
    return np.maximum(s, 0.0)


class Kujawski2021:
    """Plane-level mSWT: damage W_D (MPa), metric the normal strain range."""

    name = "KUJAWSKI2021"

    def required_params(self) -> set[str]:
        return set()

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
        """All planes at once (normal = local z): DP = mSWT(sigma_max, sigma_m) eps_a."""
        check_planes(S0r, E0r, S1r, E1r)
        sigma_max = np.maximum(S0r[:, 2, 2], S1r[:, 2, 2])
        sigma_m = 0.5 * (S0r[:, 2, 2] + S1r[:, 2, 2])
        d_eps = np.abs(E0r[:, 2, 2] - E1r[:, 2, 2])
        return mswt_stress(sigma_max, sigma_m) * d_eps / 2.0, d_eps
