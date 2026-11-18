"""Ninic (2006) stress-based multiaxial HCF Damage Function, a critical-plane criterion (N11).

D. Ninic, Int. J. Fatigue 28 (2006) 103-113, doi:10.1016/j.ijfatigue.2005.04.014, Eqs. (1)-(4)::

    D(l, m, n) = [ (tau_a / T_e)^2 + k (sigma_n,max / S_e)^2 ]^(1/2)                    (1)
    D = max over all planes (l, m, n) <= 1                                               (2)-(3)
    k = (S_e / T_e)^2 [1 - (S_e / (2 T_e))^2]                                            (4)

T_e and S_e are the fully reversed torsion and tension (bending) endurance limits; the plane of
maximum D is the critical plane, in all three dimensions (planes tilted into the surface
included). The error index is I = 100 (D - 1), Eq. (14).

Eq. (4) puts the uniaxial limit at D = 1 only when S_e/T_e >= sqrt(2): the uniaxial maximum of
Eq. (1) lies at cos^2(theta) = 2 (T_e/S_e)^2, which leaves [0, 1] for S_e/T_e < sqrt(2). Below
that ratio the maximum is on the principal plane and the uniaxial limit gives D = sqrt(k) < 1
(cast iron of Table 4, row 5: I = -10.4 %). The torsion limit always gives D = 1.

Repository choices:
  - tau_a is the shear amplitude on the plane by the default ``rss`` rule (the paper gives no
    out-of-phase definition; Tables 2-4 are reproduced with it and not with the longest chord);
  - sigma_n,max is clipped at zero, so a compressive normal stress, squared in Eq. (1), never adds
    damage (the paper's tables have tensile maxima only);
  - k follows Eq. (4) unless a workbook value ``k_NIN`` is supplied; a warning is issued when
    S_e/T_e < sqrt(2), where Eq. (4) does not calibrate the uniaxial limit;
  - D is both the plane metric and the damage, so NINIC and NINIC_ext coincide; life from the
    shared master curve DP = A N^b (an extension of the endurance criterion).

User guide: docs/method-guides/NINIC.md.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult


def ninic_k(f_m1: float, t_m1: float) -> float:
    """Normal stress sensitivity factor k of Eq. (4) from S_e = f_-1 and T_e = t_-1."""
    kappa = float(f_m1) / float(t_m1)
    return kappa**2 * (1.0 - (kappa / 2.0) ** 2)


class Ninic:
    """Ninic Damage Function; the critical plane maximises D itself."""

    name = "NINIC"
    optional_params = {"k_NIN": "Eq. (4) from Sigm1 and Taum1"}

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1"}

    def _constants(self, params: Mapping[str, float]) -> tuple[float, float, float]:
        f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
        if f_m1 <= 0.0 or t_m1 <= 0.0:
            raise ValueError(f"NINIC needs positive limits, got Sigm1 = {f_m1}, Taum1 = {t_m1}")
        k_value = params.get("k_NIN")
        if k_value is None:
            if f_m1 / t_m1 < np.sqrt(2.0):
                warnings.warn(
                    f"NINIC: S_e/T_e = {f_m1 / t_m1:.3f} < sqrt(2); Eq. (4) then gives "
                    "D = sqrt(k) < 1 at the uniaxial limit",
                    RuntimeWarning,
                    stacklevel=3,
                )
            k_value = ninic_k(f_m1, t_m1)
        k_value = float(k_value)
        if k_value < 0.0:
            raise ValueError(f"NINIC: k must be >= 0, got {k_value} (S_e/T_e > 2 in Eq. 4)")
        return f_m1, t_m1, k_value

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
        """All planes at once: D of Eq. (1) as damage and metric."""
        check_planes(S0r, E0r, S1r, E1r)
        f_m1, t_m1, k_value = self._constants(params)
        tau_a = 0.5 * np.hypot(S0r[:, 0, 2] - S1r[:, 0, 2], S0r[:, 1, 2] - S1r[:, 1, 2])
        sigma_n_max = np.maximum(np.maximum(S0r[:, 2, 2], S1r[:, 2, 2]), 0.0)
        damage = np.sqrt((tau_a / t_m1) ** 2 + k_value * (sigma_n_max / f_m1) ** 2)
        return damage, damage
