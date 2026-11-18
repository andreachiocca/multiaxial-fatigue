"""VSE: Liu (1993) virtual strain energy critical-plane parameters (roadmap G15).

K. C. Liu, ASTM STP 1191 (1993) 67-84, doi:10.1520/STP24796S, Eqs. (7)-(11)::

    dW = dW_n + dW_s = dsigma_n de_n + dtau d(gamma/2)                       (7)
    mode I  (Eq. 8): dW on the plane maximising dW_n
    mode II (Eq. 9): dW on the plane maximising dW_s

Under uniaxial loading dW = dsigma de (Eq. 11); in torsion the two modes coincide (p.73). The paper
gives no rule to choose the mode (it depends on material and temperature) and fits the life
curve, Eq. (6), to the data. ``VSE_ext`` is Eq. (7) maximised over all planes, which in general
differs from Eq. (8).

Repository choices:
  - one registered model, mode I unless the parameter ``vse_mode = 2`` selects mode II;
  - ranges are endpoint differences of the repository cycle; out-of-phase shear uses the
    longest-chord rule (Liu's out-of-phase rules, pp.80-81, are verbal only);
  - workbook strains are elastic, so the plastic term of Eq. (6) is absent; life from the shared
    master curve DP = A N^b instead of Eq. (6).

User guide: docs/method-guides/VSE.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
from multifat.models.utils.batch import check_planes

from ..base import CPPlaneResult


def _energy_parts(S0r, E0r, S1r, E1r):
    d_sigma_n = S0r[:, 2, 2] - S1r[:, 2, 2]
    d_eps_n = E0r[:, 2, 2] - E1r[:, 2, 2]
    dW_n = np.abs(d_sigma_n * d_eps_n)
    dW_s = np.abs(
        (S0r[:, 0, 2] - S1r[:, 0, 2]) * (E0r[:, 0, 2] - E1r[:, 0, 2])
        + (S0r[:, 1, 2] - S1r[:, 1, 2]) * (E0r[:, 1, 2] - E1r[:, 1, 2])
    )
    return dW_n, dW_s


class LiuVSE:
    """VSE critical-plane method; Eq. (8) (mode I) is used unless
    ``params['vse_mode'] = 2`` requests Eq. (9) (mode II)."""

    name = "VSE"
    optional_params = {"vse_mode": "1 (mode I, Eq. 8); 2 = mode II (Eq. 9)"}
    np_shear_rule = "longest_chord"

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
        check_planes(S0r, E0r, S1r, E1r)
        mode = params.get("vse_mode", None)
        mode = 1 if mode is None else int(mode)
        if mode not in (1, 2):
            raise ValueError(f"VSE: vse_mode must be 1 or 2, got {mode}")
        dW_n, dW_s = _energy_parts(S0r, E0r, S1r, E1r)
        return dW_n + dW_s, (dW_n if mode == 1 else dW_s)
