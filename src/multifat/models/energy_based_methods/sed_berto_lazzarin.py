"""Averaged strain energy density (SED) criterion of Berto & Lazzarin (2011).

Published in *International Journal of Fatigue* 33 (2011), pp. 1055-1065,
doi:10.1016/j.ijfatigue.2010.11.019 (roadmap G30), Eqs. (7)-(8)::

    Delta W = (1 / E) * (e1 * Delta K1^2 / R1^(2 * (1 - lambda1))
                         + e3 * Delta K3^2 / R3^(2 * (1 - lambda3)))
    R1 = (sqrt(2 * e1) * Delta K1A / Delta sigma1A)^(1 / (1 - lambda1))
    R3 = (sqrt(e3 / (1 + nu)) * Delta K3A / Delta tau3A)^(1 / (1 - lambda3))

where ``Delta K`` are Mode I and Mode III notch stress intensity factor ranges,
``lambda1, lambda3`` are the Williams eigenvalues, and ``e1, e3`` are opening-angle SED
shape constants calibrated at shared life ``N_A``.

Repository choices:
- Deliberately exposed as standalone notch kernels (``control_radii``, ``averaged_sed``)
  and not registered as a workbook model because smooth specimens carry no singular stress
  field (``Delta K1, Delta K3``), causing control radii ``R1, R3`` to vanish.
- Does not add a separate phase-shift term or apply mean-stress factor ``c_W`` automatically
  beyond explicit caller inputs.
"""

from __future__ import annotations

import numpy as np


def control_radii(
    dK1A: float,
    dsigma1A: float,
    dK3A: float,
    dtau3A: float,
    e1: float,
    e3: float,
    lambda1: float,
    lambda3: float,
    nu: float,
) -> tuple[float, float]:
    """Return the control radii R1 and R3 defined by Eqs. (8a) and (8b)."""
    R1 = (np.sqrt(2.0 * e1) * dK1A / dsigma1A) ** (1.0 / (1.0 - lambda1))
    R3 = (np.sqrt(e3 / (1.0 + nu)) * dK3A / dtau3A) ** (1.0 / (1.0 - lambda3))
    return float(R1), float(R3)


def averaged_sed(
    dK1: float,
    dK3: float,
    R1: float,
    R3: float,
    e1: float,
    e3: float,
    lambda1: float,
    lambda3: float,
    E: float,
) -> float:
    """Delta W from Eq. (7) (energy per unit volume, i.e. NSIF^2 / (E length))."""
    mode1_contribution = e1 * dK1**2 / R1 ** (2.0 * (1.0 - lambda1))
    mode3_contribution = e3 * dK3**2 / R3 ** (2.0 * (1.0 - lambda3))
    return float((mode1_contribution + mode3_contribution) / E)
