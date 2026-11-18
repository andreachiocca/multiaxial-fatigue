"""Volumetric energy fatigue criterion of Banvillet, Palin-Luc & Lasserre (2003).

Published in *International Journal of Fatigue* 25 (2003), pp. 755-769,
doi:10.1016/S0142-1123(03)00048-3 (roadmap G26), Eqs. (3)-(29)::

    g(C) = mean over V* of [W_g - W*] <= g^D,   V* = {M : W_g(M) >= W*}
    W* = sigma*^2 / E,   sigma* = sqrt(2 * s_Ten^2 - s_RB^2),   g^D_uniax = (s_RB^2 - s_Ten^2) / E
    F(d_T, beta) = (1 - ln(1 + (e^beta - 1) * d_T) / beta) / (1 - d_T),   d_T = W_g^s / W_g

with cyclic strain work ``W_g = (sigma_m^2 + sigma_a^2) / E`` for ``sigma_a >= |sigma_m|``
and ``2 |sigma_m| sigma_a / E`` otherwise (Eq. 29), rotating-bending annulus integration,
and smooth-cylinder combined tension-torsion limits (Eqs. 20-28, phase shift invariant).

Repository choices:
- Provided only as closed-form point/axisymmetric kernels rather than a registered model
  because general multiaxial stress states require integrating ``W_g`` over ``V*`` on an
  FE mesh.
- ``tension_torsion_index`` enforces ``sigma_a >= |sigma_m|`` and ``tau_a >= |tau_m|``, and
  returns the homogeneous reciprocal scale ``1 / lambda`` so the index varies continuously
  across the ``r* = 0`` transition (Eqs. 27-28).
- Retains ``d_T,uniax = (1 - 2 * nu) / 3`` from Eq. (22) (fixing the omitted ``1/3`` in
  Eq. 13 and inverted ``(y*/R)^2`` print in Eq. 10).
"""

from __future__ import annotations

import numpy as np


def banvillet_thresholds(s_ten: float, s_rb: float, E: float) -> tuple[float, float, float]:
    """Reference stress sigma*, threshold work W* (Eq. 8) and uniaxial limit g^D_uniax (Eq. 9)."""
    sigma_star = float(np.sqrt(2.0 * s_ten**2 - s_rb**2))
    return (
        sigma_star,
        sigma_star**2 / E,
        (s_rb**2 - s_ten**2) / E,
    )


def uniaxial_work(sigma_m: float, sigma_a: float, E: float) -> float:
    """Cyclic W_g for a sinusoidal uniaxial stress, Eq. (29)."""
    if sigma_a >= abs(sigma_m):
        return (sigma_m**2 + sigma_a**2) / E
    return 2.0 * abs(sigma_m) * sigma_a / E


def uniaxial_limit_amplitude(sigma_m: float, s_ten: float) -> float:
    """Endurance amplitude at prescribed mean stress for homogeneous uniaxial loading (W_g = s_Ten^2 / E)."""
    if abs(sigma_m) <= s_ten / np.sqrt(2.0):
        return float(np.sqrt(s_ten**2 - sigma_m**2))
    return float(s_ten**2 / (2.0 * abs(sigma_m)))


def rotating_bending_gbar(sigma_a: float, sigma_star: float, E: float) -> float:
    """Volume average of W_g - W* over V* for fully reversed rotating bending of a smooth cylinder."""
    return max(sigma_a**2 - sigma_star**2, 0.0) / (2.0 * E)


def triaxiality_function(d_T: float, beta: float) -> float:
    """F(d_T, beta) of Eq. (14); F(0) = 1 (fully reversed torsion)."""
    log_term = np.log1p(np.expm1(beta) * d_T) / beta if beta != 0.0 else d_T
    return float((1.0 - log_term) / (1.0 - d_T))


def beta_from_limits(s_rb: float, t_m1: float, nu: float) -> float:
    """beta obtained from Eq. (15) via the rotating-bending and torsion limits (requires s_RB/t_-1 < sqrt(3))."""
    target = (s_rb / t_m1) ** 2 / 3.0
    if not 0.0 < target < 1.0:
        raise ValueError("Eq. (15) admits a root only for 0 < s_RB / t_-1 < sqrt(3)")
    d_T_uniax = (1.0 - 2.0 * nu) / 3.0
    lower, upper = -60.0, 60.0  # 1 - F(d_uni)(1 - d_uni) grows monotonically with beta
    for _ in range(200):
        midpoint = 0.5 * (lower + upper)
        if triaxiality_function(d_T_uniax, midpoint) * (1.0 - d_T_uniax) > target:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def tension_torsion_index(
    sigma_m: float,
    sigma_a: float,
    tau_m: float,
    tau_a: float,
    s_ten: float,
    s_rb: float,
    beta: float,
    nu: float,
) -> float:
    """Reciprocal of lambda, where lambda scales the entire stress state of a solid smooth cylinder until
    the endurance limit is attained, Eq. (27) (r* > 0) or Eq. (28) (r* = 0).

    Both d_T and F are invariant under that scaling, and the regime of Eq. (26) is the one compatible
    with the limit load, so the index varies continuously and is homogeneous of degree 1 in the stresses
    (equal to 1 at the limit). Forming the ratio of the two sides of Eq. (27) or (28) at the applied load
    would instead jump where r* reaches 0.
    """
    if sigma_a < abs(sigma_m) or tau_a < abs(tau_m):
        raise ValueError("Eqs. (21)-(28) require sigma_a >= |sigma_m| and tau_a >= |tau_m|")
    normal_part = sigma_m**2 + sigma_a**2
    shear_part = 2.0 * (1.0 + nu) * (tau_m**2 + tau_a**2)
    if normal_part + shear_part == 0.0:
        return 0.0
    d_T_uniax = (1.0 - 2.0 * nu) / 3.0
    scale_factor = triaxiality_function(
        d_T_uniax * normal_part / (normal_part + shear_part), beta
    ) / triaxiality_function(d_T_uniax, beta)
    ref_work = scale_factor * (2.0 * s_ten**2 - s_rb**2)  # E W* of Eq. (24)
    lambda_sq = scale_factor * s_rb**2 / (shear_part + normal_part)  # lambda^2 at Eq. (27)
    if lambda_sq * normal_part >= ref_work:  # r* = 0 at that load: Eq. (28)
        lambda_sq = scale_factor * s_ten**2 / (0.5 * shear_part + normal_part)
    return float(1.0 / np.sqrt(lambda_sq))
