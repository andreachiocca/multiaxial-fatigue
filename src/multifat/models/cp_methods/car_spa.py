"""Carpinteri & Spagnoli (2001) critical-plane criterion for hard metals (CARSPA).

A. Carpinteri, A. Spagnoli, Int. J. Fatigue 23 (2001) 135-145, doi:10.1016/S0142-1123(00)00075-X,
Eqs. (22)-(27), with the Euler-angle averaging of A. Carpinteri, E. Macha, R. Brighenti,
A. Spagnoli, Int. J. Fatigue 21 (1999) 83-88, doi:10.1016/S0142-1123(98)00046-2, Eqs. (1)-(9),
Fig. 3 and Table 1, in the notation of A. Carpinteri, R. Brighenti, A. Spagnoli, Fatigue Fract.
Engng Mater. Struct. 23 (2000) 355-364, doi:10.1046/j.1460-2695.2000.00265.x::

    phi^ = sum W(t_k) phi(t_k) / W,  theta^, psi^ alike,  W = sum W(t_k)          (22)
    W(t_k) = 0 if sigma_1(t_k) < c sigma_af,  (sigma_1(t_k)/sigma_af)^m_sigma else  (23)
    delta = 45 (3/2) [1 - (tau_af/sigma_af)^2]   [deg], angle between 1^ and w     (24)
    sqrt(N_max^2 + (sigma_af/tau_af)^2 C_a^2) <= sigma_af,  N_max = N_m + N_a      (25)-(27)

phi, theta, psi are the principal Euler angles of the principal directions 1, 2, 3
(sigma_1 >= sigma_2 >= sigma_3) at time t_k, reduced by the two-stage procedure of the 1999
paper to 0 <= phi, theta <= 90 deg and -90 < psi <= 90 deg. 1^, 2^, 3^ are the directions of the
averaged angles. N_m, N_a are the mean and amplitude of the normal stress on the critical plane
of normal w, C_a the amplitude of the shear stress vector (major semi-axis of its ellipse,
2001 Eqs. 9-13). The equivalent stress sigma_eq (Eq. 27) is the returned damage parameter;
sigma_eq = sigma_af on the fatigue limit, and the error index of the paper is
I = 100 (sigma_eq - sigma_af) / sigma_af (Eq. 30).

Repository choices:
  - Euler angles: matrix (2) of the 1999 paper with the 2000/2001 notation, in which direction 1
    is the third column (1 = (cos phi sin theta, sin phi sin theta, cos theta)), 2 the first and
    3 the second. The 1st reduction (Fig. 3) is done by enumerating the column-sign changes with
    det = +1 and the two equivalent Euler triples, keeping the one with all angles in
    [-90, 90] deg; the 2nd reduction is Table 1 (phi, theta -> |phi|, |theta|; psi keeps its
    sign when phi psi >= 0, else changes it). psi = -90 deg is mapped to +90 deg (range
    -90 < psi <= 90 of the 2000 paper). At theta = 0 only phi + psi is defined; it is put in phi
    (psi = 0). These choices reproduce the explicit plane-stress cases (a1)-(c2) of the 2000 paper;
  - mirror and return (maintainer decision): the 2nd reduction folds directions into the first
    octant of the fixed frame, i.e. it mirrors them, so with the fixed frame alone a load and its
    mirror image (e.g. tau -> -tau) would give different planes. Before averaging, the fixed axes
    are mirrored so that direction 1 at the instant of largest weight lies in the first octant;
    the critical-plane normal is mirrored back afterwards. For the published tests (positive
    shear) the mirror is the identity and the tables are unchanged;
  - w = cos(delta) 1^ + sin(delta) 3^, a rotation in the plane 1^3^. The 2001 paper gives only
    the angle; this direction reproduces Tables 2-4, the opposite one gives errors down to -15 %.
    It is the clockwise rotation about 2^ stated by Marques, Benasciutti, Carpinteri & Spagnoli
    (2020), Fatigue Fract Engng Mater Struct 43:1978-1993, doi:10.1111/ffe.13273, Section 3.1, for
    the spectral CSV version (same frame: 1^ = Z', 2^ = X', 3^ = Y');
  - c = 0.5 as in Tables 2-4 (``c_CARSPA``). m_sigma = -1/m is the inverse slope of the R = -1
    tension S-N curve (``m_sigma``); without it, -1/b_fat of the Basquin fit is used. It only
    matters when the weighted instants have different principal directions: if neither is
    available such a case raises MissingParamsError (the model is skipped on that material);
  - W uses the cycle sampled at N_TIME instants of one harmonic period. Where two principal
    stresses are equal (e.g. sigma_2 = sigma_3 = 0 at the peak of an out-of-phase bending with
    tau = 0) the directions are those of the right limit in time. When sigma_1 never
    reaches c sigma_af (loads far below the limit, e.g. in a load-scaling search) the
    published average is undefined; the instants of largest sigma_1 are used with equal weight,
    which is the limit of Eq. (23) as the window shrinks to the peak;
  - scope: tau_af/sigma_af in [1/sqrt(3), 1] (hard metals, delta in [0, 45] deg), otherwise
    MaterialScopeError; a compressive mean normal stress N_m < -N_a on the critical plane is
    excluded by the source and raises LoadingScopeError. The mean shear stress is not used.
  - sigma_eq is the damage parameter of the shared master curve DP = A N^b for finite-life
    comparisons (not part of the criterion).

User guide: docs/method-guides/CARSPA.md.
"""

from __future__ import annotations

import itertools
from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError, MaterialScopeError, MissingParamsError
from ..invariant_methods.stress_invariants import harmonic_parts

N_TIME = 720  # instants per period used for the weighted average (Eq. 22)
C_DEFAULT = 0.5  # c of Eq. (23), the value used for Tables 2-4
ANGLE_TOL = 1e-9  # rad, tolerance of the range checks of the two reductions
SCOPE_TOL = 1e-9  # relative tolerance of the scope checks
TIE_TOL = 1e-9  # relative gap below which two principal stresses are equal
TIE_STEP = 1e-6  # rad, time offset of the right limit taken at equal principal stresses
HALF_PI = 0.5 * np.pi


def euler_matrix(phi, theta, psi) -> np.ndarray:
    """Matrix (2) of Carpinteri et al. (1999): rotations phi, theta, psi about Z, Y', 3."""
    cf, sf = np.cos(phi), np.sin(phi)
    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(psi), np.sin(psi)
    return np.array(
        [
            [cf * ct * cp - sf * sp, -cf * ct * sp - sf * cp, cf * st],
            [sf * ct * cp + cf * sp, -sf * ct * sp + cf * cp, sf * st],
            [-st * cp, st * sp, ct],
        ]
    )


def principal_axes(e1, e2, e3) -> np.ndarray:
    """Matrix (2) slots of the principal directions 1, 2, 3 in the 2000/2001 notation."""
    return np.stack([e2, e3, e1], axis=-1)


def _zyz(M: np.ndarray):
    """Euler angles (phi, theta, psi) of matrix (2) for a stack of rotation matrices."""
    theta = np.arccos(np.clip(M[..., 2, 2], -1.0, 1.0))
    lock = np.sin(theta) < 1e-12
    phi = np.where(
        lock, np.arctan2(M[..., 1, 0], M[..., 0, 0]), np.arctan2(M[..., 1, 2], M[..., 0, 2])
    )
    psi = np.where(lock, 0.0, np.arctan2(M[..., 2, 1], -M[..., 2, 0]))
    return phi, theta, psi


def _wrap(a):
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def reduced_euler_angles(axes: np.ndarray) -> np.ndarray:
    """Reduced principal Euler angles of Carpinteri et al. (1999), both stages.

    Parameters
    ----------
    axes : (..., 3, 3) array
        Columns: principal directions 1, 2, 3 (sigma_1 >= sigma_2 >= sigma_3), any signs.

    Returns
    -------
    (..., 3) array of (phi, theta, psi) in rad, 0 <= phi, theta <= pi/2, -pi/2 < psi <= pi/2.
    """
    axes = np.asarray(axes, dtype=float)
    M0 = principal_axes(axes[..., 0], axes[..., 1], axes[..., 2])
    out = np.full(M0.shape[:-2] + (3,), np.nan)
    found = np.zeros(M0.shape[:-2], dtype=bool)
    # 1st reduction: column-sign changes with det = +1, and the two equivalent Euler triples
    for signs in itertools.product((1.0, -1.0), repeat=3):
        M = M0 * np.asarray(signs)
        proper = np.linalg.det(M) > 0.0
        f, t, p = _zyz(M)
        for a, b, c in ((f, t, p), (f + np.pi, -t, p + np.pi)):
            a, b, c = _wrap(a), _wrap(b), _wrap(c)
            ok = (
                proper
                & ~found
                & (np.abs(a) <= HALF_PI + ANGLE_TOL)
                & (np.abs(b) <= HALF_PI + ANGLE_TOL)
                & (np.abs(c) <= HALF_PI + ANGLE_TOL)
            )
            out[ok] = np.stack([a, b, c], axis=-1)[ok]
            found |= ok
    if not np.all(found):
        raise RuntimeError("CARSPA: no reduced Euler triple found (non-orthogonal axes?)")
    phi, theta, psi = out[..., 0], out[..., 1], out[..., 2]
    # 2nd reduction (Table 1): cases (1)-(4) -> (+, +, +), cases (5)-(8) -> (+, +, -)
    psi_r = np.where(phi * psi >= 0.0, np.abs(psi), -np.abs(psi))
    psi_r = np.where(psi_r <= -HALF_PI + ANGLE_TOL, HALF_PI, psi_r)
    return np.stack([np.abs(phi), np.abs(theta), psi_r], axis=-1)


def carspa_delta(f_m1: float, t_m1: float) -> float:
    """Angle delta (rad) between 1^ and the critical-plane normal, Eq. (24)."""
    return np.radians(45.0 * 1.5 * (1.0 - (float(t_m1) / float(f_m1)) ** 2))


def check_material_scope(f_m1: float, t_m1: float) -> None:
    """Raise MaterialScopeError outside 1/sqrt(3) <= tau_af/sigma_af <= 1 (hard metals)."""
    if not (f_m1 > 0.0 and t_m1 > 0.0):
        raise ValueError(f"CARSPA needs sigma_af > 0 and tau_af > 0, got {f_m1}, {t_m1}")
    r = t_m1 / f_m1
    if r < 1.0 / np.sqrt(3.0) - SCOPE_TOL or r > 1.0 + SCOPE_TOL:
        raise MaterialScopeError(
            f"CARSPA: tau_af/sigma_af = {r:.3f} outside [1/sqrt(3), 1], the hard-metal range "
            "of Carpinteri & Spagnoli (2001) Eq. (24)"
        )


def _principal(S_mean, S_sin, S_cos, t):
    """Principal stresses (descending) and directions (columns) at the instants t."""
    S = (
        np.asarray(S_mean, float)[None]
        + np.sin(t)[:, None, None] * np.asarray(S_sin, float)[None]
        + np.cos(t)[:, None, None] * np.asarray(S_cos, float)[None]
    )
    vals, vecs = np.linalg.eigh(S)
    order = np.argsort(-vals, axis=1, kind="stable")
    return np.take_along_axis(vals, order, axis=1), np.take_along_axis(
        vecs, order[:, None, :], axis=2
    )


def averaged_directions(
    S_mean, S_sin, S_cos, f_m1: float, m_sigma: float | None, c: float = C_DEFAULT
) -> np.ndarray:
    """Weighted mean principal directions 1^, 2^, 3^ (columns), Eqs. (22)-(23).

    Includes the mirror-and-return choice of the module docstring. m_sigma may be None when the
    weighted instants share one principal frame; otherwise MissingParamsError is raised.
    """
    t = 2.0 * np.pi * np.arange(N_TIME) / N_TIME
    vals, vecs = _principal(S_mean, S_sin, S_cos, t)
    # equal principal stresses leave the directions undefined: take the right limit in time
    scale = max(float(np.max(np.abs(vals))), 1.0)
    tied = np.min(np.abs(np.diff(vals, axis=1)), axis=1) <= TIE_TOL * scale
    if np.any(tied):
        _, vecs[tied] = _principal(S_mean, S_sin, S_cos, t[tied] + TIE_STEP)
    s1 = vals[:, 0]
    active = s1 >= c * f_m1
    if not np.any(active):
        # Eq. (23) undefined: the window shrunk to the peak (repository choice)
        active = s1 >= s1.max() - 1e-9 * max(abs(s1.max()), 1.0)
    vecs, s1 = vecs[active], s1[active]
    # mirror the fixed axes: direction 1 at the instant of largest sigma_1 in the first octant
    # (ties between instants are broken by the direction itself, not by the time order)
    peak = s1 >= s1.max() - 1e-9 * max(abs(s1.max()), 1.0)
    leads = vecs[peak, :, 0]
    big = np.take_along_axis(leads, np.argmax(np.abs(leads), axis=1)[:, None], axis=1)
    leads = np.round(leads * np.where(big >= 0.0, 1.0, -1.0), 9)
    lead = leads[np.lexsort(leads.T[::-1])[-1]]
    mirror = np.where(lead < 0.0, -1.0, 1.0)
    angles = reduced_euler_angles(mirror[None, :, None] * vecs)
    if np.ptp(angles, axis=0).max() <= 1e-9:
        mean = angles[0]
    else:
        if m_sigma is None:
            raise MissingParamsError(
                "CARSPA",
                {"m_sigma"},
                "the principal directions rotate within the weighted window, so W(t_k) of "
                "Eq. (23) needs m_sigma (or b_fat)",
            )
        w = (s1 / s1.max()) ** float(m_sigma)  # (s1/sigma_af)^m up to a constant factor
        mean = (w[:, None] * angles).sum(axis=0) / w.sum()
    A = euler_matrix(*mean)
    e1, e2, e3 = A[:, 2], A[:, 0], A[:, 1]
    return mirror[:, None] * np.stack([e1, e2, e3], axis=-1)


def plane_stresses(S_mean, S_sin, S_cos, w) -> tuple[float, float, float]:
    """(N_m, N_a, C_a) on the plane of normal w for one harmonic (2001 Eqs. 6-13)."""
    w = np.asarray(w, float)
    Sm, Ss, Sc = (np.asarray(T, float) for T in (S_mean, S_sin, S_cos))
    n_m = float(w @ Sm @ w)
    n_a = float(np.hypot(w @ Ss @ w, w @ Sc @ w))
    P = np.eye(3) - np.outer(w, w)
    cs, cc = P @ Ss @ w, P @ Sc @ w
    gram = np.array([[cs @ cs, cs @ cc], [cs @ cc, cc @ cc]])
    c_a = float(np.sqrt(max(np.linalg.eigvalsh(gram).max(), 0.0)))  # major semi-axis
    return n_m, n_a, c_a


def carspa_critical_plane(S_mean, S_sin, S_cos, params: Mapping[str, float]) -> dict:
    """Critical-plane normal and stresses; keys: w, axes, delta, N_m, N_a, C_a, sigma_eq."""
    f_m1, t_m1 = float(params["Sigm1"]), float(params["Taum1"])
    check_material_scope(f_m1, t_m1)
    m_sigma = params.get("m_sigma")
    if m_sigma is None and params.get("b_fat") is not None:
        m_sigma = -1.0 / float(params["b_fat"])
    c = float(params.get("c_CARSPA") if params.get("c_CARSPA") is not None else C_DEFAULT)
    axes = averaged_directions(S_mean, S_sin, S_cos, f_m1, m_sigma, c)
    delta = carspa_delta(f_m1, t_m1)
    w = np.cos(delta) * axes[:, 0] + np.sin(delta) * axes[:, 2]
    n_m, n_a, c_a = plane_stresses(S_mean, S_sin, S_cos, w)
    scale = max(float(np.max(np.abs(S_sin))), float(np.max(np.abs(S_cos))), 1.0)
    if n_m < -n_a - SCOPE_TOL * scale:
        raise LoadingScopeError(
            "CARSPA: compressive mean normal stress N_m < -N_a on the critical plane, excluded "
            "by Carpinteri & Spagnoli (2001) below Eq. (25)"
        )
    sigma_eq = float(np.sqrt((n_m + n_a) ** 2 + (f_m1 / t_m1) ** 2 * c_a**2))
    return dict(w=w, axes=axes, delta=delta, N_m=n_m, N_a=n_a, C_a=c_a, sigma_eq=sigma_eq)


class CarSpa:
    name = "CARSPA"
    optional_params = {
        "m_sigma": "-1/b_fat; needed only when the principal directions rotate within the window",
        "b_fat": "only to estimate m_sigma",
        "c_CARSPA": "0.5 (Tables 2-4 of the source)",
    }

    def required_params(self) -> set[str]:
        # Sigm1: sigma_af, fully reversed axial limit; Taum1: tau_af, fully reversed torsion limit
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
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        res = carspa_critical_plane(S_mean, S_sin, S_cos, params)
        return CaseResult(values={"": res["sigma_eq"]})
