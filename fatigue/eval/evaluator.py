# fatigue/eval/evaluator.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping, Tuple
import numpy as np
import math
from fatigue.models.base import CriticalPlaneMethod
from fatigue.models.utils.plane_history import shear_maximum, split_strain_ranges

def cp_rot(
    S0: np.ndarray,
    E0: np.ndarray,
    S1: np.ndarray,
    E1: np.ndarray,
    R: np.ndarray,
    *,
    check_orthonormal: bool = True,
    atol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rotate stress/strain tensors into the rotated frame defined by R.

    Same passive change-of-basis as your old `cp_fun.cp_rot`:
        T_rot = R.T @ T @ R
    """

    def _as33(name: str, A: np.ndarray) -> np.ndarray:
        A = np.asarray(A, dtype=float)
        if A.shape != (3, 3):
            raise ValueError(f"{name} must be shape (3,3); got {A.shape}.")
        return A

    S0 = _as33("S0", S0)
    E0 = _as33("E0", E0)
    S1 = _as33("S1", S1)
    E1 = _as33("E1", E1)
    R = _as33("R", R)

    if check_orthonormal:
        I = np.eye(3)
        RtR = R.T @ R
        if not np.allclose(RtR, I, atol=atol):
            raise ValueError("R is not orthonormal within tolerance (R.T @ R != I).")

    def rot(T: np.ndarray) -> np.ndarray:
        return R.T @ T @ R

    return rot(S0), rot(E0), rot(S1), rot(E1)


@dataclass(frozen=True)
class BestCP:
    ext: float        # max damage over all planes
    by_metric: float  # damage on plane that maximizes metric


def best_cp_for_case(
    *,
    S0: np.ndarray,
    E0: np.ndarray,
    S1: np.ndarray,
    E1: np.ndarray,
    R_list: np.ndarray,
    method: CriticalPlaneMethod,
    params: Mapping[str, float],
    S_sin: np.ndarray | None = None,
    S_cos: np.ndarray | None = None,
    E_sin: np.ndarray | None = None,
    E_cos: np.ndarray | None = None,
) -> BestCP:
    best_ext = -np.inf
    best_metric = -np.inf
    best_by_metric = np.nan

    if len(R_list) == 0:
        raise ValueError("Critical-plane evaluation requires at least one plane")

    for Ri in R_list:
        S0r, E0r, S1r, E1r = cp_rot(S0, E0, S1, E1, Ri)

        # ------------------------------------------------------------
        # Out-of-phase support (single-harmonic) for CP methods
        # ------------------------------------------------------------
        # If harmonic coefficients are available, we override ONLY the
        # shear traction / shear strain components on the candidate plane
        # so that methods that rely on shear RANGE computed as:
        #   Δτ = sqrt((ΔSxz)^2 + (ΔSyz)^2)
        # automatically get an ellipse-aware equivalent range.
        # Normal scalar extrema are also reconstructed from their harmonics.
        if S_sin is not None and S_cos is not None:
            S_sin_r = Ri.T @ np.asarray(S_sin, dtype=float) @ Ri
            S_cos_r = Ri.T @ np.asarray(S_cos, dtype=float) @ Ri
            d_tau_eq = _equivalent_shear_range_from_harmonics(S_sin_r, S_cos_r)
            S0r, S1r = _override_plane_shear_components(
                S0r, S1r,
                delta_eq=d_tau_eq,
                T_sin_r=S_sin_r,
                T_cos_r=S_cos_r,
            )
            # Also update plane-normal (scalar) stress extrema (no ellipse needed)
            S0r, S1r = _override_plane_normal_scalar(
                S0r, S1r,
                T_sin_r=S_sin_r,
                T_cos_r=S_cos_r,
                idx=(2, 2),
            )

        if E_sin is not None and E_cos is not None:
            E_sin_r = Ri.T @ np.asarray(E_sin, dtype=float) @ Ri
            E_cos_r = Ri.T @ np.asarray(E_cos, dtype=float) @ Ri
            d_gam_eq = _equivalent_shear_range_from_harmonics(E_sin_r, E_cos_r)
            E0r, E1r = _override_plane_shear_components(
                E0r, E1r,
                delta_eq=d_gam_eq,
                T_sin_r=E_sin_r,
                T_cos_r=E_cos_r,
            )
            # Also update plane-normal (scalar) strain extrema (no ellipse needed)
            E0r, E1r = _override_plane_normal_scalar(
                E0r, E1r,
                T_sin_r=E_sin_r,
                T_cos_r=E_cos_r,
                idx=(2, 2),
            )

        # Only methods needing full-path quantities consume this internal data.
        # Never infer shear maxima or elastic/plastic strains from independently
        # synthesized normal/shear extrema: these do not form a physical tensor.
        plane_params = dict(params)
        if method.name in {"MGSE_YU", "MGSE_ZHU", "GSE", "GSA"} and S_sin is not None and S_cos is not None:
            mean_r = Ri.T @ (0.5 * (np.asarray(S0) + np.asarray(S1))) @ Ri
            plane_params["_tau_max"] = shear_maximum(mean_r, S_sin_r, S_cos_r)
        if method.name in {"GSE", "GSA"} and all(x is not None for x in (S_sin, S_cos, E_sin, E_cos)):
            plane_params["_strain_ranges"] = split_strain_ranges(
                S_sin_r, S_cos_r, E_sin_r, E_cos_r, float(params["E"]), float(params["nu"]))

        r = method.evaluate_on_plane(
            S0r=S0r, E0r=E0r,
            S1r=S1r, E1r=E1r,
            params=plane_params,
        )

        if not (math.isfinite(r.damage) and math.isfinite(r.metric)):
            raise ValueError(f"{method.name}: nonfinite plane damage or selection metric")

        if r.damage > best_ext:
            best_ext = r.damage

        if (r.metric > best_metric and not math.isclose(r.metric, best_metric, rel_tol=1e-10, abs_tol=1e-14)) or (
            math.isclose(r.metric, best_metric, rel_tol=1e-10, abs_tol=1e-14) and r.damage > best_by_metric
        ):
            best_metric = r.metric
            best_by_metric = r.damage

    return BestCP(ext=float(best_ext), by_metric=float(best_by_metric))


def _equivalent_shear_range_from_harmonics(
    T_sin_r: np.ndarray,
    T_cos_r: np.ndarray,
    *,
    a: Tuple[int, int] = (0, 2),
    b: Tuple[int, int] = (1, 2),
) -> float:
    """Equivalent shear RANGE on the candidate plane from single-harmonic coeffs.

    After rotation, the shear traction vector in the plane has two components
    (local xz and yz). With out-of-phase loading these two components trace
    an ellipse over one cycle. For single-harmonic loading, that ellipse can be
    written as:

        τ(t) = v_sin*sin(ωt) + v_cos*cos(ωt)

    where v_sin and v_cos are 2D vectors. The ellipse semi-axes lengths are
    sqrt(eigenvalues of M M^T) with M=[v_sin v_cos]. The user-requested shear
    amplitude is the vector-sum of major and minor semi-axes:

        τ_a,eq = sqrt(a^2 + b^2) = sqrt(||v_sin||^2 + ||v_cos||^2)

    This function returns the corresponding RANGE:

        Δτ_eq = 2 * τ_a,eq
    """
    T_sin_r = np.asarray(T_sin_r, dtype=float)
    T_cos_r = np.asarray(T_cos_r, dtype=float)

    v_sin = np.array([float(T_sin_r[a]), float(T_sin_r[b])], dtype=float)
    v_cos = np.array([float(T_cos_r[a]), float(T_cos_r[b])], dtype=float)

    amp = math.sqrt(float(v_sin @ v_sin + v_cos @ v_cos))
    return float(2.0 * amp)


def _override_plane_shear_components(
    T0r: np.ndarray,
    T1r: np.ndarray,
    *,
    delta_eq: float,
    T_sin_r: np.ndarray,
    T_cos_r: np.ndarray,
    a: Tuple[int, int] = (0, 2),
    b: Tuple[int, int] = (1, 2),
    eps: float = 1e-14,
) -> Tuple[np.ndarray, np.ndarray]:
    """Override the (xz,yz) components so that their 2-norm range equals delta_eq.

    This keeps the 2-loadstep interface for CP methods intact.
    """
    T0r = np.asarray(T0r, dtype=float).copy()
    T1r = np.asarray(T1r, dtype=float).copy()

    if not np.isfinite(delta_eq):
        raise ValueError("Nonfinite harmonic shear range")
    if abs(delta_eq) <= eps:
        # A zero physical range can coexist with nonzero legacy component-wise
        # endpoints (e.g. cancelling antiphase loads). Encode the mean exactly.
        for idx in (a, b):
            mean = 0.5 * (T0r[idx] + T1r[idx])
            T0r[idx] = T1r[idx] = mean
            T0r[idx[::-1]] = T1r[idx[::-1]] = mean
        return T0r, T1r

    v_sin = np.array([float(T_sin_r[a]), float(T_sin_r[b])], dtype=float)
    v_cos = np.array([float(T_cos_r[a]), float(T_cos_r[b])], dtype=float)

    # Pick a stable direction for the shear-range vector.
    w = v_sin
    nw = float(np.linalg.norm(w))
    if nw <= eps:
        w = v_cos
        nw = float(np.linalg.norm(w))
    if nw <= eps:
        # shear is essentially zero
        return T0r, T1r

    u = w / nw
    half = 0.5 * float(delta_eq) * u

    # Preserve the mean shear level already present in T0r/T1r.
    m_a = 0.5 * float(T0r[a] + T1r[a])
    m_b = 0.5 * float(T0r[b] + T1r[b])

    T0r[a] = m_a + half[0]
    T1r[a] = m_a - half[0]
    T0r[b] = m_b + half[1]
    T1r[b] = m_b - half[1]

    # enforce symmetry
    T0r[a[1], a[0]] = T0r[a]
    T1r[a[1], a[0]] = T1r[a]
    T0r[b[1], b[0]] = T0r[b]
    T1r[b[1], b[0]] = T1r[b]

    return T0r, T1r

def _override_plane_normal_scalar(
    T0r: np.ndarray,
    T1r: np.ndarray,
    *,
    T_sin_r: np.ndarray,
    T_cos_r: np.ndarray,
    idx: tuple[int, int] = (2, 2),
) -> Tuple[np.ndarray, np.ndarray]:
    """Override the normal (scalar) component on the plane using harmonic extrema.

    For any scalar quantity:
        q(t) = q_mean + q_sin*sin(ωt) + q_cos*cos(ωt)
    the amplitude is:
        q_a = sqrt(q_sin^2 + q_cos^2)
    and extrema are q_mean ± q_a.

    Here we treat the plane-normal component as T[idx] (default zz in the plane frame).
    """
    T0r = np.asarray(T0r, dtype=float).copy()
    T1r = np.asarray(T1r, dtype=float).copy()

    i, j = idx
    q_mean = 0.5 * float(T0r[i, j] + T1r[i, j])
    q_sin = float(np.asarray(T_sin_r, dtype=float)[i, j])
    q_cos = float(np.asarray(T_cos_r, dtype=float)[i, j])
    q_amp = math.sqrt(q_sin * q_sin + q_cos * q_cos)

    T0r[i, j] = q_mean + q_amp
    T1r[i, j] = q_mean - q_amp
    return T0r, T1r
