"""Liu & Mahadevan (2005) multiaxial high-cycle fatigue criterion and life model (N16).

Y. Liu, S. Mahadevan, Int. J. Fatigue 27 (2005) 790-800, doi:10.1016/j.ijfatigue.2005.01.003,
Eqs. (6), (19)-(24), Table 1::

    sqrt( [sigma_a,c (1 + eta sigma_m,c / f_-1) / f_-1]^2 + (tau_a,c / t_-1)^2
          + k (sigma^H_a / f_-1)^2 ) = beta                                               (20)
    s = t_-1 / f_-1 <= 1:  cos 2 alpha = [-2 + sqrt(4 - 4 (1/s^2 - 3)(5 - 1/s^2 - 4 s^2))]
                                         / [2 (5 - 1/s^2 - 4 s^2)],
                           beta = [cos^2(2 alpha) s^2 + sin^2(2 alpha)]^(1/2),  k = 0      (19)
    s > 1:                 alpha = 0,  beta = s,  k = 9 (s^2 - 1)                          (Table 1)
    eta = 3/4 + (sqrt(3) - f_-1/t_-1) / (4 (sqrt(3) - 1))  (s <= 1),  eta = 1 (s > 1)      (21)
    (1/beta) sqrt( [sigma_a,c (1 + eta_N sigma_m,c / f_N)]^2 + (f_N/t_N)^2 tau_a,c^2
                   + k (sigma^H_a)^2 ) = f_N,   s = t_N / f_N                             (23)-(24)

The fracture plane is the plane of maximum normal stress amplitude; the critical plane lies at
the material angle alpha from it; sigma^H_a is the hydrostatic stress amplitude. Eq. (19) makes
both fully reversed limits (bending f_-1, torsion t_-1) satisfy Eq. (20) exactly. The life follows
from Eq. (23) with the uniaxial and torsional S-N curves f_N, t_N, solved by iteration; alpha, beta,
k and eta are re-evaluated at s = t_N / f_N (Section 3.4).

Repository choices:
  - maintainer decision: the paper's rotation is 2D; in 3D the critical plane is the worst plane on
    the cone of planes at angle alpha around the fracture-plane normal (checked: both reversed limits
    still give exactly beta for t_-1/f_-1 >= 0.5; below, a tilted plane of the cone exceeds the
    torsion limit, so s < 0.5 (under the Tresca ratio, not met by metals) raises ValueError).
    Several fracture planes of equal amplitude: the worst one is kept;
  - maintainer decision: the shear amplitude tau_a,c is the rss of the harmonic shear components
    (the paper is silent); it reproduces Fig. 2 (hard steel) at 0, 60 and 90 deg, the longest chord
    does not. Figs. 2-5 are not reproduced for 42CrMo4, 34Cr4 and 30NCD16 (unresolved);
  - Eq. (19) is evaluated in its rationalised form cos 2 alpha = -C / (1 + sqrt(1 - A C)),
    A = 5 - 1/s^2 - 4 s^2, C = 1/s^2 - 3, which removes the 0/0 at s = 1/2 and s = 1;
  - eta from Eq. (21) (the simplified function; experimental eta values are not used); the
    bracket (1 + eta sigma_m,c / f) is clipped at zero for large compressive means;
  - maintainer decision: native life of Eq. (23) with Basquin curves f_N = sigma'_f (2N)^b and
    t_N = tau'_f (2N)^b0 (the paper allows any S-N form; Table 5 uses other forms), solved by
    bisection in log N; DP is the fatigue-limit index sqrt(...)/beta of Eq. (20), a diagnostic;
  - the life is searched in 1-1e12 cycles, restricted to the lives where the S-N curves give
    t_N/f_N >= 0.5; a solution outside raises LoadingScopeError, none below 1e12 gives an
    infinite life;
  - the cone is sampled at CONE_POINTS azimuths and the fracture plane found on a Fibonacci grid
    refined by a shrinking local search.

User guide: docs/method-guides/LIUMAHADEVAN.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..comparison import ComparisonCapabilities
from ..errors import LoadingScopeError, MissingParamsError
from ..invariant_methods.stress_invariants import harmonic_parts

CONE_POINTS = 720
GRID_POINTS = 1500
TIE_TOL = 1e-7
S_MIN = 0.5  # validity of the alpha-cone rule (see the module docstring)
LOG_N_MAX = 12.0  # beyond 1e12 cycles the life is reported as infinite
_LIFE_KEYS = ("sigma_f", "b_fat", "tau_f", "b0")


def material_angle(s: float) -> tuple[float, float, float]:
    """alpha (rad), beta and k of Eq. (19) / Table 1 for s = t_-1 / f_-1."""
    if not s >= S_MIN:
        raise ValueError(
            f"LIUMAHADEVAN: t_-1/f_-1 = {s:.3f} < {S_MIN}; below the Tresca ratio the worst plane of "
            "the alpha-cone no longer meets the torsion limit (repository restriction)"
        )
    if s >= 1.0:
        return 0.0, float(s), 9.0 * (s * s - 1.0)
    A = 5.0 - 1.0 / s**2 - 4.0 * s**2
    C = 1.0 / s**2 - 3.0
    c2a = float(np.clip(-C / (1.0 + np.sqrt(1.0 - A * C)), -1.0, 1.0))
    beta = float(np.sqrt(c2a**2 * s**2 + 1.0 - c2a**2))
    return 0.5 * float(np.arccos(c2a)), beta, 0.0


def mean_stress_factor(s: float) -> float:
    """eta of Eq. (21) (simplified function) for s = t / f."""
    if s >= 1.0:
        return 1.0
    return 0.75 + 0.25 * (np.sqrt(3.0) - 1.0 / s) / (np.sqrt(3.0) - 1.0)


def _hemisphere(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    z = 1.0 - i / n
    phi = np.pi * (1.0 + 5.0**0.5) * i
    r = np.sqrt(1.0 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def _normal_amplitude(N, S_sin, S_cos):
    a = np.einsum("ni,ij,nj->n", N, S_sin, N)
    b = np.einsum("ni,ij,nj->n", N, S_cos, N)
    return np.hypot(a, b)


def _tangents(n):
    t = np.cross(n, [1.0, 0.0, 0.0])
    if np.linalg.norm(t) < 1e-6:
        t = np.cross(n, [0.0, 1.0, 0.0])
    t /= np.linalg.norm(t)
    return t, np.cross(n, t)


def fracture_normals(S_sin, S_cos) -> list[np.ndarray]:
    """Normals of the planes of maximum normal stress amplitude (all ties)."""
    grid = _hemisphere(GRID_POINTS)
    amp = _normal_amplitude(grid, S_sin, S_cos)
    azimuths = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    refined = []
    for idx in np.argsort(-amp)[:8]:
        n, best, step = grid[idx], amp[idx], 0.05
        while step > 1e-8:
            t, u = _tangents(n)
            trial = n + step * (np.cos(azimuths)[:, None] * t + np.sin(azimuths)[:, None] * u)
            trial /= np.linalg.norm(trial, axis=1)[:, None]
            values = _normal_amplitude(trial, S_sin, S_cos)
            j = int(np.argmax(values))
            if values[j] > best:
                n, best = trial[j], values[j]
            else:
                step *= 0.5
        refined.append((best, n))
    top = max(a for a, _ in refined)
    out: list[np.ndarray] = []
    for a, n in refined:
        if a >= top * (1.0 - TIE_TOL) and all(abs(abs(float(n @ m)) - 1.0) > 1e-6 for m in out):
            out.append(n)
    return out


def _cone(n, alpha):
    if alpha == 0.0:
        return n[None]
    t, u = _tangents(n)
    p = np.linspace(0.0, 2.0 * np.pi, CONE_POINTS, endpoint=False)
    return np.cos(alpha) * n + np.sin(alpha) * (np.cos(p)[:, None] * t + np.sin(p)[:, None] * u)


def plane_stresses(N, S_mean, S_sin, S_cos):
    """sigma_a, sigma_m and rss tau_a on the planes with normals N."""
    t_sin, t_cos, t_mean = N @ S_sin.T, N @ S_cos.T, N @ S_mean.T
    n_sin = np.einsum("ni,ni->n", t_sin, N)
    n_cos = np.einsum("ni,ni->n", t_cos, N)
    n_mean = np.einsum("ni,ni->n", t_mean, N)
    shear_sin = t_sin - n_sin[:, None] * N
    shear_cos = t_cos - n_cos[:, None] * N
    tau_a = np.sqrt(np.sum(shear_sin**2, axis=1) + np.sum(shear_cos**2, axis=1))
    return np.hypot(n_sin, n_cos), n_mean, tau_a


class _Geometry:
    """Fracture normals and hydrostatic amplitude of one cycle, shared by all evaluations."""

    def __init__(self, S_mean, S_sin, S_cos):
        self.S = tuple(np.asarray(x, float) for x in (S_mean, S_sin, S_cos))
        self.normals = fracture_normals(self.S[1], self.S[2])
        self.sigma_H_a = float(np.hypot(np.trace(self.S[1]), np.trace(self.S[2])) / 3.0)

    def equivalent(self, f: float, t: float) -> float:
        """(1/beta) sqrt(...) of Eqs. (20)/(23) at strengths f, t (MPa)."""
        s = t / f
        alpha, beta, k = material_angle(s)
        eta = mean_stress_factor(s)
        worst = 0.0
        for n in self.normals:
            sigma_a, sigma_m, tau_a = plane_stresses(_cone(n, alpha), *self.S)
            bracket = np.maximum(1.0 + eta * sigma_m / f, 0.0)
            value = np.sqrt(
                (sigma_a * bracket) ** 2 + (f / t) ** 2 * tau_a**2 + k * self.sigma_H_a**2
            )
            worst = max(worst, float(np.max(value)))
        return worst / beta


def liu_mahadevan_index(S_mean, S_sin, S_cos, f_m1: float, t_m1: float) -> float:
    """Fatigue-limit index of Eq. (20): sqrt(...)/beta (1 on the limit)."""
    return _Geometry(S_mean, S_sin, S_cos).equivalent(f_m1, t_m1) / f_m1


def liu_mahadevan_life(
    S_mean, S_sin, S_cos, sigma_f, b_fat, tau_f, b0, *, geometry: _Geometry | None = None
) -> float:
    """N_f of Eq. (23) with f_N = sigma'_f (2N)^b, t_N = tau'_f (2N)^b0 (bisection in log10 N)."""
    geometry = geometry or _Geometry(S_mean, S_sin, S_cos)

    def excess(log_n: float) -> float:
        two_n = 2.0 * 10.0**log_n
        f_n, t_n = sigma_f * two_n**b_fat, tau_f * two_n**b0
        return geometry.equivalent(f_n, t_n) / f_n - 1.0

    lo, hi = _valid_log_lives(sigma_f, b_fat, tau_f, b0)
    if excess(hi) < 0.0:
        if hi < LOG_N_MAX:
            raise LoadingScopeError(_SN_SCOPE)
        return float("inf")
    if excess(lo) >= 0.0:
        if lo > 0.0:
            raise LoadingScopeError(_SN_SCOPE)
        return 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if excess(mid) < 0.0 else (lo, mid)
    return float(10.0 ** (0.5 * (lo + hi)))


_SN_SCOPE = (
    "LIUMAHADEVAN: the life lies where the S-N curves give t_N/f_N < 0.5 (outside the model's "
    "repository validity)"
)


def _valid_log_lives(sigma_f, b_fat, tau_f, b0) -> tuple[float, float]:
    """log10 N interval of [0, LOG_N_MAX] where s_N = t_N / f_N >= S_MIN (s_N is monotone)."""
    grid = np.linspace(0.0, LOG_N_MAX, 1201)
    s_n = (tau_f / sigma_f) * (2.0 * 10.0**grid) ** (b0 - b_fat)
    ok = grid[s_n >= S_MIN]
    if ok.size == 0:
        raise LoadingScopeError(_SN_SCOPE)
    return float(ok[0]), float(ok[-1])


class LiuMahadevan:
    """Liu & Mahadevan (2005) criterion with its native life law (see the module docstring)."""

    name = "LIUMAHADEVAN"
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_comparison=False, dp_reference="none"
    )

    def required_params(self) -> set[str]:
        return {"Sigm1", "Taum1", *_LIFE_KEYS}

    def _life_constants(self, params):
        values = []
        for key in _LIFE_KEYS:
            raw = params.get(key)
            if raw is None or not np.isfinite(float(raw)):
                raise MissingParamsError(self.name, [key], "axial and torsional Basquin curves")
            values.append(float(raw))
        sigma_f, b_fat, tau_f, b0 = values
        if sigma_f <= 0.0 or tau_f <= 0.0 or b_fat >= 0.0 or b0 >= 0.0:
            raise ValueError("LIUMAHADEVAN requires sigma_f, tau_f > 0 and b_fat, b0 < 0")
        return values

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        geometry = _Geometry(S_mean, S_sin, S_cos)
        f_m1 = float(params["Sigm1"])
        index = geometry.equivalent(f_m1, float(params["Taum1"])) / f_m1
        Nf = liu_mahadevan_life(
            S_mean, S_sin, S_cos, *self._life_constants(params), geometry=geometry
        )
        return CaseResult(values={"": index}, aux={"Nf": Nf})

    def life_fit(self, aux, N, params: Mapping[str, float]):
        """Nothing to fit: the S-N curves are the material's Basquin constants (returned as record)."""
        return {key: float(params[key]) for key in _LIFE_KEYS}

    def native_life(self, aux, params: Mapping[str, float], fit=None):
        """N_f of Eq. (23), computed in ``evaluate_case``."""
        return aux["Nf"]
