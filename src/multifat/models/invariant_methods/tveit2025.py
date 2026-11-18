"""Continuum fatigue model of Tveit, Ottosen, Erduran & Reyes (2025), Mech. Mater. 209, 105452.

Source: Tveit, Ottosen, Erduran & Reyes, Continuum mechanics fatigue model for arbitrary ratios of
shear and normal stress fatigue strengths, Mechanics of Materials 209 (2025) 105452,
doi:10.1016/j.mechmat.2025.105452, Eqs. (1)-(12), p. 2-4, and Eq. (24), p. 5:

    beta = (sigma_e + A I1 - f_m1) / f_m1                                         (1)
    sigma_e = sqrt(3 J2) - lambda (1 - 27/4 J3^2 / J2^3),  J2, J3 of s - alpha    (2)-(5)
    lambda = sqrt(3) t_m1 - f_m1                                                   (24)
    d alpha = C (s - alpha) d beta                                                 (6)
    dD = K exp(L beta) d beta                                                      (7)
    alpha and D evolve only while beta > 0 and d beta > 0                          (8)
    d beta = (G : ds + A dI1) / (f_m1 + C G : (s - alpha)),  G = d sigma_e / ds    (9)-(12)

The back-stress alpha and the damage D start at zero; crack initiation occurs at D = 1, and the
number of load cycles needed to reach it is the native life. A path that the moving endurance
surface ends up enclosing accumulates no further damage and has an infinite life.

Repository choices:
  - The cycle S(t) = S_mean + S_sin sin(2 pi t) + S_cos cos(2 pi t) is sampled at 1001 points of
    one closed period and repeated as blocks (one block = one cycle). Increments between
    consecutive points are explicit, and the first point of every block is not integrated.
  - Blocks are repeated until one of these occurs:
      failure: D >= 1 inside a block. The life is interpolated linearly in the increment that
        crosses D = 1: block n fails after (n - 1) + (i - 1 + f)/(npts - 1) cycles when the
        increment from point i - 1 to point i crosses at its fraction f.
      shakedown: a block (from the second on) without any damage increment. The back-stress
        moves only together with the damage, so every later block is the same: infinite life.
      cyclic steady state: from block 400 on, the back-stress shows no net drift over the last
        200 blocks (distance to its position 200 blocks earlier <= 0.01 times its path length in
        those blocks; in a steady cycle it moves back and forth, in a shakedown one way only), and
        the mean damage per block of the last 200 blocks is within 0.1 % of the 200 blocks before.
        A back-stress that swings in every block can still creep one way, and the peak of beta
        over the block then sinks until a shakedown. So the steady state is not accepted while the
        trend of that peak (difference of its means over the two windows) reaches zero before the
        extrapolated failure and within the remaining blocks; the integration goes on instead.
        The life is extrapolated, N = n + (1 - D) / mean dD.
    After 60000 blocks the life is extrapolated in the same way, with a warning. Comparing only
    two consecutive blocks, as before November 2026, stopped during slow shakedowns and returned
    lives of 9e7 to 2.4e9 cycles instead of infinity.
  - These numerical settings are optional constants: npts_TVEIT (points per cycle, 1001),
    min_blocks_TVEIT (blocks before the steady-state test, 400; at least 2 x 200),
    max_blocks_TVEIT (60000), tol_alpha_TVEIT (drift tolerance, 0.01) and tol_damage_TVEIT (mean
    damage tolerance, 1e-3). A constants file may set them. A life extrapolated after
    max_blocks_TVEIT blocks raises a warning, which the pipeline reports in the notes of the run.
  - With 1001 points per cycle the life is about 0.02 decade shorter than the converged value
    (Table 1 NiCrMo steel, bending at sigma_a = 680 MPa).
  - G uses the deviatoric projection s s - (2/3) J2 I of dJ3/ds. Eq. (10) prints a different
    isotropic part, which has no effect because G is only contracted with deviatoric tensors.
  - The damage parameter passed to the universal curve is the fatigue index, the maximum over the
    cycle of (sigma_e + A I1) / f_m1 with alpha = 0 (beta + 1 of Eq. (1) for the virgin material).
    The model's own prediction is the native life.
  - A, C, K and L are read from a constants file as A_TVEIT, C_TVEIT, K_TVEIT and L_TVEIT.
    The calibration of C, K and L (Eq. (25), after Ottosen et al. 2008) is not implemented here.
  - The constants hold only together with the f_m1 and t_m1 they were calibrated with, so the
    constants file may also set Sigm1_TVEIT and Taum1_TVEIT. These replace the workbook Sigm1 and
    Taum1 for this model only. Without them, the workbook values are used.

The name under which the registry exposes the method is "TVEIT2025".
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..comparison import ComparisonCapabilities
from .stress_invariants import harmonic_parts

N_POINTS = 1001  # samples per closed cycle
WINDOW = 200  # blocks per window of the cyclic steady-state test
MIN_BLOCKS = 2 * WINDOW  # blocks before the steady-state test
MAX_BLOCKS = 60000
TOL_ALPHA = 0.01  # net back-stress drift over a window, relative to its movement in the window
TOL_DAMAGE = 1e-3  # relative change of the mean damage per block between two windows
_NORM_WEIGHTS = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)  # tensor norm of the six stored components
_SQRT3 = math.sqrt(3.0)

# Constants-file names of the numerical settings, the keyword of tveit2025_life and the default.
SETTINGS = {
    "npts_TVEIT": ("n_points", N_POINTS),
    "min_blocks_TVEIT": ("min_blocks", MIN_BLOCKS),
    "max_blocks_TVEIT": ("max_blocks", MAX_BLOCKS),
    "tol_alpha_TVEIT": ("tol_alpha", TOL_ALPHA),
    "tol_damage_TVEIT": ("tol_damage", TOL_DAMAGE),
}


@dataclass(frozen=True)
class TveitLife:
    """Native life in cycles (inf = no predicted failure), its termination rule and block count."""

    Nf: float
    termination: str
    blocks: int


def _cycle_components(S_mean, S_sin, S_cos, n_points) -> list[np.ndarray]:
    """Sampled cycle as six component series (xx, yy, zz, xy, yz, zx)."""
    phase = 2.0 * math.pi * np.linspace(0.0, 1.0, n_points)
    sin_t, cos_t = np.sin(phase), np.cos(phase)
    S_mean, S_sin, S_cos = (np.asarray(T, dtype=float) for T in (S_mean, S_sin, S_cos))
    index = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0))
    return [S_mean[i, j] + S_sin[i, j] * sin_t + S_cos[i, j] * cos_t for i, j in index]


def _effective_stress(hxx, hyy, hzz, hxy, hyz, hzx, lam):
    """Eq. (5) and its gradient G (Eq. 10, deviatoric projection) for the shifted deviator s - alpha.

    Returns (sigma_e, Gxx, Gyy, Gzz, Gxy, Gyz, Gzx); all zero when J2 < 1e-20.
    """
    tr3 = (hxx + hyy + hzz) / 3.0
    sxx = hxx - tr3
    syy = hyy - tr3
    szz = hzz - tr3
    sxy = hxy
    syz = hyz
    szx = hzx
    J2 = 0.5 * (sxx * sxx + syy * syy + szz * szz) + (sxy * sxy + syz * syz + szx * szx)
    if J2 < 1e-20:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    J3 = (
        sxx * (syy * szz - syz * syz)
        - sxy * (sxy * szz - syz * szx)
        + szx * (sxy * syz - syy * szx)
    )
    J2_3 = J2 * J2 * J2
    sigma_e = math.sqrt(3.0 * J2) - lam * (1.0 - 6.75 * J3 * J3 / J2_3)
    a = _SQRT3 / (2.0 * math.sqrt(J2)) - lam * 20.25 * J3 * J3 / (J2_3 * J2)
    b = lam * 13.5 * J3 / J2_3
    c = (2.0 / 3.0) * J2
    qxx = sxx * sxx + sxy * sxy + szx * szx
    qyy = sxy * sxy + syy * syy + syz * syz
    qzz = szx * szx + syz * syz + szz * szz
    qxy = sxx * sxy + sxy * syy + szx * syz
    qyz = sxy * szx + syy * syz + syz * szz
    qzx = sxx * szx + sxy * syz + szx * szz
    return (
        sigma_e,
        a * sxx + b * (qxx - c),
        a * syy + b * (qyy - c),
        a * szz + b * (qzz - c),
        a * sxy + b * qxy,
        a * syz + b * qyz,
        a * szx + b * qzx,
    )


def _block_life(
    comps, f_m1, A, C, K, L, lam, min_blocks, max_blocks, tol_alpha, tol_damage
) -> TveitLife:
    """Integration of Eqs. (6)-(12) over repeated blocks of the sampled cycle.

    The integration ends at failure (D >= 1 within a block), at shakedown (a block without any
    damage increment: the back-stress no longer moves, so every later block is the same) or in a
    cyclic steady state (the back-stress returns to where it was a window earlier, the mean
    damage per block no longer changes, and the peak of beta is not heading for zero within reach);
    otherwise after max_blocks blocks.
    """
    n = comps[0].size
    I1 = (comps[0] + comps[1] + comps[2]).tolist()
    tr3 = np.asarray(I1) / 3.0
    dev = [(comps[k] - tr3).tolist() for k in range(3)] + [comps[k].tolist() for k in range(3, 6)]
    s_xx, s_yy, s_zz, s_xy, s_yz, s_zx = dev
    w = _NORM_WEIGHTS

    alpha = [0.0] * 6
    D = 0.0
    block_damage: list[float] = []
    block_movement: list[float] = []  # path length of the back-stress within each block
    block_beta: list[float] = []  # peak of beta over the integrated points of each block
    alpha_end = [list(alpha)]  # back-stress at the end of each block (index 0: the start)
    block = 0
    for block in range(1, max_blocks + 1):
        D_start = D
        movement = 0.0
        beta_peak = -math.inf
        for i in range(n):
            hxx = s_xx[i] - alpha[0]
            hyy = s_yy[i] - alpha[1]
            hzz = s_zz[i] - alpha[2]
            hxy = s_xy[i] - alpha[3]
            hyz = s_yz[i] - alpha[4]
            hzx = s_zx[i] - alpha[5]
            sigma_e, gxx, gyy, gzz, gxy, gyz, gzx = _effective_stress(
                hxx, hyy, hzz, hxy, hyz, hzx, lam
            )
            beta = (sigma_e + A * I1[i] - f_m1) / f_m1
            D_previous = D
            if i > 0 and beta > beta_peak:
                beta_peak = beta
            if i > 0 and beta >= 0 and sigma_e > 1e-14:
                dI1 = I1[i] - I1[i - 1]
                dxx = s_xx[i] - s_xx[i - 1]
                dyy = s_yy[i] - s_yy[i - 1]
                dzz = s_zz[i] - s_zz[i - 1]
                dxy = s_xy[i] - s_xy[i - 1]
                dyz = s_yz[i] - s_yz[i - 1]
                dzx = s_zx[i] - s_zx[i - 1]
                G_s = gxx * hxx + gyy * hyy + gzz * hzz + 2.0 * (gxy * hxy + gyz * hyz + gzx * hzx)
                G_ds = gxx * dxx + gyy * dyy + gzz * dzz + 2.0 * (gxy * dxy + gyz * dyz + gzx * dzx)
                denominator = f_m1 + C * G_s
                if abs(denominator) > 1e-14:
                    d_beta = (G_ds + A * dI1) / denominator
                    if d_beta > 0:
                        step = d_beta * C
                        alpha[0] += step * hxx
                        alpha[1] += step * hyy
                        alpha[2] += step * hzz
                        alpha[3] += step * hxy
                        alpha[4] += step * hyz
                        alpha[5] += step * hzx
                        movement += step * math.sqrt(
                            w[0] * hxx * hxx
                            + w[1] * hyy * hyy
                            + w[2] * hzz * hzz
                            + w[3] * hxy * hxy
                            + w[4] * hyz * hyz
                            + w[5] * hzx * hzx
                        )
                        D = D + d_beta * K * math.exp(L * beta)
            if D_previous < 1 and D >= 1:
                # The increment from point i - 1 to point i crosses D = 1, at the fraction `part`
                # of it; i - 1 increments of the n - 1 in the cycle precede it.
                part = (1 - D_previous) / (D - D_previous) if D > D_previous else 1.0
                return TveitLife((block - 1) + (i - 1 + part) / (n - 1), "explicit_failure", block)
        damage = D - D_start
        block_damage.append(damage)
        block_movement.append(movement)
        block_beta.append(beta_peak)
        alpha_end.append(list(alpha))
        if block >= 2 and damage <= 0.0:
            return TveitLife(math.inf, "shakedown", block)
        if block >= max(min_blocks, 2 * WINDOW):
            drift = math.sqrt(
                sum(w[k] * (alpha_end[-1][k] - alpha_end[-1 - WINDOW][k]) ** 2 for k in range(6))
            )
            moved = sum(block_movement[-WINDOW:])
            last = sum(block_damage[-WINDOW:]) / WINDOW
            before = sum(block_damage[-2 * WINDOW : -WINDOW]) / WINDOW
            if drift <= tol_alpha * moved and abs(last - before) <= tol_damage * last:
                # A back-stress that swings in every block can still creep one way and end in a
                # shakedown. Keep integrating when the trend of the peak of beta reaches zero before
                # the extrapolated failure and within the remaining blocks.
                slope = (sum(block_beta[-WINDOW:]) - sum(block_beta[-2 * WINDOW : -WINDOW])) / (
                    WINDOW * WINDOW
                )
                remaining = (1 - D) / last
                if not (
                    slope < 0.0 and block_beta[-1] / -slope < min(remaining, max_blocks - block)
                ):
                    return TveitLife(block + remaining, "cyclic_steady_state", block)
    last = sum(block_damage[-WINDOW:]) / min(WINDOW, len(block_damage))
    if last <= 0.0:
        return TveitLife(math.inf, "max_block_extrapolation", block)
    return TveitLife(block + max(0.0, (1 - D) / last), "max_block_extrapolation", block)


def _check_constants(f_m1, t_m1, C, K):
    if not (f_m1 > 0.0 and t_m1 > 0.0):
        raise ValueError(f"TVEIT2025 needs f_m1 > 0 and t_m1 > 0, got f_m1 = {f_m1}, t_m1 = {t_m1}")
    if not (C >= 0.0 and K > 0.0):
        raise ValueError(f"TVEIT2025 needs C >= 0 and K > 0, got C = {C}, K = {K}")


def _whole_number(value, name: str, smallest: int) -> int:
    if not (float(value) == int(value) and int(value) >= smallest):
        raise ValueError(f"TVEIT2025 needs a whole number {name} >= {smallest}, got {value}")
    return int(value)


def tveit2025_life(
    S_mean,
    S_sin,
    S_cos,
    f_m1,
    t_m1,
    A,
    C,
    K,
    L,
    *,
    n_points=N_POINTS,
    min_blocks=MIN_BLOCKS,
    max_blocks=MAX_BLOCKS,
    tol_alpha=TOL_ALPHA,
    tol_damage=TOL_DAMAGE,
) -> TveitLife:
    """Native life of the single-harmonic cycle (Eqs. 1-12 and 24), in cycles."""
    f_m1, t_m1, A, C, K, L = (float(v) for v in (f_m1, t_m1, A, C, K, L))
    _check_constants(f_m1, t_m1, C, K)
    n_points = _whole_number(n_points, "n_points", 3)
    min_blocks = _whole_number(min_blocks, "min_blocks", 2)
    max_blocks = _whole_number(max_blocks, "max_blocks", 2)
    if not (tol_alpha > 0.0 and tol_damage > 0.0):
        raise ValueError(
            f"TVEIT2025 needs tol_alpha > 0 and tol_damage > 0, got {tol_alpha} and {tol_damage}"
        )
    lam = _SQRT3 * t_m1 - f_m1
    comps = _cycle_components(S_mean, S_sin, S_cos, n_points)
    return _block_life(
        comps, f_m1, A, C, K, L, lam, min_blocks, max_blocks, float(tol_alpha), float(tol_damage)
    )


def fatigue_index(S_mean, S_sin, S_cos, f_m1, t_m1, A, *, n_points=N_POINTS) -> float:
    """Maximum of (sigma_e + A I1) / f_m1 over the sampled cycle with alpha = 0 (Eqs. 1, 5, 24)."""
    f_m1, t_m1, A = float(f_m1), float(t_m1), float(A)
    lam = _SQRT3 * t_m1 - f_m1
    comps = _cycle_components(S_mean, S_sin, S_cos, _whole_number(n_points, "n_points", 3))
    I1 = comps[0] + comps[1] + comps[2]
    sigma_e = [
        _effective_stress(*(float(c[i]) for c in comps), lam)[0] for i in range(comps[0].size)
    ]
    return float(np.max((np.asarray(sigma_e) + A * I1) / f_m1))


class Tveit2025:
    """Tveit et al. (2025) continuum model: fatigue index for the universal curve, native life."""

    name = "TVEIT2025"
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_comparison=False, dp_reference="none"
    )
    optional_params = {
        "Sigm1_TVEIT": "Sigm1 of the workbook",
        "Taum1_TVEIT": "Taum1 of the workbook",
        **{key: str(default) for key, (_, default) in SETTINGS.items()},
    }

    def required_params(self) -> set[str]:
        return {"A_TVEIT", "C_TVEIT", "K_TVEIT", "L_TVEIT"}

    @staticmethod
    def _limit(params: Mapping[str, float], key: str) -> float:
        """The model's own fatigue limit from the constants file, else the workbook value."""
        value = params.get(f"{key}_TVEIT")
        if value is None:
            value = params.get(key)
        if value is None:
            raise ValueError(
                f"TVEIT2025 needs {key}_TVEIT in the constants file or {key} in the workbook"
            )
        return float(value)

    def evaluate_case(
        self, *, S0, E0, S1, E1, params: Mapping[str, float], R_list=None, harmonics=None
    ) -> CaseResult:
        S_mean, S_sin, S_cos = harmonic_parts(S0, S1, harmonics)
        f_m1, t_m1 = self._limit(params, "Sigm1"), self._limit(params, "Taum1")
        A = float(params["A_TVEIT"])
        settings = {
            keyword: default if params.get(key) is None else params[key]
            for key, (keyword, default) in SETTINGS.items()
        }
        prediction = tveit2025_life(
            S_mean,
            S_sin,
            S_cos,
            f_m1,
            t_m1,
            A,
            float(params["C_TVEIT"]),
            float(params["K_TVEIT"]),
            float(params["L_TVEIT"]),
            **settings,
        )
        if prediction.termination == "max_block_extrapolation":
            warnings.warn(
                "TVEIT2025: lives extrapolated after max_blocks_TVEIT blocks, before the damage "
                "per block became stationary. Raise max_blocks_TVEIT to check them",
                stacklevel=2,
            )
        index = fatigue_index(S_mean, S_sin, S_cos, f_m1, t_m1, A, n_points=settings["n_points"])
        return CaseResult(
            values={"": index},
            aux={"Nf": prediction.Nf, "blocks": float(prediction.blocks)},
        )

    def life_fit(self, aux, N, params: Mapping[str, float]) -> dict:
        """No fit: C, K and L come from the calibration record (A_TVEIT ... L_TVEIT)."""
        return {}

    def native_life(self, aux, params: Mapping[str, float], fit):
        return aux["Nf"]
