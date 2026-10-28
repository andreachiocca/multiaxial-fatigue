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
  - Blocks are repeated until D >= 1 inside a block (the life is then interpolated linearly in the
    last increment) or until the back-stress and the damage per block are stationary (relative
    changes below 1e-6 and 1e-5 after at least 8 blocks). The life is then extrapolated with the
    damage of the last block, N = n + (1 - D) / dD, and the same extrapolation is applied after
    10000 blocks. A non-positive dD gives an infinite life.
  - These numerical settings are optional constants with the defaults above: npts_TVEIT (points
    per cycle), min_blocks_TVEIT, max_blocks_TVEIT, tol_alpha_TVEIT and tol_damage_TVEIT. A
    constants file may set them. A life extrapolated after max_blocks_TVEIT blocks raises a
    warning, which the pipeline reports in the notes of the run.
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
from .stress_invariants import harmonic_parts

N_POINTS = 1001  # samples per closed cycle
MIN_BLOCKS = 8
MAX_BLOCKS = 10000
TOL_ALPHA = 1e-6  # relative change of the back-stress between blocks
TOL_DAMAGE = 1e-5  # relative change of the damage per block
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
    """Integration of Eqs. (6)-(12) over repeated blocks of the sampled cycle."""
    n = comps[0].size
    I1 = (comps[0] + comps[1] + comps[2]).tolist()
    tr3 = np.asarray(I1) / 3.0
    dev = [(comps[k] - tr3).tolist() for k in range(3)] + [comps[k].tolist() for k in range(3, 6)]
    s_xx, s_yy, s_zz, s_xy, s_yz, s_zx = dev

    alpha = [0.0] * 6
    D = 0.0
    previous_alpha = None
    block_damage: list[float] = []
    alpha_change: list[float] = []
    stable_block = math.nan
    fail_block = math.nan
    fail_fraction = math.nan
    block = 0
    for block in range(1, max_blocks + 1):
        D_start = D
        failed = False
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
                        D = D + d_beta * K * math.exp(L * beta)
            if D_previous < 1 and D >= 1:
                failed = True
                part = (1 - D_previous) / (D - D_previous) if D > D_previous else 1.0
                fail_fraction = 0.0 if i == 0 else (i + part) / (n - 1)
                break
        block_damage.append(D - D_start)
        if previous_alpha is not None:
            weights = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
            change = math.sqrt(
                sum((alpha[k] - previous_alpha[k]) ** 2 * weights[k] for k in range(6))
            )
            size = max(math.sqrt(sum(alpha[k] ** 2 * weights[k] for k in range(6))), 1e-12)
            alpha_change.append(change / size)
        else:
            alpha_change.append(math.nan)
        previous_alpha = list(alpha)
        if failed:
            fail_block = block
            break
        if block >= max(2, min_blocks):
            damage_change = abs(block_damage[-1] - block_damage[-2]) / max(
                abs(block_damage[-1]), 1e-14
            )
            if (
                not math.isnan(alpha_change[-1])
                and alpha_change[-1] < tol_alpha
                and damage_change < tol_damage
            ):
                stable_block = block
                break
    if not math.isnan(fail_block):
        return TveitLife((fail_block - 1) + fail_fraction, "explicit_failure", block)
    if math.isnan(stable_block):
        stable_block = block
        termination = "max_block_extrapolation"
    else:
        termination = "stabilized_extrapolation"
    damage_per_block = block_damage[int(stable_block) - 1]
    if damage_per_block <= 0:
        return TveitLife(math.inf, termination, block)
    return TveitLife(stable_block + max(0.0, (1 - D) / damage_per_block), termination, block)


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
    max_blocks = _whole_number(max_blocks, "max_blocks", min_blocks)
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
