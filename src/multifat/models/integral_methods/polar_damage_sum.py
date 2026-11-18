"""Albinmousa & Al Hussain (2021) polar damage sum of the SWT and Fatemi-Socie parameters (N59).

J. Albinmousa, M. Al Hussain, Forces in Mechanics 4 (2021) 100025,
doi:10.1016/j.finmec.2021.100025, Eqs. (7)-(14), (17)-(19)::

    D_SWT(phi) = sigma_max(phi) Delta eps_max(phi) / 2                                    (8)
    D_FS(phi)  = Delta gamma_max(phi)/2 (1 + k_phi sigma_n,max(phi) / sigma_yt)           (9)
    D^s = 1/2 int_0^2pi D(phi)^2 dphi                                                       (12)-(14)
    D^s = alpha N_f^beta                                                                     (17)-(19)

phi is the orientation of a plane perpendicular to the free surface (the x-y plane), the stress and
strain on it follow from the in-plane transformation of Eqs. (7a)-(7d), and D^s is the area enclosed
by the polar curve D(phi). Every plane of the surface contributes, instead of the critical plane
only. k_phi is fitted so that the D^s - life curves of the axial and shear tests coincide (Fig. 4),
and alpha, beta by a power fit (Eqs. 18-19).

For one harmonic, on the plane with unit normal n = (cos phi, sin phi, 0) and in-plane tangent t:
sigma(phi) = n.S.n, eps(phi) = n.E.n, gamma(phi) = 2 t.E.n, each m + a sin + b cos; the maximum over
the cycle is m + sqrt(a^2 + b^2) and half the range is sqrt(a^2 + b^2).

Repository choices:
  - elastic version (maintainer decision): the source is a strain-controlled low-cycle method built on measured stabilised hysteresis
    loops (elastic and plastic Poisson ratios, Eq. 7e). Here the strains are the linear-elastic
    strains of the imported harmonic stresses, so the model is a high-cycle, elastic reading of the
    concept, not the published LCF procedure;
  - k_phi (`k_PDS`) is calibrated like k_FS, by minimising the life scatter of the R = -1 tension
    and torsion failures (cli/calibrate.py), which is the paper's collapse of the axial and shear
    curves; sigma_yt is the workbook yield strength `Sy`. On a fatigue-limit dataset k_phi is set,
    like k_FS, so that the R = -1 tension and torsion limits give the same D^s;
  - the life law D^s = alpha N^beta is the shared master curve, fitted on the R = -1 uniaxial
    tension failures (the paper fits axial and shear tests together);
  - D(phi) is taken as max(D(phi), 0) before squaring: a negative D (a compressive sigma_max for
    SWT, or 1 + k_phi sigma_n,max/sigma_yt < 0 for FS) would otherwise add damage;
  - plane stress on the x-y surface: a stress with a z component raises LoadingScopeError. An
    in-phase equibiaxial state has no in-plane shear on any surface plane, so D_FS = 0 there (the
    two-dimensional polar sum does not see the out-of-plane shear) and its load factor is undefined;
  - the integral uses N_PHI equally spaced planes over [0, 2 pi) (exact for the trigonometric
    polynomials of one harmonic up to high order).

User guides: docs/method-guides/PDS_SWT.md and docs/method-guides/PDS_FS.md.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from ..base import CaseResult
from ..errors import LoadingScopeError
from ..invariant_methods.stress_invariants import harmonic_parts

N_PHI = 720
SCOPE_TOL = 1e-9
ZERO_RTOL = 1e-9  # in-plane shear strain below this fraction of the strain scale is zero
_PHI = np.linspace(0.0, 2.0 * np.pi, N_PHI, endpoint=False)
_N = np.column_stack([np.cos(_PHI), np.sin(_PHI), np.zeros(N_PHI)])
_T = np.column_stack([-np.sin(_PHI), np.cos(_PHI), np.zeros(N_PHI)])


def _normal(T) -> np.ndarray:
    return np.einsum("ni,ij,nj->n", _N, np.asarray(T, float), _N)


def _shear(T) -> np.ndarray:
    return np.einsum("ni,ij,nj->n", _T, np.asarray(T, float), _N)


def polar_damage(
    kind: str, S_mean, S_sin, S_cos, E_sin, E_cos, k_PDS=None, sigma_yt=None
) -> np.ndarray:
    """D(phi) of Eq. (8) (``kind="SWT"``) or Eq. (9) (``kind="FS"``) on the N_PHI surface planes."""
    sigma_max = _normal(S_mean) + np.hypot(_normal(S_sin), _normal(S_cos))
    if kind == "SWT":
        eps_amp = np.hypot(_normal(E_sin), _normal(E_cos))  # Delta eps_max(phi)/2
        return sigma_max * eps_amp
    if kind == "FS":
        gam_amp = 2.0 * np.hypot(_shear(E_sin), _shear(E_cos))  # Delta gamma_max(phi)/2
        # round-off of a state without in-plane shear (in-phase equibiaxial) is no damage
        scale = max(float(np.max(np.abs(E_sin))), float(np.max(np.abs(E_cos))))
        if np.isfinite(scale):  # undefined strains (no E, nu) stay NaN
            gam_amp = np.where(gam_amp > ZERO_RTOL * scale, gam_amp, 0.0)
        return gam_amp * (1.0 + float(k_PDS) * sigma_max / float(sigma_yt))
    raise ValueError(f"polar damage kind must be 'SWT' or 'FS', got {kind!r}")


def polar_damage_sum(
    kind: str, S_mean, S_sin, S_cos, E_sin, E_cos, k_PDS=None, sigma_yt=None
) -> float:
    """D^s = 1/2 int_0^2pi max(D(phi), 0)^2 dphi, Eqs. (12)-(14)."""
    D = np.maximum(polar_damage(kind, S_mean, S_sin, S_cos, E_sin, E_cos, k_PDS, sigma_yt), 0.0)
    return float(0.5 * np.mean(D**2) * 2.0 * np.pi)


def _check_plane_stress(S_mean, S_sin, S_cos) -> None:
    scale = max(float(np.max(np.abs(T))) for T in (S_mean, S_sin, S_cos)) or 1.0
    if max(float(np.max(np.abs(T[2, :]))) for T in (S_mean, S_sin, S_cos)) > SCOPE_TOL * scale:
        raise LoadingScopeError("PDS: plane stress on the x-y surface only, as in the source")


class _PolarDamageSum:
    name = ""
    kind = ""

    def _damage(self, params, S_mean, S_sin, S_cos, E_sin, E_cos) -> float:
        raise NotImplementedError

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
        _check_plane_stress(S_mean, S_sin, S_cos)
        if harmonics is not None and harmonics.get("E_sin") is not None:
            E_sin, E_cos = (np.asarray(harmonics[k], float) for k in ("E_sin", "E_cos"))
        else:
            E_sin, E_cos = 0.5 * (np.asarray(E0, float) - np.asarray(E1, float)), np.zeros((3, 3))
        return CaseResult(values={"": self._damage(params, S_mean, S_sin, S_cos, E_sin, E_cos)})


class PolarDamageSumSWT(_PolarDamageSum):
    name = "PDS_SWT"
    kind = "SWT"
    fitted_constants = ()
    input_classes = ("elastic constants",)

    def required_params(self) -> set[str]:
        return set()

    def _damage(self, params, S_mean, S_sin, S_cos, E_sin, E_cos) -> float:
        return polar_damage_sum("SWT", S_mean, S_sin, S_cos, E_sin, E_cos)


class PolarDamageSumFS(_PolarDamageSum):
    name = "PDS_FS"
    kind = "FS"
    fitted_constants = ("k_PDS",)
    input_classes = ("elastic constants", "static strength")

    def required_params(self) -> set[str]:
        return {"Sy", "k_PDS"}

    def _damage(self, params, S_mean, S_sin, S_cos, E_sin, E_cos) -> float:
        k_PDS, Sy = float(params["k_PDS"]), float(params["Sy"])
        if not Sy > 0.0:
            raise ValueError(f"PDS_FS needs Sy > 0, got Sy = {Sy}")
        if k_PDS < 0.0:
            raise ValueError(f"PDS_FS needs k_PDS >= 0, got k_PDS = {k_PDS}")
        return polar_damage_sum("FS", S_mean, S_sin, S_cos, E_sin, E_cos, k_PDS, Sy)
