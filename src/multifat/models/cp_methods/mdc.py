# multifat/models/cp_methods/mdc.py
"""Critical-plane criterion of Margetin, Ďurka & Chmelko (2016), Frattura ed
Integrità Strutturale 37:146-152, doi:10.3221/IGF-ESIS.37.20 (roadmap G06),
Eqs. (4)-(8):

    tau_r = sqrt(k sigma_a,cr^2 + tau_a,cr^2) = tau'_f (2 N_f)^b_tau,
    k = (2 tau'_f / sigma'_f)^2 (2 sigma_a,cr / sigma'_f)^(2 (b_tau - b_sigma) / b_sigma) - 1

The critical plane is the one carrying the largest shear-stress amplitude, with the
axial and torsional Basquin curves sigma_a = sigma'_f (2N)^b_sigma and
tau_a = tau'_f (2N)^b_tau (Eqs. 5, 6). The factor k is stress-level dependent: it is
built so that pure axial loading (sigma_a,cr = tau_a,cr = sigma_a/2) reproduces the
axial curve, while pure torsion reproduces the torsional one. Eqs. (7)-(8) use
sigma_max in the place of the sigma_a,cr appearing in Eq. (4); taking
sigma_max = sigma_a,cr reproduces the Nf_com values of Table 4 (the paper tests fully
reversed loading only), so the criterion carries no mean-stress effect. If the
radicand becomes negative (k < 0 at low normal stress), it is clamped to zero. The
damage parameter is tau_r (MPa); equal plane metrics are broken in favour of the
larger tau_r.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from multifat.models.errors import MissingParamsError
from multifat.models.utils.batch import check_planes, pypow2

from ..base import CPPlaneResult
from ..comparison import ComparisonCapabilities

_PARAM_KEYS = ("sigma_f", "b_fat", "tau_f", "b0")


def _extract_constants(
    params: Mapping[str, float],
) -> tuple[float, float, float, float]:
    collected = []
    for key in _PARAM_KEYS:
        raw = params.get(key, None)
        if raw is None or not np.isfinite(float(raw)):
            raise MissingParamsError("MDC", [key], "Axial and torsional Basquin constants")
        collected.append(float(raw))
    sigma_f, b_sigma, tau_f, b_tau = collected
    if sigma_f <= 0.0 or tau_f <= 0.0 or b_sigma >= 0.0 or b_tau >= 0.0:
        raise ValueError("MDC requires sigma_f, tau_f > 0 and b_fat, b0 < 0")
    return sigma_f, b_sigma, tau_f, b_tau


def mdc_equivalent_shear(
    sigma_a, tau_a, sigma_f: float, b_sigma: float, tau_f: float, b_tau: float
):
    """Return tau_r of Eq. (8) from the normal and shear amplitudes on the critical plane."""
    sigma_a = np.abs(np.asarray(sigma_a, dtype=float))
    tau_a = np.asarray(tau_a, dtype=float)
    b_ratio = 2.0 * (b_tau - b_sigma) / b_sigma
    guarded = np.where(sigma_a > 0.0, sigma_a, 1.0)
    k_sigma_a_sq = np.where(
        sigma_a > 0.0,
        ((2.0 * tau_f / sigma_f) ** 2 * (2.0 * guarded / sigma_f) ** b_ratio - 1.0) * sigma_a**2,
        0.0,
    )
    return np.sqrt(np.maximum(k_sigma_a_sq + tau_a**2, 0.0))


def mdc_life(tau_r, tau_f: float, b_tau: float):
    """Return N_f from tau_r = tau'_f (2 N_f)^b_tau, Eq. (8)."""
    return 0.5 * (np.asarray(tau_r, dtype=float) / tau_f) ** (1.0 / b_tau)


class MDC:
    """MDC critical-plane criterion: maximum shear-amplitude plane, damage tau_r (Eq. 8)."""

    name = "MDC"
    comparison_capabilities = ComparisonCapabilities(
        life_prediction="native", dp_reference="model_law"
    )

    def required_params(self) -> set[str]:
        return set(_PARAM_KEYS)

    def native_life(self, aux, params: Mapping[str, float], fit=None):
        """Return N_f from tau_r through the torsional Basquin curve, Eq. (8)."""
        _, _, tau_f, b_tau = _extract_constants(params)
        return mdc_life(aux["dp"], tau_f, b_tau)

    def expected_dp_at_life(self, N_exp, params: Mapping[str, float]):
        """Eq. (8) equivalent shear strength at experimental life (MPa)."""
        _, _, tau_f, b_tau = _extract_constants(params)
        n = np.asarray(N_exp, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return np.where(n > 0, tau_f * np.power(2.0 * n, b_tau), np.nan)

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
        sigma_f, b_sigma, tau_f, b_tau = _extract_constants(params)
        d_tau_xz = S0r[:, 0, 2] - S1r[:, 0, 2]
        d_tau_yz = S0r[:, 1, 2] - S1r[:, 1, 2]
        d_tau = np.sqrt(pypow2(d_tau_xz) + pypow2(d_tau_yz))
        sigma_a = 0.5 * np.abs(S0r[:, 2, 2] - S1r[:, 2, 2])
        return (
            mdc_equivalent_shear(sigma_a, 0.5 * d_tau, sigma_f, b_sigma, tau_f, b_tau),
            d_tau,
        )
