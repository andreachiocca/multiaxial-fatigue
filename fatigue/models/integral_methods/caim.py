"""Critical Amplitude–Integral Mean, Böhme et al. (2026).

Coefficients and clipping follow the supplied damage_criteria_CAIM.m.
Sphere quadrature replaces its uniform angular sum; a numerical plane search
refines the amplitude maximum. Single-harmonic stress histories only.
"""
import numpy as np

from .bohme_papuga import BohmePapuga, _plane_amplitudes


class CAIM(BohmePapuga):
    """Hybrid criterion, normalized by the reversed axial fatigue strength."""

    name = "CAIM"

    def _amplitude_coefficients(self, k2):
        return k2, 1.0 if k2 < 2.0 else k2 - k2**2 / 4.0

    def _combine(self, a, b, c, d, sigma_a, sigma_m, tau_a, tau_m,
                 weights, normals, mean, sin, cos):
        amplitude = a * tau_a**2 + b * sigma_a**2
        # Refine several grid maxima to reduce orientation-dependent grid error.
        seeds = normals[np.argsort(amplitude)[-12:]]
        seeds = np.concatenate((seeds, np.eye(3)))
        offsets = np.concatenate((np.eye(3), -np.eye(3), np.zeros((1, 3))))
        step = 0.15
        for _ in range(24):
            candidates = seeds[:, None, :] + step * offsets[None, :, :]
            candidates /= np.linalg.norm(candidates, axis=2, keepdims=True)
            sa, _, ta, _ = _plane_amplitudes(mean, sin, cos, candidates.reshape(-1, 3))
            values = (a * ta**2 + b * sa**2).reshape(len(seeds), -1)
            seeds = candidates[np.arange(len(seeds)), np.argmax(values, axis=1)]
            step *= 0.65
        critical = max(float(np.max(amplitude)), float(np.max(values)))
        mean_integral = np.dot(weights, c * tau_a * tau_m + d * sigma_a * sigma_m)
        return critical + mean_integral
