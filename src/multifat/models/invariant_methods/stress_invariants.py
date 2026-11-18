"""Amplitudes of stress invariants that the invariant-based criteria rely on (Papadopoulos et al. 1997).

The quantity sqrt(J2,a) gives the radius of the smallest hypersphere that encloses the deviatoric
load path when measured with the J2 norm ||s|| = sqrt(s:s / 2) (Eqs. 38-39). For a single-harmonic
history S(t) = M + A sin t + B cos t the path is a centred ellipse, hence the enclosing hypersphere
shares the centre dev(M) and its radius equals the semi-major axis, i.e. the largest eigenvalue of
the Gram matrix built from dev(A) and dev(B). In the bending-torsion case this reduces to Eq. (44).
"""

from __future__ import annotations

import numpy as np


def _deviatoric(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    return T - np.trace(T) / 3.0 * np.eye(3)


def sqrt_j2_amplitude(S_sin: np.ndarray, S_cos: np.ndarray) -> float:
    S_sin_dev, S_cos_dev = _deviatoric(S_sin), _deviatoric(S_cos)
    g_ss, g_cc, g_sc = (
        0.5 * np.sum(x * y)
        for x, y in ((S_sin_dev, S_sin_dev), (S_cos_dev, S_cos_dev), (S_sin_dev, S_cos_dev))
    )
    return float(np.sqrt(max(0.0, 0.5 * (g_ss + g_cc + np.hypot(g_ss - g_cc, 2.0 * g_sc)))))


def hydrostatic_amplitude_mean(S_mean, S_sin, S_cos) -> tuple[float, float]:
    """Return sigma_H,a and sigma_H,m for the single-harmonic history (Eq. 45: max = a + m)."""
    trace_mean = np.trace(np.asarray(S_mean, float))
    sigma_H_a = (
        np.hypot(
            np.trace(np.asarray(S_sin, float)),
            np.trace(np.asarray(S_cos, float)),
        )
        / 3.0
    )
    return float(sigma_H_a), float(trace_mean / 3.0)


PEAK_SAMPLES = 720  # time samples bracketing the peak before the golden-section refinement


def peak_principal_stress(S_mean, S_sin, S_cos) -> float:
    """max over eta of the largest principal stress of S_mean + S_sin sin(eta) + S_cos cos(eta)."""

    def largest_eig(phase):
        phase = np.atleast_1d(phase)
        tensor = (
            S_mean + np.sin(phase)[:, None, None] * S_sin + np.cos(phase)[:, None, None] * S_cos
        )
        return np.linalg.eigvalsh(tensor)[:, -1]

    samples = np.linspace(0.0, 2.0 * np.pi, PEAK_SAMPLES, endpoint=False)
    largest = largest_eig(samples)
    # coincident maxima give the same value, which is all this needs
    top = int(np.argmax(largest))
    spacing = samples[1] - samples[0]
    lo, hi = samples[top] - spacing, samples[top] + spacing
    for _ in range(60):  # golden-section search inside the bracket
        left = hi - 0.618034 * (hi - lo)
        right = lo + 0.618034 * (hi - lo)
        if largest_eig(left)[0] >= largest_eig(right)[0]:
            hi = right
        else:
            lo = left
    return max(float(largest[top]), float(largest_eig(0.5 * (lo + hi))[0]))


def harmonic_parts(S0, S1, harmonics):
    if harmonics is None:
        start, end = np.asarray(S0, dtype=float), np.asarray(S1, dtype=float)
        return 0.5 * (start + end), 0.5 * (start - end), np.zeros((3, 3))
    return tuple(np.asarray(harmonics[key], dtype=float) for key in ("S_mean", "S_sin", "S_cos"))
