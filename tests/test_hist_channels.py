"""Scalar channels of a tensor history and the proportionality verdict."""

from __future__ import annotations

import numpy as np
import pytest

from multifat.history.channels import (
    component,
    dominant,
    normal_stress,
    resolve_channel,
    signed_von_mises,
)
from multifat.history.proportionality import Proportionality, proportionality
from multifat.history.signal import TimeHistory

TIME = np.linspace(0.0, 4.0, 401)
LAM = np.sin(2.0 * np.pi * TIME) + 0.3 * np.sin(6.0 * np.pi * TIME + 1.0)


def _sym(rng):
    a = rng.normal(size=(3, 3)) * 40.0
    return 0.5 * (a + a.T)


def _proportional(seed=0, static=True):
    rng = np.random.default_rng(seed)
    A = _sym(rng)
    static_part = _sym(rng) if static else np.zeros((3, 3))
    return LAM[:, None, None] * A + static_part, A


def _weights_dot(A, B):
    return float(np.sum(A * B))


def test_dominant_channel_is_affine_in_the_load_factor():
    S, A = _proportional()
    channel = dominant(S)
    slope, intercept = np.polyfit(LAM, channel, 1)
    assert abs(slope) > 1e-3
    assert np.max(np.abs(channel - (slope * LAM + intercept))) < 1e-9 * np.max(np.abs(channel))


def test_dominant_channel_sign_is_fixed_by_the_largest_component():
    S, A = _proportional(seed=3, static=False)
    c = dominant(S)
    assert np.allclose(dominant(-S), -c)
    weights = np.array([1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)])
    six = np.array([A[0, 0], A[1, 1], A[2, 2], A[0, 1], A[1, 2], A[0, 2]]) * weights
    idx = int(np.argmax(np.abs(six)))
    assert np.corrcoef(c, LAM * six[idx])[0, 1] > 0.0


def test_von_mises_is_signed_for_uniaxial_tension_and_compression():
    S = np.zeros((3, 3, 3))
    S[:, 2, 2] = [100.0, -100.0, 0.0]
    assert np.allclose(signed_von_mises(S), [100.0, -100.0, 0.0])


def test_von_mises_zero_trace_uses_principal_sign():
    S = np.zeros((2, 3, 3))
    S[0, 0, 2] = S[0, 2, 0] = 50.0  # pure shear: trace zero, max principal positive
    S[1, 0, 2] = S[1, 2, 0] = -50.0
    assert np.allclose(signed_von_mises(S), np.sqrt(3.0) * 50.0)


def test_component_and_normal_channels_match_direct_projection():
    S, _ = _proportional(seed=5)
    assert np.array_equal(component(S, 1, 2), S[:, 1, 2])
    n = np.array([1.0, 2.0, -2.0])
    unit = n / np.linalg.norm(n)
    expected = np.einsum("i,tij,j->t", unit, S, unit)
    assert np.allclose(normal_stress(S, n), expected)
    with pytest.raises(ValueError, match="zero"):
        normal_stress(S, np.zeros(3))


def test_resolve_channel_forms_and_unknown_rejected():
    S, _ = _proportional(seed=7)
    h = TimeHistory.from_arrays(TIME, S, stress_unit="MPa")
    assert np.allclose(resolve_channel("dominant", h), dominant(S))
    assert np.allclose(resolve_channel("von_mises", h), signed_von_mises(S))
    assert np.allclose(resolve_channel(("component", 0, 1), h), S[:, 0, 1])
    assert np.allclose(resolve_channel(("normal", (0, 0, 1)), h), S[:, 2, 2])
    assert np.allclose(resolve_channel(lambda s: s[:, 0, 0] * 2.0, h), 2.0 * S[:, 0, 0])
    with pytest.raises(ValueError, match="channel"):
        resolve_channel("nonsense", h)
    with pytest.raises(ValueError, match="shape"):
        resolve_channel(lambda s: s[:, 0], h)


def test_proportional_with_static_hydrostatic_offset_is_proportional():
    S, _ = _proportional(seed=1)
    S = S + 123.0 * np.eye(3)
    verdict = proportionality(S)
    assert (
        isinstance(verdict, Proportionality) and verdict.proportional and verdict.rank_ratio < 1e-12
    )


def test_noise_1e9_does_not_flip_proportionality():
    S, _ = _proportional(seed=2)
    rng = np.random.default_rng(0)
    noise = rng.normal(size=S.shape)
    noise = 0.5 * (noise + np.swapaxes(noise, -1, -2))
    verdict = proportionality(S + 1e-9 * np.max(np.abs(S)) * noise)
    assert verdict.proportional and 0.0 < verdict.rank_ratio < 1e-6


def test_elliptical_shear_path_is_not_proportional():
    S = np.zeros((len(TIME), 3, 3))
    S[:, 0, 2] = S[:, 2, 0] = 80.0 * np.sin(2.0 * np.pi * TIME)
    S[:, 1, 2] = S[:, 2, 1] = 40.0 * np.cos(2.0 * np.pi * TIME)
    verdict = proportionality(S)
    assert not verdict.proportional and verdict.rank_ratio == pytest.approx(0.5, rel=1e-2)


def test_constant_history_is_proportional_ratio_zero():
    S = np.tile(np.diag([10.0, 20.0, 30.0]), (50, 1, 1))
    verdict = proportionality(S)
    assert verdict.proportional and verdict.rank_ratio == 0.0
    assert np.all(dominant(S) == 0.0)


def test_tolerance_is_a_parameter():
    S = np.zeros((len(TIME), 3, 3))
    S[:, 2, 2] = 100.0 * np.sin(2.0 * np.pi * TIME)
    S[:, 0, 2] = S[:, 2, 0] = 1e-3 * np.cos(2.0 * np.pi * TIME)
    assert not proportionality(S, tol=1e-6).proportional
    assert proportionality(S, tol=1e-4).proportional
