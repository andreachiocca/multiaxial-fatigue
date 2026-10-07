"""Analytical limiting cases for the BP sphere integral."""

import unittest

import numpy as np

from fatigue.models.integral_methods.bohme_papuga import BohmePapuga, _strengths


class BohmePapugaTests(unittest.TestCase):
    def setUp(self):
        self.model = BohmePapuga()
        # Values are illustrative; the four strengths represent maximum stress
        # for R=0 and amplitude for R=-1, respectively.
        self.params = {"Sigm1": 240.0, "Taum1": 145.0, "Sig0": 400.0, "Tau0": 265.0}

    def damage(self, mean, sin, cos=None, params=None):
        if cos is None:
            cos = np.zeros((3, 3))
        return self.model.evaluate_case(
            S0=mean + sin, S1=mean - sin, E0=None, E1=None,
            params=params or self.params,
            harmonics={"S_mean": mean, "S_sin": sin, "S_cos": cos},
        ).values[""]

    def test_pure_loading_limits(self):
        zero = np.zeros((3, 3))
        axial = np.diag([1.0, 0.0, 0.0])
        torsion = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        s, t, s0, t0 = (_strengths(self.params))
        for mean, amplitude in (
            (zero, s * axial), (zero, t * torsion),
            (s0 / 2 * axial, s0 / 2 * axial),
            (t0 / 2 * torsion, t0 / 2 * torsion),
        ):
            with self.subTest(mean=mean.tolist(), amplitude=amplitude.tolist()):
                self.assertAlmostEqual(self.damage(mean, amplitude), 1.0, places=9)

    def test_phase_origin_does_not_change_damage(self):
        mean = np.diag([20.0, -10.0, 0.0])
        sin = np.array([[80.0, 40.0, 0.0], [40.0, -20.0, 0.0], [0.0, 0.0, 0.0]])
        cos = np.array([[15.0, -55.0, 0.0], [-55.0, 5.0, 0.0], [0.0, 0.0, 0.0]])
        # Shifting the time origin by 90 degrees rotates sin/cos coefficients.
        self.assertAlmostEqual(self.damage(mean, sin, cos), self.damage(mean, cos, -sin), places=12)

    def test_estimated_zero_to_maximum_strengths(self):
        params = {"Sigm1": 240.0, "Taum1": 145.0, "Su": 580.0}
        s, t, s0, t0 = _strengths(params)
        self.assertAlmostEqual(s0, 4 * s * 580 / (s + 2 * 580))
        self.assertAlmostEqual(4 * t / t0 - 2 * s / s0, 1.0)


if __name__ == "__main__":
    unittest.main()
