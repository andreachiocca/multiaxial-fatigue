"""Analytical limits and independent angular implementation of MATLAB equations."""
import unittest
import numpy as np
from fatigue.models.integral_methods.caim import CAIM
from fatigue.models.integral_methods.bohme_papuga import BohmePapuga
from fatigue.models.registry import get_models


class CAIMTests(unittest.TestCase):
    params = dict(Sigm1=240., Taum1=145., Sig0=400., Tau0=265.)

    def damage(self, mean, sin, cos=None, params=None, model=None):
        return (model or CAIM()).evaluate_case(
            S0=mean+sin, S1=mean-sin, E0=None, E1=None,
            params=params or self.params,
            harmonics=dict(S_mean=mean, S_sin=sin,
                           S_cos=np.zeros((3, 3)) if cos is None else cos)).values['']

    def test_four_loading_limits_and_coefficient_branches(self):
        axial = np.diag([1., 0., 0.])
        torsion = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]])
        for t in (145., 180.):
            params = dict(self.params, Taum1=t)
            for mean, amp in ((0*axial, 240*axial), (0*axial, t*torsion),
                              (200*axial, 200*axial), (132.5*torsion, 132.5*torsion)):
                with self.subTest(t=t, amp=amp.tolist()):
                    self.assertAlmostEqual(self.damage(mean, amp, params=params), 1., places=6)

    def test_dense_matlab_equations(self):
        mean = np.diag([35., -12., 5.])
        sin = np.array([[110., 45., 0.], [45., -20., 8.], [0., 8., 15.]])
        cos = np.array([[15., -65., 5.], [-65., 10., 0.], [5., 0., -5.]])
        # Independent midpoint angular quadrature, with SVD for shear ellipse.
        step = np.pi/180
        gamma, phi = np.meshgrid((np.arange(180)+.5)*step, np.arange(360)*step)
        n = np.stack((np.sin(gamma)*np.cos(phi), np.sin(gamma)*np.sin(phi),
                      np.cos(gamma)), axis=-1).reshape(-1, 3)
        normal, shear = [], []
        for tensor in (mean, sin, cos):
            traction = np.einsum('ij,nj->ni', tensor, n)
            normal.append(np.einsum('ni,ni->n', traction, n))
            shear.append(traction-normal[-1][:, None]*n)
        sm = normal[0]
        sa = np.hypot(normal[1], normal[2])
        tm = np.linalg.norm(shear[0], axis=1)
        ta = np.linalg.svd(np.stack(shear[1:], axis=-1), compute_uv=False)[:, 0]
        a = (240/145)**2
        b = a-a*a/4
        c = 2.5*a*((290/265)**2-1)
        d = 5*((480/400)**2-1-2*c/15)
        fc = np.max(a*ta**2+b*sa**2)
        fi = np.sum((c*ta*tm+d*sa*sm)*np.sin(gamma).ravel())*step**2/(4*np.pi)
        expected = np.sqrt(max(0., fc+fi))/240
        self.assertAlmostEqual(self.damage(mean, sin, cos), expected, delta=2e-4)
        self.assertAlmostEqual(self.damage(mean, sin, cos), self.damage(mean, cos, -sin), places=12)
        q, _ = np.linalg.qr(np.array([[1., 2., 3.], [2., -1., 1.], [1., 1., -1.]]))
        rotated = [q @ x @ q.T for x in (mean, sin, cos)]
        self.assertAlmostEqual(self.damage(*rotated), self.damage(mean, sin, cos), delta=2e-4)

    def test_distinct_registered_models_and_clipping(self):
        self.assertEqual([m.name for m in get_models(['BP', 'CAIM'])], ['BP', 'CAIM'])
        amp = 50*np.eye(3)
        self.assertNotAlmostEqual(self.damage(0*amp, amp), self.damage(0*amp, amp, model=BohmePapuga()))
        self.assertEqual(self.damage(-10000*np.eye(3), amp), 0.)
        self.assertEqual(self.damage(amp, 0*amp), 0.)


if __name__ == '__main__':
    unittest.main()
