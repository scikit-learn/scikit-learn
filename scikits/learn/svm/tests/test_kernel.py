from unittest import TestCase
from numpy.testing import assert_array_equal
import numpy as N

from ..kernel import Linear, Polynomial,Sigmoid, RBF, Custom

def kernelfunc(x, y):
    return 8 * N.dot(x, y.T)

class TestKernel(TestCase):
    def test_linear_kernel(self):
        kernel = Linear()
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 4.)

    def test_polynomial_kernel(self):
        kernel = Polynomial(degree=6, gamma=1.0, coef0=1.0)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 15625.)

    def test_sigmoid_kernel(self):
        kernel = Sigmoid(gamma=0.2, coef0=0.3)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 0.80049902)

    def test_rbf_kernel(self):
        kernel = RBF(gamma=1.0)
        x, y = N.array([2.]), N.array([3.])
        self.assertAlmostEqual(kernel(x, y), N.exp(-1.))

    def test_custom_kernel(self):
        kernel = Custom(kernelfunc)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 32.0)

    def test_multidim_input(self):
        kernels = [
            Linear(),
            Polynomial(degree=6, gamma=1.0, coef0=1.0),
            Sigmoid(gamma=0.2, coef0=0.3),
            RBF(gamma=1.0),
            Custom(kernelfunc)
            ]
        args = [
            N.random.randn(10),
            N.random.randn(1, 10),
            N.random.randn(5, 10)
            ]
        for kernel in kernels:
            self.assert_(type(repr(kernel)) is str)
            for i, x in enumerate(args):
                zshape0 = N.atleast_2d(x).shape[0]
                for y in args[i:]:
                    zshape1 = N.atleast_2d(y).shape[0]
                    z = kernel(x, y)
                    self.assertEqual(z.shape[0], zshape0)
                    self.assertEqual(z.shape[1], zshape1)
                    u = kernel(y, x)
                    assert_array_equal(u.squeeze(), z.squeeze())
