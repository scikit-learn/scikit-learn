from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.kernel import *
restore_path()

def kernelfunc(x, y):
    return 8 * N.dot(x, y.T)

class test_kernel(NumpyTestCase):
    def check_linear_kernel(self):
        kernel = LinearKernel()
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 4.)

    def check_polynomial_kernel(self):
        kernel = PolynomialKernel(degree=6, gamma=1.0, coef0=1.0)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 15625.)

    def check_sigmoid_kernel(self):
        kernel = SigmoidKernel(gamma=0.2, coef0=0.3)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 0.80049902)

    def check_rbf_kernel(self):
        kernel = RBFKernel(gamma=1.0)
        x, y = N.array([2.]), N.array([3.])
        self.assertAlmostEqual(kernel(x, y), N.exp(-1.))

    def check_custom_kernel(self):
        kernel = CustomKernel(kernelfunc)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x), 32.0)

    def check_multidim_input(self):
        kernels = [
            LinearKernel(),
            PolynomialKernel(degree=6, gamma=1.0, coef0=1.0),
            SigmoidKernel(gamma=0.2, coef0=0.3),
            RBFKernel(gamma=1.0),
            CustomKernel(kernelfunc)
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

if __name__ == '__main__':
    NumpyTest().run()
