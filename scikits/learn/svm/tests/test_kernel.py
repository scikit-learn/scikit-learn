from numpy.testing import *
import numpy as N

set_local_path('../..')
from svm.kernel import *
restore_path()

class test_kernel(NumpyTestCase):
    def check_linear_kernel(self):
        kernel = LinearKernel()
        dot = N.dot
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x, dot), 4.)

    def check_polynomial_kernel(self):
        kernel = PolynomialKernel(degree=6, gamma=1.0, coef0=1.0)
        dot = N.dot
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x, dot), 15625.)

    def check_sigmoid_kernel(self):
        kernel = SigmoidKernel(gamma=0.2, coef0=0.3)
        dot = N.dot
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x, dot), 0.80049902)

    def check_rbf_kernel(self):
        kernel = RBFKernel(gamma=1.0)
        dot = N.dot
        x, y = N.array([2.]), N.array([3.])
        self.assertAlmostEqual(kernel(x, y, dot), N.exp(-1.))

    def check_custom_kernel(self):
        def f(x, y, dot):
            return 4 * dot(x, y)
        kernel = CustomKernel(f)
        def dot(x, y):
            return 2 * N.dot(x, y)
        x = N.array([2.])
        self.assertAlmostEqual(kernel(x, x, dot), 32.0)

if __name__ == '__main__':
    NumpyTest().run()
