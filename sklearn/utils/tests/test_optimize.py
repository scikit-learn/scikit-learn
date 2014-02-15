""" Test the optimize module
"""
import numpy as np

from sklearn.utils.optimize import newton_cg


def test_newton_cg():
    rng = np.random.RandomState(0)
    A = rng.normal(size=(10, 10))

    def func(x):
        Ax = A.dot(x)
        return .5*(Ax).dot(Ax)

    def func_grad_hess(x):
        return func(x), A.T.dot(A.dot(x)), lambda x: A.T.dot(A.dot(x))

    x0 = np.ones(10)
    out = newton_cg(func_grad_hess, func, x0)
    np.testing.assert_array_almost_equal(np.zeros_like(x0), out)


