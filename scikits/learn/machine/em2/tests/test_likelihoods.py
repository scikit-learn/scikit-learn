import numpy as np
from numpy.testing import assert_array_almost_equal
from unittest import TestCase

from scikits.learn.machine.em2.likelihoods import logsumexp, mnormalik, normalik

class TestLogsumexp(TestCase):
    """Class to compare logsumexp vs naive implementation."""
    def test_underlow(self):
        """This function checks that logsumexp works as expected."""
        # We check wether naive implementation would underflow, to be sure we
        # are actually testing something here.
        errst = np.seterr(under='raise')
        try:
            try:
                a = np.array([[-1000.]])
                self.naive_logsumexp(a)
                raise AssertionError("expected to catch underflow, we should"\
                                     "not be here")
            except FloatingPointError:
                pass
            assert logsumexp(a) == -1000.
            try:
                a = np.array([[-1000., -1000, -1000]])
                self.naive_logsumexp(a)
                raise AssertionError("expected to catch underflow, we should"\
                                     "not be here")
            except FloatingPointError, e:
                pass
            assert_array_almost_equal(logsumexp(a), 
                                      -998.90138771)
        finally:
            np.seterr(under=errst['under'])

    def naive_logsumexp(self, data):
        return np.log(np.sum(np.exp(data), 1)) 

    def test_1d(self):
        data = np.random.randn(1e1)[:, np.newaxis]
        mu = np.array([[-5], [-6]])
        va = np.array([[0.1], [0.1]])
        y = mnormalik(data, mu, va, log=True)
        a1 = logsumexp(y)
        a2 = self.naive_logsumexp(y)
        assert_array_almost_equal(a1, a2)

    def test_2d_diag(self):
        data = np.random.randn(1e1, 2)
        mu = np.array([[-3, -1], [3, 3]])
        va = np.array([[1.1, 0.4], [0.6, 0.8]])
        y = mnormalik(data, mu, va, log=True)
        a1 = logsumexp(y)
        a2 = self.naive_logsumexp(y)
        assert_array_almost_equal(a1, a2)

X = {}
X[1] = np.linspace(-2, 2, 10)[:, np.newaxis]
X[2] = np.concatenate((np.linspace(-2, 2, 10)[:, np.newaxis], 
                       np.linspace(-1, 3, 10)[:, np.newaxis]), axis=1)
X[3] = np.concatenate((np.linspace(-2, 2, 10)[:, np.newaxis], 
                       np.linspace(-3, 3, 10)[:, np.newaxis]), axis=1)

MU = {}
VA = {}
MU[1] = np.array([1.])
VA[1] = np.array([2.])
MU[2] = np.array([-1., 2.])
VA[2] = np.array([2., 3.])
MU[3] = np.array([0.2, -1.0])
VA[3] = np.array([[1.2, 0.1], [0.1, 0.5]])

Y = {}
Y[1] = np.array([0.02973257230591, 0.05512079811082, 0.09257745306945,
    0.14086453882683, 0.19418015562214, 0.24250166773127, 0.27436665745048,
    0.28122547107069, 0.26114678964743, 0.21969564473386])
Y[2] = np.array([0.01129091565384, 0.02025416837152, 0.03081845516786,
    0.03977576221540, 0.04354490552910, 0.04043592581117, 0.03184994053539,
    0.02127948225225, 0.01205937178755, 0.00579694938623 ])
Y[3] = np.array([0.00096157109751, 0.01368908714856, 0.07380823191162,
    0.15072050533842, 0.11656739937861, 0.03414436965525, 0.00378789836599,
    0.00015915297541, 0.00000253261067, 0.00000001526368])

class TestNormalLikelihood(TestCase):
    def _test(self, id, decimal=10, log=False):
        #from scikits.learn.machine.em.densities import gauss_den
        y = normalik(X[id], MU[id], VA[id], log=log)
        #assert_array_almost_equal(gauss_den(X[id], MU[id], VA[id], log=log), y)
        if log:
            assert_array_almost_equal(y, np.log(Y[id]), decimal)
        else:
            assert_array_almost_equal(y, Y[id], decimal)

    def test_oned(self):
        return self._test(1)

    def test_twod(self):
        return self._test(2)

    def test_log_one1d(self):
        return self._test(1, log=True)

    def test_log_one2d(self):
        return self._test(2, log=True)

class TestMnormalLikelihood(TestCase):
    def _test(self, n, d, k, decimal=10, log=False):
        x = np.random.randn(n, d)
        mu = np.random.randn(k, d)
        va = np.random.randn(k, d) ** 2
        yr = np.empty((n, k))

        y = mnormalik(x, mu, va, log=log)
        for c in range(k):
            yr[:, c] = normalik(x, mu[c], va[c], log=log)

        assert_array_almost_equal(y, yr, decimal)

    def test_oned(self):
        return self._test(1e3, 1, 2)

    def test_twod(self):
        return self._test(1e3, 2, 2)

    def test_log_one1d(self):
        return self._test(1e3, 1, 2, log=True)

    def test_log_one2d(self):
        return self._test(1e3, 2, 2, log=True)


