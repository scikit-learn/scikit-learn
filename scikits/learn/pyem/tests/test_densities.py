#! /usr/bin/env python
# Last Change: Thu Jul 12 05:00 PM 2007 J

# TODO:
#   - having "fake tests" to check that all mode (scalar, diag and full) are
#   executables
#   - having a dataset to check against

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from pyem.densities import gauss_den
import pyem.densities
restore_path()

#Optional:
set_local_path()
# import modules that are located in the same directory as this file.
from testcommon import DEF_DEC
restore_path()

class TestDensities(NumpyTestCase):
    def _generate_test_data_1d(self):
        self.va     = 2.0
        self.mu     = 1.0
        self.X      = N.linspace(-2, 2, 10)[:, N.newaxis]

        self.Yt     = N.array([0.02973257230591, 0.05512079811082,
            0.09257745306945, 0.14086453882683, 0.19418015562214,
            0.24250166773127, 0.27436665745048, 0.28122547107069,
            0.26114678964743, 0.21969564473386])

    def _generate_test_data_2d_diag(self):
        #============================
        # Small test in 2d (diagonal)
        #============================
        self.mu  = N.atleast_2d([-1.0, 2.0])
        self.va  = N.atleast_2d([2.0, 3.0])
        
        self.X  = N.zeros((10, 2))
        self.X[:,0] = N.linspace(-2, 2, 10)
        self.X[:,1] = N.linspace(-1, 3, 10)

        self.Yt  = N.array([0.01129091565384, 0.02025416837152,
            0.03081845516786, 0.03977576221540, 0.04354490552910,
            0.04043592581117, 0.03184994053539, 0.02127948225225,
            0.01205937178755, 0.00579694938623 ])


    def _generate_test_data_2d_full(self):
        #============================
        # Small test in 2d (full mat)
        #============================
        self.mu = N.array([[0.2, -1.0]])
        self.va = N.array([[1.2, 0.1], [0.1, 0.5]])
        X1      = N.linspace(-2, 2, 10)[:, N.newaxis]
        X2      = N.linspace(-3, 3, 10)[:, N.newaxis]
        self.X  = N.concatenate(([X1, X2]), 1)
        
        self.Yt = N.array([0.00096157109751, 0.01368908714856,
            0.07380823191162, 0.15072050533842, 
            0.11656739937861, 0.03414436965525,
            0.00378789836599, 0.00015915297541, 
            0.00000253261067, 0.00000001526368])

#=====================
# Basic accuracy tests
#=====================
class test_py_implementation(TestDensities):
    def _test(self, level, decimal = DEF_DEC):
        Y   = gauss_den(self.X, self.mu, self.va)
        assert_array_almost_equal(Y, self.Yt, decimal)

    def _test_log(self, level, decimal = DEF_DEC):
        Y   = gauss_den(self.X, self.mu, self.va, log = True)
        assert_array_almost_equal(N.exp(Y), self.Yt, decimal)

    def test_2d_diag(self, level = 0):
        self._generate_test_data_2d_diag()
        self._test(level)

    def test_2d_full(self, level = 0):
        self._generate_test_data_2d_full()
        self._test(level)
    
    def test_1d(self, level = 0):
        self._generate_test_data_1d()
        self._test(level)

    def test_2d_diag_log(self, level = 0):
        self._generate_test_data_2d_diag()
        self._test_log(level)

    def test_2d_full_log(self, level = 0):
        self._generate_test_data_2d_full()
        self._test_log(level)

    def test_1d_log(self, level = 0):
        self._generate_test_data_1d()
        self._test_log(level)

#=====================
# Basic speed tests
#=====================
class test_speed(NumpyTestCase):
    def __init__(self, *args, **kw):
        NumpyTestCase.__init__(self, *args, **kw)
        import sys
        import re
        try:
            a = open('/proc/cpuinfo').readlines()
            b = re.compile('cpu MHz')
            c = [i for i in a if b.match(i)]
            fcpu = float(c[0].split(':')[1])
            self.fcpu = fcpu * 1e6
            self.hascpu = True
        except:
            print "Could not read cpu frequency"
            self.hascpu = False
            self.fcpu = 0.

    def _prepare(self, n, d, mode):
        niter = 10
        x = 0.1 * N.random.randn(n, d)
        mu = 0.1 * N.random.randn(d)
        if mode == 'diag':
            va = 0.1 * N.random.randn(d) ** 2
        elif mode == 'full':
            a = N.random.randn(d, d)
            va = 0.1 * N.dot(a.T, a)
        st = self.measure("gauss_den(x, mu, va)", niter)
        return st / niter

    def _bench(self, n, d, mode):
        st = self._prepare(n, d, mode)
        print "%d dimension, %d samples, %s mode: %8.2f " % (d, n, mode, st)
        if self.hascpu:
            print "Cost per frame is %f; cost per sample is %f" % \
                    (st * self.fcpu / n, st * self.fcpu / n / d)
    
    def test1(self, level = 5):
        cls = self.__class__
        for n, d in [(1e5, 1), (1e5, 5), (1e5, 10), (1e5, 30), (1e4, 100)]:
            self._bench(n, d, 'diag')
        for n, d in [(1e4, 2), (1e4, 5), (1e4, 10), (5000, 40)]:
            self._bench(n, d, 'full')

#================
# Logsumexp tests
#================
class test_py_logsumexp(TestDensities):
    """Class to compare logsumexp vs naive implementation."""
    def test_underlow(self):
        """This function checks that logsumexp works as expected."""
        # We check wether naive implementation would underflow, to be sure we
        # are actually testing something here.
        errst = N.seterr(under='raise')
        try:
            try:
                a = N.array([[-1000]])
                self.naive_logsumexp(a)
                raise AssertionError("expected to catch underflow, we should"\
                                     "not be here")
            except FloatingPointError, e:
                print "Catching underflow, as expected"
            assert pyem.densities.logsumexp(a) == -1000.
            try:
                a = N.array([[-1000, -1000, -1000]])
                self.naive_logsumexp(a)
                raise AssertionError("expected to catch underflow, we should"\
                                     "not be here")
            except FloatingPointError, e:
                print "Catching underflow, as expected"
            assert_array_almost_equal(pyem.densities.logsumexp(a), 
                                      -998.90138771)
        finally:
            N.seterr(under=errst['under'])

    def naive_logsumexp(self, data):
        return N.log(N.sum(N.exp(data), 1)) 

    def test_1d(self):
        data = N.random.randn(1e1)[:, N.newaxis]
        mu = N.array([[-5], [-6]])
        va = N.array([[0.1], [0.1]])
        y = pyem.densities.multiple_gauss_den(data, mu, va, log = True)
        a1 = pyem.densities.logsumexp(y)
        a2 = self.naive_logsumexp(y)
        assert_array_equal(a1, a2)

    def test_2d_full(self):
        data = N.random.randn(1e1, 2)
        mu = N.array([[-3, -1], [3, 3]])
        va = N.array([[1.1, 0.4], [0.6, 0.8], [0.4, 0.2], [0.3, 0.9]])
        y = pyem.densities.multiple_gauss_den(data, mu, va, log = True)
        a1 = pyem.densities.logsumexp(y)
        a2 = self.naive_logsumexp(y)
        assert_array_almost_equal(a1, a2, DEF_DEC)

    def test_2d_diag(self):
        data = N.random.randn(1e1, 2)
        mu = N.array([[-3, -1], [3, 3]])
        va = N.array([[1.1, 0.4], [0.6, 0.8]])
        y = pyem.densities.multiple_gauss_den(data, mu, va, log = True)
        a1 = pyem.densities.logsumexp(y)
        a2 = self.naive_logsumexp(y)
        assert_array_almost_equal(a1, a2, DEF_DEC)

#=======================
# Test C implementation
#=======================
class test_c_implementation(TestDensities):
    def _test(self, level, decimal = DEF_DEC):
        try:
            from pyem._c_densities import gauss_den as c_gauss_den
            Y   = c_gauss_den(self.X, self.mu, self.va)
            assert_array_almost_equal(Y, self.Yt, decimal)
        except Exception, inst:
            print "Error while importing C implementation, not tested"
            print " -> (Import error was %s)" % inst 

    def test_1d(self, level = 0):
        self._generate_test_data_1d()
        self._test(level)

    def test_2d_diag(self, level = 0):
        self._generate_test_data_2d_diag()
        self._test(level)

    def test_2d_full(self, level = 0):
        self._generate_test_data_2d_full()
        self._test(level)

class test_gauss_ell(NumpyTestCase):
    def test_dim(self):
        pyem.densities.gauss_ell([0, 1], [1, 2.], [0, 1])
        try:
            pyem.densities.gauss_ell([0, 1], [1, 2.], [0, 2])
            raise AssertionError("this call should not succeed, bogus dim.")
        except ValueError, e:
            print "Call with bogus dim did not succeed, OK"


if __name__ == "__main__":
    NumpyTest().run()
