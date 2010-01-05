#! /usr/bin/env python
# Last Change: Fri Sep 29 06:00 PM 2006 J

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from pyem.densities import gauss_den
from pyem._c_densities import gauss_den as c_gauss_den

restore_path()

#Optional:
set_local_path()
# import modules that are located in the same directory as this file.
restore_path()

class test_densities(NumpyTestCase):
    def _generate_test_data_1d(self):
        self.va     = 2.0
        self.mu     = 1.0
        self.X      = N.linspace(-2, 2, 10)[:, N.newaxis]

        self.Yt     = N.array([0.02973257230591, 0.05512079811082, 0.09257745306945, 
                0.14086453882683,
                0.19418015562214, 0.24250166773127, 0.27436665745048, 0.28122547107069,
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

        self.Yt  = N.array([0.01129091565384, 0.02025416837152, 0.03081845516786, 
                0.03977576221540, 0.04354490552910, 0.04043592581117, 
                0.03184994053539, 0.02127948225225, 0.01205937178755, 
                0.00579694938623 ])


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

    def _check_py(self, level, decimal = 12):
        Y   = gauss_den(self.X, self.mu, self.va)
        assert_array_almost_equal(Y, self.Yt, decimal)

    def _check_c(self, level, decimal = 12):
        Y   = c_gauss_den(self.X, self.mu, self.va)
        assert_array_almost_equal(Y, self.Yt, decimal)

    def check_py_1d(self, level = 1):
        self._generate_test_data_1d()
        self._check_py(level)

    def check_py_2d_diag(self, level = 1):
        self._generate_test_data_2d_diag()
        self._check_py(level)

    def check_py_2d_full(self, level = 1):
        self._generate_test_data_2d_full()
        self._check_py(level)

    def check_c_1d(self, level = 1):
        self._generate_test_data_1d()
        self._check_c(level)

    def check_c_2d_diag(self, level = 1):
        self._generate_test_data_2d_diag()
        self._check_c(level)

    def check_c_2d_full(self, level = 1):
        self._generate_test_data_2d_full()
        self._check_c(level)

if __name__ == "__main__":
    NumpyTest().run()

# def generate_test_data(n, d, mode = 'diag', file='test.dat'):
#     """Generate a set of data of dimension d, with n frames,
#     that is input data, mean, var and output of gden, so that
#     other implementations can be tested against"""
#     mu  = randn(1, d)
#     if mode == 'diag':
#         va  = abs(randn(1, d))
#     elif mode == 'full':
#         va  = randn(d, d)
#         va  = dot(va, va.transpose())
# 
#     input   = randn(n, d)
#     output  = gauss_den(input, mu, va)
# 
#     import tables
#     h5file  = tables.openFile(file, "w")
# 
#     h5file.createArray(h5file.root, 'input', input)
#     h5file.createArray(h5file.root, 'mu', mu)
#     h5file.createArray(h5file.root, 'va', va)
#     h5file.createArray(h5file.root, 'output', output)
# 
#     h5file.close()
# 
# def test_gauss_den():
#     """"""
#     # import tables
#     # import numpy as N
#     # 
#     # filename    = 'dendata.h5'
# 
#     # # # Dimension 1
#     # # d   = 1
#     # # mu  = 1.0
#     # # va  = 2.0
# 
#     # # X   = randn(1e3, 1)
# 
#     # # Y   = gauss_den(X, mu, va)
# 
#     # # h5file      = tables.openFile(filename, "w")
# 
#     # # h5file.createArray(h5file.root, 'X', X)
#     # # h5file.createArray(h5file.root, 'mu', mu)
#     # # h5file.createArray(h5file.root, 'va', va)
#     # # h5file.createArray(h5file.root, 'Y', Y)
# 
#     # # h5file.close()
# 
#     # # # Dimension 2, diag
#     # # d   = 2
#     # # mu  = N.array([1.0, -2.0])
#     # # va  = N.array([1.0, 2.0])
# 
#     # # X   = randn(1e3, 2)
# 
#     # # Y   = gauss_den(X, mu, va)
# 
#     # # h5file      = tables.openFile(filename, "w")
# 
#     # # h5file.createArray(h5file.root, 'X', X)
#     # # h5file.createArray(h5file.root, 'mu', mu)
#     # # h5file.createArray(h5file.root, 'va', va)
#     # # h5file.createArray(h5file.root, 'Y', Y)
# 
#     # # Dimension 2, full
#     # d   = 2
#     # mu  = N.array([[0.2, -1.0]])
#     # va  = N.array([[1.2, 0.1], [0.1, 0.5]])
# 
#     # X   = randn(1e3, 2)
# 
#     # Y   = gauss_den(X, mu, va)
# 
#     # h5file      = tables.openFile(filename, "w")
# 
#     # h5file.createArray(h5file.root, 'X', X)
#     # h5file.createArray(h5file.root, 'mu', mu)
#     # h5file.createArray(h5file.root, 'va', va)
#     # h5file.createArray(h5file.root, 'Y', Y)
# 
#     # h5file.close()
# 
