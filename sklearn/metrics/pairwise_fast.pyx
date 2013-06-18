# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
import cython

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _chi2_kernel_fast(cython.floating[:,:] X,
                      cython.floating[:,:] Y,
                      cython.floating[:, :] result):
    cdef int i, j, k
    cdef int n_samples_X = X.shape[0]
    cdef int n_samples_Y = Y.shape[0]
    cdef int n_features = X.shape[1]
    cdef double res, nom, denom
    for i in xrange(n_samples_X):
        for j in xrange(n_samples_Y):
            res = 0
            for k in xrange(n_features):
                denom = (X[i, k] - Y[j, k])
                nom = (X[i, k] + Y[j, k])
                if nom != 0:
                    res  += denom * denom / nom
            result[i, j] = -res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _euclidean_distances_fast(cython.floating[::1] XX,
                              cython.floating[::1] YY,
                              cython.floating[:, ::1] distances):
    cdef int i, j
    cdef int n_samples_XX = XX.shape[0]
    cdef int n_samples_YY = YY.shape[0]
    cdef double res
    for i in xrange(n_samples_XX):
        for j in xrange(n_samples_YY):
            res = XX[i] + YY[j] - 2 * distances[i, j]
            if res < 0:
                res = 0
            distances[i, j] = res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _euclidean_distances_fast_sparse(cython.floating[:, :] XX,
                                     cython.floating[:, :] YY,
                                     cython.floating[:, :] distances):
    cdef int i, j
    cdef int n_samples_XX = XX.shape[0]
    cdef int n_samples_YY = YY.shape[1]
    cdef double res
    for i in xrange(n_samples_XX):
        for j in xrange(n_samples_YY):
            res = XX[i, 0] + YY[0, j] - 2 * distances[i, j]
            if res < 0:
                res = 0
            distances[i, j] = res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _euclidean_pdistances_fast(cython.floating[::1] XX,
                               cython.floating[:, ::1] distances):
    cdef int i, j
    cdef int n_samples_XX = XX.shape[0]
    cdef double res
    for i in xrange(n_samples_XX):
        for j in xrange(i):
            res = 2 * XX[i] - 2 * distances[i, j]
            if res < 0:
                res = 0
            distances[i, j] = distances[j, i] = res
        distances[i, i] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _euclidean_pdistances_fast_sparse(cython.floating[:, :] XX,
                                      cython.floating[:, :] distances):
    cdef int i, j
    cdef int n_samples_XX = XX.shape[0]
    cdef double res
    for i in xrange(n_samples_XX):
        for j in xrange(i):
            res = 2 * XX[i, 0] - 2 * distances[i, j]
            if res < 0:
                res = 0
            distances[i, j] = distances[j, i] = res
        distances[i, i] = 0

