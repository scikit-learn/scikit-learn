# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
import cython

np.import_array()

ctypedef np.float64_t double_t
ctypedef unsigned int uint_t
ctypedef float [:, :] float_array_2d_t 
ctypedef double [:, :] double_array_2d_t 

cdef fused floating_array_2d_t:
    float_array_2d_t
    double_array_2d_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _chi2_kernel_fast(floating_array_2d_t X,
                      floating_array_2d_t Y,
                      floating_array_2d_t result):
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


def _manhattan_distances_sparse(X, Y,
        np.ndarray[double_t, ndim=2] result):
    cdef:
        int i, r
        int n_x = len(X.data)
        int n_y = len(Y.data)
        int n_clusters = Y.shape[0]
        int n_documents = X.shape[0]
        
        np.ndarray[double_t, ndim=1] xdata = X.data
        np.ndarray[uint_t, ndim=1] xrow = X.row
        np.ndarray[uint_t, ndim=1] xcol = X.col
        np.ndarray[double_t, ndim=1] ydata = Y.data
        np.ndarray[uint_t, ndim=1] yrow = Y.row
        np.ndarray[uint_t, ndim=1] ycol = Y.col

    # for x, i, j in zip(X.data, X.row, X.col):
    for i in xrange(n_x):
        r = xrow[i] * n_clusters
        for y in xrange(n_clusters):
            r += y
            result[r, xcol[i]] = xdata[i]  # result[r, j] = x

    # for y, i, j in zip(Y.data, Y.row, Y.col):
    for i in xrange(n_y):
        for x in xrange(n_documents):
            r = x * n_clusters + yrow[i]
            result[r, ycol[i]] -= ydata[i]  # result[r, j] -= y

    result = np.abs(result)

