# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# Licence: BSD 3 clause

from libc.stdint cimport int64_t, int32_t
cimport numpy as np
import numpy as np
import cython


np.import_array()

ctypedef np.float64_t double_t
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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _manhattan_distances_sparse(X, Y):
    cdef:
        int i, r_, r, idx
        int n_x = len(X.data)
        int n_y = len(Y.data)
        int n_clusters = Y.shape[0]
        int n_documents = X.shape[0]

        np.ndarray[double_t, ndim=1] xdata = X.data
        np.ndarray[int, ndim=1] xrow = X.row
        np.ndarray[int, ndim=1] xcol = X.col
        np.ndarray[double_t, ndim=1] ydata = Y.data
        np.ndarray[int, ndim=1] yrow = Y.row
        np.ndarray[int, ndim=1] ycol = Y.col

        int64_t nnz = X.nnz * Y.shape[0] + Y.nnz * X.shape[0]
        np.ndarray[double_t, ndim=1] data = np.empty(nnz, dtype=np.float64)
        np.ndarray[int32_t, ndim=1] rows = np.empty(nnz, dtype=np.int32)
        np.ndarray[int32_t, ndim=1] cols = np.empty(nnz, dtype=np.int32)

    idx = 0
    for i in xrange(n_x):
        r_ = xrow[i] * n_clusters
        for y in xrange(n_clusters):
            r = r_ + y
            rows[idx] = r
            cols[idx] = xcol[i]
            data[idx] = xdata[i]
            idx += 1

    for i in xrange(n_y):
        for x in xrange(n_documents):
            r = x * n_clusters + yrow[i]
            rows[idx] = r
            cols[idx] = ycol[i]
            data[idx] = -ydata[i]
            idx += 1
    return data, rows, cols