#cython: boundscheck=False
#cython: cdivision=True
#cython: wraparound=False
# cython: language_level=3
#
# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#         Lars Buitinck
#
# License: BSD 3 clause

import numpy as np
cimport numpy as np
from cython cimport floating
from libc.string cimport memset
from libc.math cimport sqrt

from ..utils._cython_blas cimport _asum


np.import_array()


def _chi2_kernel_fast(floating[:, :] X,
                      floating[:, :] Y,
                      floating[:, :] result):
    cdef np.npy_intp i, j, k
    cdef np.npy_intp n_samples_X = X.shape[0]
    cdef np.npy_intp n_samples_Y = Y.shape[0]
    cdef np.npy_intp n_features = X.shape[1]
    cdef double res, nom, denom

    with nogil:
        for i in range(n_samples_X):
            for j in range(n_samples_Y):
                res = 0
                for k in range(n_features):
                    denom = (X[i, k] - Y[j, k])
                    nom = (X[i, k] + Y[j, k])
                    if nom != 0:
                        res  += denom * denom / nom
                result[i, j] = -res


def _sparse_manhattan(floating[::1] X_data, int[:] X_indices, int[:] X_indptr,
                      floating[::1] Y_data, int[:] Y_indices, int[:] Y_indptr,
                      np.npy_intp n_features, double[:, ::1] D):
    """Pairwise L1 distances for CSR matrices.

    Usage:

    >>> D = np.zeros(X.shape[0], Y.shape[0])
    >>> sparse_manhattan(X.data, X.indices, X.indptr,
    ...                  Y.data, Y.indices, Y.indptr,
    ...                  X.shape[1], D)
    """
    cdef double[::1] row = np.empty(n_features)
    cdef np.npy_intp ix, iy, j

    with nogil:
        for ix in range(D.shape[0]):
            for iy in range(D.shape[1]):
                # Simple strategy: densify current row of X, then subtract the
                # corresponding row of Y.
                memset(&row[0], 0, n_features * sizeof(double))
                for j in range(X_indptr[ix], X_indptr[ix + 1]):
                    row[X_indices[j]] = X_data[j]
                for j in range(Y_indptr[iy], Y_indptr[iy + 1]):
                    row[Y_indices[j]] -= Y_data[j]

                D[ix, iy] = _asum(n_features, &row[0], 1)


cpdef _euclidean_distances_exact(floating[:, ::1] X,
                                 floating[:, ::1] Y,
                                 floating[:, ::1] D,
                                 bint squared):
    cdef:
        int n_samples_X = X.shape[0]
        int n_samples_Y = Y.shape[0]
        int n_features = X.shape[1]
        int i, j, k

    for i in range(n_samples_X):
        for j in range(n_samples_Y):
            D[i, j] = _euclidean_distance_dense_dense(
                &X[i, 0], &Y[j, 0], n_features, squared) 


cdef floating _euclidean_distance_dense_dense(floating *x,
                                              floating *y,
                                              int n_features,
                                              bint squared) nogil:
    """Euclidean distance between x dense and y dense"""
    cdef:
        int i
        floating tmp
        floating result = 0

    # We manually unroll the loop for better cache optimization.
    for i in range(n_features):
        tmp = x[i] - y[i]
        result += tmp * tmp

    return result if squared else sqrt(result)


cdef floating _euclidean_sparse_dense(floating *x_data,
                                      int *x_indices,
                                      floating *y,
                                      floating y_squared_norm,
                                      int nnz,
                                      bint squared) nogil:
    """Euclidean distance between x sparse and y dense"""
    cdef:
        int i
        floating yi
        floating tmp = 0.0
        floating result = 0.0

    for i in range(nnz):
        yi = y[x_indices[i]]
        tmp = x_data[i] - yi
        result += (tmp * tmp) - (yi * yi)

    result += y_squared_norm

    if result < 0: result = 0.0

    return result if squared else sqrt(result)


cdef floating _euclidean_sparse_sparse(floating *x_data,
                                       int *x_indices,
                                       floating *y_data,
                                       floating *y_indices,
                                       int x_nnz,
                                       int y_nnz,
                                       bint squared) nogil:
    """Euclidean distance between x sparse and y sparse"""
    cdef:
        int i
        floating yi
        floating tmp = 0.0
        floating result = 0.0

    for i in range(x_nnz):
        # yi = y[x_indices[i]]
        x_data[i] - y_data[j]          j = x_indices[i] if in y_indices
        tmp = x_data[i] - yi
        result += (tmp * tmp) - (yi * yi)

    result += y_squared_norm

    if result < 0: result = 0.0

    return result if squared else sqrt(result)