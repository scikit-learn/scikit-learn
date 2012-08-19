# Optimized euclidean distances computation between multiple points
#
# Author: Vlad Niculae <vlad@vene.ro>
# License: Simple BSD.
#
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

import numpy as np
cimport numpy as np

cimport cython
ctypedef np.int32_t INTEGER
ctypedef np.float64_t DTYPE_t

cdef extern from "math.h":
    double sqrt "sqrt"(double x)

cdef extern from "cblas.h":
    void dgemm "cblas_dgemm"(char storage, char transa, char transb, int m,
                             int n, int k, double alpha, double *A, int lda,
                             double *B, int ldb, double beta, double *C,
                             int ldc)
    void dsyrk "cblas_dsyrk"(char storage, char uplo, char trans, int n, int k,
                             double alpha, double *A, int lda, double beta,
                             double *C, int ldc)


cdef inline _euclidean_distances(int X_rows,
                                 int Y_rows,
                                 np.ndarray[DTYPE_t, ndim=1] X_norm_squared,
                                 np.ndarray[DTYPE_t, ndim=1] Y_norm_squared,
                                 np.ndarray[DTYPE_t, ndim=2] XYt,
                                 bint squared=False,
                                 bint sym=False):
    cdef DTYPE_t tmp
    for ii from 0 <= ii < X_rows:
        for jj from 0 <= jj < Y_rows:
            if ii == jj and sym:
                XYt[ii, jj] = 0
            else:
                tmp = XYt[jj, ii] if sym and ii < jj else XYt[ii, jj]
                tmp = X_norm_squared[ii] + Y_norm_squared[jj] - 2 * tmp
                XYt[ii, jj] = (tmp if squared else sqrt(tmp)) if tmp > 0 else 0


cdef inline dense_sum_sq(int X_rows, int X_cols,
                         np.ndarray[DTYPE_t, ndim=2] X,
                         np.ndarray[DTYPE_t, ndim=1] out):
    cdef unsigned int ii, jj
    for ii from 0 <= ii < X_rows:
        out[ii] = 0
        for jj from 0 <= jj < X_cols:
            out[ii] += X[ii, jj] * X[ii, jj]


cdef inline sparse_sum_sq(int X_rows, np.ndarray[DTYPE_t, ndim=1] X_data,
                   np.ndarray[int, ndim=1] X_indices,
                   np.ndarray[int, ndim=1] X_indptr,
                   np.ndarray[DTYPE_t, ndim=1] out):
    cdef unsigned int ii, jj
    for ii from 0 <= ii < X_rows:
        out[ii] = 0
        for jj from X_indptr[ii] <= jj < X_indptr[ii + 1]:
            out[ii] += X_data[jj] * X_data[jj]


def dense_euclidean_distances(int X_rows,
                              int Y_rows,
                              int X_cols,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] X,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] Y,
                              np.ndarray[DTYPE_t, ndim=1] X_norm_squared,
                              np.ndarray[DTYPE_t, ndim=1] Y_norm_squared,
                              bint X_norm_precomputed,
                              bint Y_norm_precomputed,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] out,
                              bint squared=False):
    if X_norm_squared is None:
        X_norm_squared = np.empty(X_rows, dtype=np.float64)
    if Y_norm_squared is None:
        Y_norm_squared = np.empty(Y_rows, dtype=np.float64)

    # 101 = c-order, 111 = no transpose, 112 = transpose (computes XY')
    dgemm(101, 111, 112, X_rows, Y_rows, X_cols, 1.0,
          <DTYPE_t *>X.data, X.strides[0] / sizeof(DTYPE_t),
          <DTYPE_t *>Y.data, Y.strides[0] / sizeof(DTYPE_t), 0.0,
          <DTYPE_t *>out.data, out.strides[0] / sizeof(DTYPE_t))

    if not X_norm_precomputed:
        dense_sum_sq(X_rows, X_cols, X, X_norm_squared)
    if not Y_norm_precomputed:
        dense_sum_sq(Y_rows, X_cols, Y, Y_norm_squared)
    _euclidean_distances(X_rows, Y_rows, X_norm_squared, Y_norm_squared, out,
                         squared)


def dense_euclidean_distances_sym(int X_rows,
                                  int X_cols,
                                  np.ndarray[DTYPE_t, ndim=2, mode="c"] X,
                                  np.ndarray[DTYPE_t, ndim=1] X_norm_squared,
                                  bint X_norm_precomputed,
                                  np.ndarray[DTYPE_t, ndim=2, mode="c"] out,
                                  bint squared=False):
    if X_norm_squared is None:
        X_norm_squared = np.empty(X_rows, dtype=np.float64)

    # 101 = c-order, 122 = lower triangular, 111 = no transpose (computes XX')
    # Attention, the upper triangle contains garbage.
    dsyrk(101, 122, 111, X_rows, X_cols, 1.0, <DTYPE_t *>X.data,
          X.strides[0] / sizeof(DTYPE_t), 0.0, <DTYPE_t *>out.data,
          out.strides[0] / sizeof(DTYPE_t))

    if not X_norm_precomputed:
        dense_sum_sq(X_rows, X_cols, X, X_norm_squared)

    _euclidean_distances(X_rows, X_rows, X_norm_squared, X_norm_squared, out,
                         squared, True)


# This function assumes the XYt parameter contains the precomputed outer dot
def sparse_euclidean_distances(int X_rows,
                               int Y_rows,
                               np.ndarray[DTYPE_t, ndim=1] X_data,
                               np.ndarray[int, ndim=1] X_indices,
                               np.ndarray[int, ndim=1] X_indptr,
                               np.ndarray[DTYPE_t, ndim=1] Y_data,
                               np.ndarray[int, ndim=1] Y_indices,
                               np.ndarray[int, ndim=1] Y_indptr,
                               np.ndarray[DTYPE_t, ndim=1] X_norm_squared,
                               np.ndarray[DTYPE_t, ndim=1] Y_norm_squared,
                               bint X_norm_precomputed,
                               bint Y_norm_precomputed,
                               np.ndarray[DTYPE_t, ndim=2] XYt,
                               bint squared=False):
    if X_norm_squared is None:
        X_norm_squared = np.empty(X_rows, dtype=np.float64)
    if Y_norm_squared is None:
        Y_norm_squared = np.empty(Y_rows, dtype=np.float64)

    if not X_norm_precomputed:
        sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, X_norm_squared)
    if not Y_norm_precomputed:
        sparse_sum_sq(Y_rows, Y_data, Y_indices, Y_indptr, Y_norm_squared)

    _euclidean_distances(X_rows, Y_rows, X_norm_squared, Y_norm_squared, XYt,
                         squared)


# This function assumes the XXt parameter contains the precomputed outer dot
def sparse_euclidean_distances_sym(int X_rows,
                                   np.ndarray[DTYPE_t, ndim=1] X_data,
                                   np.ndarray[int, ndim=1] X_indices,
                                   np.ndarray[int, ndim=1] X_indptr,
                                   np.ndarray[DTYPE_t, ndim=1] X_norm_squared,
                                   bint X_norm_precomputed,
                                   np.ndarray[DTYPE_t, ndim=2] XXt,
                                   bint squared=False):
    if X_norm_squared is None:
        X_norm_squared = np.empty(X_rows, dtype=np.float64)

    if not X_norm_precomputed:
        sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, X_norm_squared)

    _euclidean_distances(X_rows, X_rows, X_norm_squared, X_norm_squared, XXt,
                         squared, True)


# This function assumes the XYt parameter contains the precomputed outer dot
def sparse_dense_euclidean_distances(int X_rows,
                                     int X_cols,
                                     int Y_rows,
                                     np.ndarray[DTYPE_t, ndim=1] X_data,
                                     np.ndarray[int, ndim=1] X_indices,
                                     np.ndarray[int, ndim=1] X_indptr,
                                     np.ndarray[DTYPE_t, ndim=2, mode="c"] Y,
                                     np.ndarray[DTYPE_t, ndim=1]
                                         X_norm_squared,
                                     np.ndarray[DTYPE_t, ndim=1]
                                         Y_norm_squared,
                                     bint X_norm_precomputed,
                                     bint Y_norm_precomputed,
                                     np.ndarray[DTYPE_t, ndim=2] XYt,
                                     bint squared=False):
    if X_norm_squared is None:
        X_norm_squared = np.empty(X_rows, dtype=np.float64)
    if Y_norm_squared is None:
        Y_norm_squared = np.empty(Y_rows, dtype=np.float64)

    if not X_norm_precomputed:
        sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, X_norm_squared)
    if not Y_norm_precomputed:
        dense_sum_sq(Y_rows, X_cols, Y, Y_norm_squared)

    _euclidean_distances(X_rows, Y_rows, X_norm_squared, Y_norm_squared, XYt,
                         squared)
