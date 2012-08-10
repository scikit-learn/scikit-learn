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
    void gemm "cblas_dgemm"(char storage, char transa, char transb, int m,
                            int n, int k, double alpha, double *A, int lda,
                            double *B, int ldb, double beta, double *C,
                            int ldc)
    void syrk "cblas_dsyrk"(char storage, char uplo, char trans, int n, int k, 
                            double alpha, double *A, int lda, double beta, 
                            double *C, int ldc)


cdef inline _euclidean_distances(int X_rows, int Y_rows,
                                 np.ndarray[DTYPE_t, ndim=1] XX,
                                 np.ndarray[DTYPE_t, ndim=1] YY,
                                 np.ndarray[DTYPE_t, ndim=2] XYt,
                                 bint squared=False, bint sym=False):
    cdef DTYPE_t tmp
    for ii from 0 <= ii < X_rows:
        for jj from 0 <= jj < Y_rows:
            if ii == jj and sym:  # maybe break sym into a separate function
                XYt[ii, jj] = 0
            else:
                tmp = XYt[jj, ii] if sym and ii < jj else XYt[ii, jj]
                tmp = XX[ii] + YY[jj] - 2 * tmp
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


def dense_euclidean_distances(int X_rows, int Y_rows, int X_cols,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] X,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] Y,
                              np.ndarray[DTYPE_t, ndim=2, mode="c"] out,
                              bint squared=False):
    cdef np.ndarray[DTYPE_t, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    cdef np.ndarray[DTYPE_t, ndim=1] YY = np.empty(Y_rows, dtype=np.float64)
    gemm(101, 111, 112, X_rows, Y_rows, X_cols, 1.0,
         <DTYPE_t *>X.data, X.strides[0] / sizeof(DTYPE_t),
         <DTYPE_t *>Y.data, Y.strides[0] / sizeof(DTYPE_t), 0.0,
         <DTYPE_t *>out.data, out.strides[0] / sizeof(DTYPE_t))
    dense_sum_sq(X_rows, X_cols, X, XX)
    dense_sum_sq(Y_rows, X_cols, Y, YY)
    _euclidean_distances(X_rows, Y_rows, XX, YY, out, squared)


def dense_euclidean_distances_sym(int X_rows, int X_cols,
                                  np.ndarray[DTYPE_t, ndim=2, mode="c"] X,
                                  np.ndarray[DTYPE_t, ndim=2, mode="c"] out,
                                  bint squared=False):
    cdef np.ndarray[DTYPE_t, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    syrk(101, 122, 111, X_rows, X_cols, 1.0, <DTYPE_t *>X.data,
         X.strides[0] / sizeof(DTYPE_t), 0.0, <DTYPE_t *>out.data,
         out.strides[0] / sizeof(DTYPE_t))
    dense_sum_sq(X_rows, X_cols, X, XX)
    _euclidean_distances(X_rows, X_rows, XX, XX, out, squared, True)


def sparse_euclidean_distances(int X_rows, int Y_rows,
                               np.ndarray[DTYPE_t, ndim=1] X_data,
                               np.ndarray[int, ndim=1] X_indices,
                               np.ndarray[int, ndim=1] X_indptr,
                               np.ndarray[DTYPE_t, ndim=1] Y_data,
                               np.ndarray[int, ndim=1] Y_indices,
                               np.ndarray[int, ndim=1] Y_indptr,
                               np.ndarray[DTYPE_t, ndim=2] XYt,
                               bint squared=False):
    cdef np.ndarray[DTYPE_t, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    cdef np.ndarray[DTYPE_t, ndim=1] YY = np.empty(Y_rows, dtype=np.float64)
    sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, XX)
    sparse_sum_sq(Y_rows, Y_data, Y_indices, Y_indptr, YY)
    _euclidean_distances(X_rows, Y_rows, XX, YY, XYt, squared)


def sparse_euclidean_distances_sym(int X_rows,
                                   np.ndarray[DTYPE_t, ndim=1] X_data,
                                   np.ndarray[int, ndim=1] X_indices,
                                   np.ndarray[int, ndim=1] X_indptr,
                                   np.ndarray[DTYPE_t, ndim=2] XXt,
                                   bint squared=False):
    cdef np.ndarray[DTYPE_t, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, XX)
    _euclidean_distances(X_rows, X_rows, XX, XX, XXt, squared, True)
