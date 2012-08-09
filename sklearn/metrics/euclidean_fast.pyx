import numpy as np
cimport numpy as np

cimport cython

from cpython cimport bool
ctypedef np.int32_t INTEGER
ctypedef np.float64_t DOUBLE

cdef extern from "math.h":
    double sqrt "sqrt"(double x)

cdef extern from "cblas.h":
    double gemm "cblas_dgemm"(char STORAGE, char TRANSA, char TRANSB, int M,
                              int N, int K, double ALPHA, double *A, int LDA,
                              double *B, int LDB, double BETA, double *C,
                              int LDC)


@cython.boundscheck(False)
@cython.wraparound(False)
def dense_euclidean_distances(np.ndarray[DOUBLE, ndim=2, mode="c"] X,
                              np.ndarray[DOUBLE, ndim=2, mode="c"] Y,
                              np.ndarray[DOUBLE, ndim=2, mode="c"] out,
                              bool squared=False):
    cdef np.ndarray[DOUBLE, ndim=1] XX = np.empty(X.shape[0], dtype=np.float64)
    cdef np.ndarray[DOUBLE, ndim=1] YY = np.empty(Y.shape[0], dtype=np.float64)
    gemm(101, 111, 112, X.shape[0], Y.shape[0], Y.shape[1], 1.0,
         <DOUBLE *>X.data, X.strides[0] // sizeof(DOUBLE),
         <DOUBLE *>Y.data, Y.strides[0] // sizeof(DOUBLE), 0.0,
         <DOUBLE *>out.data, out.strides[0] // sizeof(DOUBLE))
    cdef unsigned int ii, jj
    for ii in range(X.shape[0]):
        XX[ii] = 0
        for jj in range(X.shape[1]):
            XX[ii] += X[ii, jj] ** 2
    for ii in range(Y.shape[0]):
        YY[ii] = 0
        for jj in range(Y.shape[1]):
            YY[ii] += Y[ii, jj] ** 2
    for ii in range(XX.shape[0]):
        for jj in range(YY.shape[0]):
            out[ii, jj] = XX[ii] + YY[jj] - 2 * out[ii, jj]
            out[ii, jj] = (out[ii, jj] if squared else sqrt(out[ii, jj])) \
                          if out[ii, jj] > 0 else 0
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
cdef sparse_sum_sq(n_samples, np.ndarray[DOUBLE, ndim=1] X_data,
                   np.ndarray[INTEGER, ndim=1] X_indices,
                   np.ndarray[INTEGER, ndim=1] X_indptr):
    cdef unsigned int ii, jj
    cdef DOUBLE tmp
    cdef np.ndarray[DOUBLE, ndim=1] out = np.empty(n_samples, dtype=np.float64)
    for ii in range(n_samples):
        tmp = 0.0
        for jj in range(X_indptr[ii], X_indptr[ii + 1]):
            tmp += X_data[jj] * X_data[jj]
        out[ii] = tmp
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def sparse_euclidean_distances(INTEGER X_rows, INTEGER Y_rows,
                               np.ndarray[DOUBLE, ndim=1] X_data,
                               np.ndarray[INTEGER, ndim=1] X_indices,
                               np.ndarray[INTEGER, ndim=1] X_indptr,
                               np.ndarray[DOUBLE, ndim=1] Y_data,
                               np.ndarray[INTEGER, ndim=1] Y_indices,
                               np.ndarray[INTEGER, ndim=1] Y_indptr,
                               np.ndarray[DOUBLE, ndim=2] XYt,
                               bool squared=False):
    cdef unsigned int ii, jj
    cdef np.ndarray[DOUBLE, ndim=1] XX = sparse_sum_sq(X_rows,
                                                       X_data, X_indices,
                                                       X_indptr)
    cdef np.ndarray[DOUBLE, ndim=1] YY = sparse_sum_sq(Y_rows,
                                                       Y_data, Y_indices,
                                                       Y_indptr)

    for ii in range(X_rows):
        for jj in range(Y_rows):
            XYt[ii, jj] = XX[ii] + YY[jj] - 2 * XYt[ii, jj]
            XYt[ii, jj] = (XYt[ii, jj] if squared else sqrt(XYt[ii, jj])) \
                          if XYt[ii, jj] > 0 else 0
    return XYt
