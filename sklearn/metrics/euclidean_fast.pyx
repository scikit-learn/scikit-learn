import numpy as np
cimport numpy as np

cimport cython
ctypedef np.int32_t INTEGER
ctypedef np.float64_t DOUBLE

cdef extern from "math.h":
    double sqrt "sqrt"(double x)

cdef extern from "cblas.h":
    double gemm "cblas_dgemm"(char storage, char transa, char transb, int m,
                              int n, int k, double alpha, double *A, int lda,
                              double *B, int ldb, double beta, double *C,
                              int ldc)
    int syrk "cblas_dsyrk"(char storage, char uplo, char trans, int n, int k, 
                           double alpha, double *A, int lda, double beta, 
                           double *C, int ldc)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline _euclidean_distances(INTEGER X_rows, INTEGER Y_rows,
                                 np.ndarray[DOUBLE, ndim=1] XX,
                                 np.ndarray[DOUBLE, ndim=1] YY,
                                 np.ndarray[DOUBLE, ndim=2] XYt,
                                 bint squared=False, bint self=False):
    cdef DOUBLE tmp
    for ii in range(X_rows):
        for jj in range(Y_rows):
            if ii == jj and self:
                XYt[ii, jj] = 0
            else:
                tmp = XYt[jj, ii] if self and ii < jj else XYt[ii, jj] 
                tmp = XX[ii] + YY[jj] - 2 * tmp
                XYt[ii, jj] = (tmp if squared else sqrt(tmp)) if tmp > 0 else 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline dense_sum_sq(INTEGER X_rows, INTEGER X_cols,
                         np.ndarray[DOUBLE, ndim=2] X,
                         np.ndarray[DOUBLE, ndim=1] out):
    cdef unsigned int ii, jj
    for ii in range(X_rows):
        out[ii] = 0
        for jj in range(X_cols):
            out[ii] += X[ii, jj] * X[ii, jj]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline sparse_sum_sq(INTEGER n_samples, np.ndarray[DOUBLE, ndim=1] X_data,
                   np.ndarray[INTEGER, ndim=1] X_indices,
                   np.ndarray[INTEGER, ndim=1] X_indptr,
                   np.ndarray[DOUBLE, ndim=1] out):
    cdef unsigned int ii, jj
    for ii in range(n_samples):
        out[ii] = 0
        for jj in range(X_indptr[ii], X_indptr[ii + 1]):
            out[ii] += X_data[jj] * X_data[jj]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dense_euclidean_distances(INTEGER X_rows, INTEGER Y_rows, INTEGER X_cols, 
                              np.ndarray[DOUBLE, ndim=2, mode="c"] X,
                              np.ndarray[DOUBLE, ndim=2, mode="c"] Y,
                              np.ndarray[DOUBLE, ndim=2, mode="c"] out,
                              bint squared=False):
    cdef np.ndarray[DOUBLE, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    cdef np.ndarray[DOUBLE, ndim=1] YY = np.empty(Y_rows, dtype=np.float64)
    gemm(101, 111, 112, X_rows, Y_rows, X_cols, 1.0,
         <DOUBLE *>X.data, X.strides[0] / sizeof(DOUBLE),
         <DOUBLE *>Y.data, Y.strides[0] / sizeof(DOUBLE), 0.0,
         <DOUBLE *>out.data, out.strides[0] / sizeof(DOUBLE))
    dense_sum_sq(X_rows, X_cols, X, XX)
    dense_sum_sq(Y_rows, X_cols, Y, YY)
    _euclidean_distances(X_rows, Y_rows, XX, YY, out, squared)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dense_euclidean_distances_self(INTEGER X_rows, INTEGER X_cols, 
                                   np.ndarray[DOUBLE, ndim=2, mode="c"] X,
                                   np.ndarray[DOUBLE, ndim=2, mode="c"] out,
                                   bint squared=False):
    cdef np.ndarray[DOUBLE, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    syrk(101, 122, 111, X_rows, X_cols, 1.0, <DOUBLE *>X.data,
         X.strides[0] / sizeof(DOUBLE), 0.0, <DOUBLE *>out.data,
         out.strides[0] / sizeof(DOUBLE))
    dense_sum_sq(X_rows, X_cols, X, XX)
    _euclidean_distances(X_rows, X_rows, XX, XX, out, squared, True)


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
                               bint squared=False):
    cdef np.ndarray[DOUBLE, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    cdef np.ndarray[DOUBLE, ndim=1] YY = np.empty(Y_rows, dtype=np.float64)
    sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, XX)
    sparse_sum_sq(Y_rows, Y_data, Y_indices, Y_indptr, YY)
    _euclidean_distances(X_rows, Y_rows, XX, YY, XYt, squared)



@cython.boundscheck(False)
@cython.wraparound(False)
def sparse_euclidean_distances_self(INTEGER X_rows,
                                    np.ndarray[DOUBLE, ndim=1] X_data,
                                    np.ndarray[INTEGER, ndim=1] X_indices,
                                    np.ndarray[INTEGER, ndim=1] X_indptr,
                                    np.ndarray[DOUBLE, ndim=2] XXt,
                                    bint squared=False):
    cdef np.ndarray[DOUBLE, ndim=1] XX = np.empty(X_rows, dtype=np.float64)
    sparse_sum_sq(X_rows, X_data, X_indices, X_indptr, XX)
    _euclidean_distances(X_rows, X_rows, XX, XX, XXt, squared)
