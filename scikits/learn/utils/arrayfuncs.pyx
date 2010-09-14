"""
Small collection of auxiliary functions that operate on arrays

"""
cimport numpy as np
import  numpy as np


cdef extern from "cblas.h":
   enum CBLAS_ORDER:
       CblasRowMajor=101
       CblasColMajor=102
   enum CBLAS_TRANSPOSE:
       CblasNoTrans=111
       CblasTrans=112
       CblasConjTrans=113
       AtlasConj=114
   enum CBLAS_UPLO:
       CblasUpper=121
       CblasLower=122
   enum CBLAS_DIAG:
       CblasNonUnit=131
       CblasUnit=132

   double cblas_ddot(int N, double *X, int incX, double *Y, int incY)

   void cblas_dtrsv(CBLAS_ORDER Order,  CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA,  CBLAS_DIAG Diag,
                 int N, double *A, int lda, double *X,
                 int incX)


ctypedef np.float64_t DOUBLE
ctypedef np.uint8_t   BOOL
# cython has no support for np.bool type. See
# http://www.mail-archive.com/cython-dev@codespeak.net/msg05913.html


def dot_over (
    np.ndarray[DOUBLE, ndim=2] X,
    np.ndarray[DOUBLE, ndim=1] y,
    np.ndarray[BOOL, ndim=1] mask,
    BOOL val,
    np.ndarray[DOUBLE, ndim=1] out):
    """
    Computes the dot product over a mask
    """

    cdef int i, j=0, N, incX, incY, incXY
    cdef double *X_ptr = <double *> X.data
    N     = X.shape[1]
    incX  = X.strides[1] / sizeof (double)
    incY  = y.strides[0] / sizeof (double)
    incXY = X.strides[0] / sizeof (double)

    for i in range(X.shape[0]):
        if mask[i] == val:
            out[j] = cblas_ddot (N, X_ptr + i*incXY, incX, <double *> y.data, incY)
            j += 1


def solve_triangular (
    np.ndarray[DOUBLE, ndim=2] X,
    np.ndarray[DOUBLE, ndim=1] y
    ):
    """
    Solves a triangular system (overwrites y)

    TODO: option for upper, lower, etc.
    """
    cdef int lda = <int> X.strides[0] / sizeof(double)

    cblas_dtrsv (CblasRowMajor, CblasLower, CblasNoTrans,
                 CblasNonUnit, <int> X.shape[0], <double *> X.data,
                 lda, <double *> y.data, 1);
