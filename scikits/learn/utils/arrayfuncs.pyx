"""
Small collection of auxiliary functions that operate on arrays

"""
cimport numpy as np
import  numpy as np


cdef extern from "cblas.h":
    double cblas_ddot(int N, double *X, int incX, double *Y, int incY)


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
