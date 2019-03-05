# Synopsis: Unified consistent BLAS wrappers
# Author: Elvis Dohmatob <gmdopp@gmail.com>
#
# Notes: A modern alternative would be to use scipy.linalg.cython_blas

from libc.math cimport fabs
from cython cimport floating


cdef void fused_copy(int N, floating *X, int incX, floating *Y,
                     int incY) nogil:
    """Copy data from one buffer to another"""
    if floating is float:
        scopy(N, X, incX, Y, incY)
    else:
        dcopy(N, X, incX, Y, incY)


cdef void fused_scal(int N, floating alpha, floating *X, int incX) nogil:
    """In-place scaling of a buffer"""
    if floating is float:
        sscal(N, alpha, X, incX)
    else:
        dscal(N, alpha, X, incX)


cdef floating fused_asum(int N, floating *X, int incX) nogil:
    """Sum of absolute values of a buffer"""
    if floating is float:
        return sasum(N, X, incX)
    else:
        return dasum(N, X, incX)


cdef floating fused_dot(int N, floating *X, int incX, floating *Y,
                        int incY) nogil:
    """Inner product of X with Y"""
    if floating is float:
        return sdot(N, X, incX, Y, incY)
    else:
        return ddot(N, X, incX, Y, incY)


cdef floating fused_nrm2(int N, floating *X, int incX) nogil:
    """L2 norm"""
    if floating is float:
        return snrm2(N, X, incX)
    else:
        return dnrm2(N, X, incX)


cdef void fused_axpy(int N, floating alpha, floating *X, int incX,
                     floating *Y, int incY) nogil:
    """Computes Y = Y + alpha * X"""
    if floating is float:
        saxpy(N, alpha, X, incX, Y, incY)
    else:
        daxpy(N, alpha, X, incX, Y, incY)

cdef void fused_ger(CBLAS_ORDER Order, int M, int N, floating alpha,
                    floating *X, int incX, floating *Y, int incY,
                    floating *A, int lda) nogil:
    """Computes A = A + alpha * outer(X, Y)"""
    if floating is float:
        sger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
    else:
        dger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
