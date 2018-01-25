# Synopsis: Unified consistent BLAS wrappers
# Author: Elvis Dohmatob <gmdopp@gmail.com>
#
# Notes: A modern alternative would be to use scipy.linalg.cython_blas

from libc.math cimport fabs
from types cimport complexing
from utils cimport cabs, cabsf


cdef void fused_copy(int N, complexing *X, int incX, complexing *Y,
                     int incY) nogil:
    """Copy data from one buffer to another"""
    if complexing is float:
        scopy(N, X, incX, Y, incY)
    elif complexing is double:
        dcopy(N, X, incX, Y, incY)
    elif complexing is complex:
        zcopy(N, X, incX, Y, incY)
    else:
        ccopy(N, X, incX, Y, incY)


cdef void fused_scal(int N, complexing alpha, complexing *X, int incX) nogil:
    """In-place scaling of a buffer"""
    if complexing is float:
        sscal(N, alpha, X, incX)
    elif complexing is double:
        dscal(N, alpha, X, incX)
    elif complexing is complex:
        zscal(N, &alpha, X, incX)
    else:
        cscal(N, &alpha, X, incX)


cdef complexing fused_asum(int N, complexing *X, int incX) nogil:
    """Sum of absolute values of a buffer"""
    cdef complexing asumX
    cdef int k
    if complexing is float:
        return sasum(N, X, incX)
    elif complexing is double:
        return dasum(N, X, incX)
    else:
        asumX = 0.
        for k in range(N):
            k *= incX
            if complexing is complex:
                asumX += cabs(X[k])
            else:
                asumX += cabsf(X[k])
        return asumX


cdef complexing fused_dotc(int N, complexing *X, int incX, complexing *Y,
                           int incY) nogil:
    """Inner product of the complex conjugate of X with Y"""
    cdef complexing dotcXY
    if complexing is float:
        dotcXY = sdot(N, X, incX, Y, incY)
    elif complexing is double:
        dotcXY = ddot(N, X, incX, Y, incY)
    elif complexing is complex:
        zdotc(N, X, incX, Y, incY, &dotcXY)
    else:
        cdotc(N, X, incX, Y, incY, &dotcXY)
    return dotcXY


cdef complexing fused_dotu(int N, complexing *X, int incX, complexing *Y,
                           int incY) nogil:
    """Inner product of X with Y"""
    cdef complexing dotuXY
    if complexing is float:
        dotuXY = sdot(N, X, incX, Y, incY)
    elif complexing is double:
        dotuXY = ddot(N, X, incX, Y, incY)
    elif complexing is complex:
        zdotu(N, X, incX, Y, incY, &dotuXY)
    else:
        cdotu(N, X, incX, Y, incY, &dotuXY)
    return dotuXY


cdef complexing fused_nrm2(int N, complexing *X, int incX) nogil:
    """L2 norm"""
    if complexing is float:
        return snrm2(N, X, incX)
    elif complexing is double:
        return dnrm2(N, X, incX)
    elif complexing is complex:
        return dznrm2(N, X, incX)
    else:
        return scnrm2(N, X, incX)


cdef void fused_axpy(int N, complexing alpha, complexing *X, int incX,
                     complexing *Y, int incY) nogil:
    """Computes Y = Y + alpha * X"""
    if complexing is float:
        saxpy(N, alpha, X, incX, Y, incY)
    elif complexing is double:
        daxpy(N, alpha, X, incX, Y, incY)
    elif complexing is complex:
        zaxpy(N, &alpha, X, incX, Y, incY)
    else:
        caxpy(N, &alpha, X, incX, Y, incY)


cdef void fused_geru(CBLAS_ORDER Order, int M, int N, complexing alpha,
                     complexing *X, int incX, complexing *Y, int incY,
                     complexing *A, int lda) nogil:
    """Computes A = A + alpha * outer(X, Y)"""
    if complexing is float:
        sger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
    elif complexing is double:
        dger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
    elif complexing is complex:
        zgeru(Order, M, N, &alpha, X, incX, Y, incY, A, lda)
    else:
        cgeru(Order, M, N, &alpha, X, incX, Y, incY, A, lda)

