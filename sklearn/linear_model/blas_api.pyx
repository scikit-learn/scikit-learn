from types cimport complexing


cdef void fused_copy(int N,
                     complexing *X,
                     int incX,
                     complexing *Y,
                     int incY) nogil:
    if complexing is float:
        scopy(N, X, incX, Y, incY)
    elif complexing is double:
        dcopy(N, X, incX, Y, incY)
    if complexing is complex:
        zcopy(N, X, incX, Y, incY)
    else:
        ccopy(N, X, incX, Y, incY)


cdef void fused_scal(int N,
                     complexing alpha,
                     complexing *X,
                     int incX) nogil:
    if complexing is float:
        sscal(N, alpha, X, incX)
    elif complexing is double:
        dscal(N, alpha, X, incX)
    elif complexing is complex:
        zscal(N, &alpha, X, incX)
    else:
        cscal(N, &alpha, X, incX)


cdef complexing fused_dotc(int N,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil:
    """Wraps sdot, ddot, cdotc_sub, and zdotc_sub"""
    cdef complexing z
    if complexing is float or complexing is double:
        return fused_dotu(N, X, incX, Y, incY)
    elif complexing is complex:
        zdotc(N, X, incX, Y, incY, &z)
        return z
    else:
        cdotc(N, X, incX, Y, incY, &z)
        return z


cdef complexing fused_dotu(int N,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil:
    """Wraps sdot, ddot, cdotu_sub, and zdotu_sub"""
    cdef complexing z
    if complexing is float:
        return sdot(N, X, incX, Y, incY)
    elif complexing is double:
        return ddot(N, X, incX, Y, incY)
    elif complexing is complex:
        zdotu(N, X, incX, Y, incY, &z)
        return z
    else:
        cdotu(N, X, incX, Y, incY, &z)
        return z


cdef double fused_nrm2(int N,
                       complexing *X,
                       int incX) nogil:
    if complexing is float:
        return snrm2(N, X, incX)
    elif complexing is double:
        return dnrm2(N, X, incX)
    elif complexing is complex:
        return dznrm2(N, X, incX)
    else:
        return scnrm2(N, X, incX)


cdef complexing fused_axpy(int N,
                           complexing alpha,
                           complexing *X,
                           int incX,
                           complexing *Y,
                           int incY) nogil:
    if complexing is float:
        saxpy(N, alpha, X, incX, Y, incY)
    elif complexing is double:
        daxpy(N, alpha, X, incX, Y, incY)
    elif complexing is complex:
        zaxpy(N, &alpha, X, incX, Y, incY)
    else:
        caxpy(N, &alpha, X, incX, Y, incY)


cdef void fused_geru(CBLAS_ORDER Order,
                     int M,
                     int N,
                     complexing alpha,
                     complexing *X,
                     int incX,
                     complexing *Y,
                     int incY,
                     complexing *A,
                     int lda) nogil:
    if complexing is float:
        sger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
    elif complexing is double:
        dger(Order, M, N, alpha, X, incX, Y, incY, A, lda)
    elif complexing is complex:
        zgeru(Order, M, N, &alpha, X, incX, Y, incY, A, lda)
    else:
        cgeru(Order, M, N, &alpha, X, incX, Y, incY, A, lda)



