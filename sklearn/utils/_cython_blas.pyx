# cython: boundscheck=False, wraparound=False, cdivision=True

from cython cimport floating

from scipy.linalg.cython_blas cimport sdot, ddot
from scipy.linalg.cython_blas cimport sasum, dasum
from scipy.linalg.cython_blas cimport saxpy, daxpy
from scipy.linalg.cython_blas cimport snrm2, dnrm2
from scipy.linalg.cython_blas cimport scopy, dcopy
from scipy.linalg.cython_blas cimport sscal, dscal
from scipy.linalg.cython_blas cimport sgemv, dgemv
from scipy.linalg.cython_blas cimport sger, dger
from scipy.linalg.cython_blas cimport sgemm, dgemm


################
# BLAS Level 1 #
################

cdef floating _xdot(int n, floating *x, int incx,
                    floating *y, int incy) nogil:
    """x.T.y"""
    if floating is float:
        return sdot(&n, x, &incx, y, &incy)
    else:
        return ddot(&n, x, &incx, y, &incy)


cpdef _xdot_memview(floating[::1] x, floating[::1] y):
    return _xdot(x.shape[0], &x[0], 1, &y[0], 1)


cdef floating _xasum(int n, floating *x, int incx) nogil:
    """sum(|x_i|)"""
    if floating is float:
        return sasum(&n, x, &incx)
    else:
        return dasum(&n, x, &incx)


cpdef _xasum_memview(floating[::1] x):
    return _xasum(x.shape[0], &x[0], 1)


cdef void _xaxpy(int n, floating alpha, floating *x, int incx,
                 floating *y, int incy) nogil:
    """y := alpha * x + y"""
    if floating is float:
        saxpy(&n, &alpha, x, &incx, y, &incy)
    else:
        daxpy(&n, &alpha, x, &incx, y, &incy)


cpdef _xaxpy_memview(floating alpha, floating[::1] x, floating[::1] y):
    _xaxpy(x.shape[0], alpha, &x[0], 1, &y[0], 1)


cdef floating _xnrm2(int n, floating *x, int incx) nogil:
    """sqrt(sum((x_i)^2))"""
    if floating is float:
        return snrm2(&n, x, &incx)
    else:
        return dnrm2(&n, x, &incx)


cpdef _xnrm2_memview(floating[::1] x):
    return _xnrm2(x.shape[0], &x[0], 1)


cdef void _xcopy(int n, floating *x, int incx, floating *y, int incy) nogil:
    """y := x"""
    if floating is float:
        scopy(&n, x, &incx, y, &incy)
    else:
        dcopy(&n, x, &incx, y, &incy)


cpdef _xcopy_memview(floating[::1] x, floating[::1] y):
    _xcopy(x.shape[0], &x[0], 1, &y[0], 1)


cdef void _xscal(int n, floating alpha, floating *x, int incx) nogil:
    """x := alpha * x"""
    if floating is float:
        sscal(&n, &alpha, x, &incx)
    else:
        dscal(&n, &alpha, x, &incx)


cpdef _xscal_memview(floating alpha, floating[::1] x):
    _xscal(x.shape[0], alpha, &x[0], 1)


################
# BLAS Level 2 #
################

cdef void _xgemv(BLAS_Order order, BLAS_Trans ta, int m, int n, floating alpha,
                 floating *A, int lda, floating *x, int incx,
                 floating beta, floating *y, int incy) nogil:
    """y := alpha * op(A).x + beta * y"""
    cdef char ta_ = ta
    if order == RowMajor:
        ta_ = NoTrans if ta == Trans else Trans 
        if floating is float:
            sgemv(&ta_, &n, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy)
        else:
            dgemv(&ta_, &n, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy)
    elif order == ColMajor:
        if floating is float:
            sgemv(&ta_, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)
        else:
            dgemv(&ta_, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)


cpdef _xgemv_memview(BLAS_Order order, BLAS_Trans ta, floating alpha,
                     floating[:, :] A, floating[::1] x, floating beta,
                     floating[::1] y):
    cdef:
        int m = A.shape[0]
        int n = A.shape[1]
        int lda = m if order == ColMajor else n

    _xgemv(order, ta, m, n, alpha, &A[0, 0], lda, &x[0], 1, beta, &y[0], 1)


cdef void _xger(BLAS_Order order, int m, int n, floating alpha, floating *x,
                int incx, floating *y, int incy, floating *A, int lda) nogil:
    """A := alpha * x.y.T + A"""
    if order == RowMajor:
        if floating is float:
            sger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
        else:
            dger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
    elif order == ColMajor:
        if floating is float:
            sger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)
        else:
            dger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)


cpdef _xger_memview(BLAS_Order order, floating alpha, floating[::1] x, floating[::] y,
                    floating[:, :] A):
    cdef:
        BLAS_Order order_ = ColMajor if order == ColMajor else RowMajor
        int m = A.shape[0]
        int n = A.shape[1]
        int lda = m if order == ColMajor else n
    
    _xger(order_, m, n, alpha, &x[0], 1, &y[0], 1, &A[0, 0], lda)


################
# BLAS Level 3 #
################

cdef void _xgemm(BLAS_Order order, BLAS_Trans ta, BLAS_Trans tb, int m, int n,
                 int k, floating alpha, floating *A, int lda, floating *B,
                 int ldb, floating beta, floating *C, int ldc) nogil:
    """C := alpha * op(A).op(B) + beta * C"""
    cdef:
        char ta_ = ta
        char tb_ = tb
    if order == RowMajor:
        if floating is float:
            sgemm(&tb_, &ta_, &n, &m, &k, &alpha, B,
                  &ldb, A, &lda, &beta, C, &ldc)
        else:
            dgemm(&tb_, &ta_, &n, &m, &k, &alpha, B,
                  &ldb, A, &lda, &beta, C, &ldc)
    elif order == ColMajor:
        if floating is float:
            sgemm(&ta_, &tb_, &m, &n, &k, &alpha, A,
                  &lda, B, &ldb, &beta, C, &ldc)
        else:
            dgemm(&ta_, &tb_, &m, &n, &k, &alpha, A,
                  &lda, B, &ldb, &beta, C, &ldc)


cpdef _xgemm_memview(BLAS_Order order, BLAS_Trans ta, BLAS_Trans tb,
                     floating alpha, floating[:, :] A, floating[:, :] B,
                     floating beta, floating[:, :] C):
    cdef:
        int m = A.shape[0] if ta == NoTrans else A.shape[1]
        int n = B.shape[1] if tb == NoTrans else B.shape[0]
        int k = A.shape[1] if ta == NoTrans else A.shape[0]
        int lda, ldb, ldc

    if order == ColMajor:
        lda = m if ta == NoTrans else k
        ldb = k if tb == NoTrans else n
        ldc = m
    else:
        lda = k if ta == NoTrans else m
        ldb = n if tb == NoTrans else k
        ldc = n

    _xgemm(order, ta, tb, m, n, k, alpha, &A[0, 0],
           lda, &B[0, 0], ldb, beta, &C[0, 0], ldc)