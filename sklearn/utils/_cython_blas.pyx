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

cdef void _xgemv(char layout, char ta, int m, int n, floating alpha,
                 floating *A, int lda, floating *x, int incx,
                 floating beta, floating *y, int incy) nogil:
    """y := alpha * op(A).x + beta * y"""
    if layout == 'C':
        ta = 'n' if ta == 't' else 't' 
        if floating is float:
            sgemv(&ta, &n, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy)
        else:
            dgemv(&ta, &n, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy)
    elif layout == 'F':
        if floating is float:
            sgemv(&ta, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)
        else:
            dgemv(&ta, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)


cpdef _xgemv_memview(layout, ta, floating alpha, floating[:, :] A,
                     floating[::1] x, floating beta, floating[::1] y):
    cdef:
        char layout_ = 'F' if layout == 'F' else 'C'
        char ta_ = 'n' if ta == 'n' else 't'
        int m = A.shape[0]
        int n = A.shape[1]
        int lda = m if layout == 'F' else n

    _xgemv(layout_, ta_, m, n, alpha, &A[0, 0], lda,
           &x[0], 1, beta, &y[0], 1)


cdef void _xger(char layout, int m, int n, floating alpha, floating *x,
                int incx, floating *y, int incy, floating *A, int lda) nogil:
    """A := alpha * x.y.T + A"""
    if layout == 'C':
        if floating is float:
            sger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
        else:
            dger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
    elif layout == 'F':
        if floating is float:
            sger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)
        else:
            dger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)


cpdef _xger_memview(layout, floating alpha, floating[::1] x, floating[::] y,
                    floating[:, :] A):
    cdef:
        char layout_ = 'F' if layout == 'F' else 'C'
        int m = A.shape[0]
        int n = A.shape[1]
        int lda = m if layout[0] == 'F' else n
    
    _xger(layout_, m, n, alpha, &x[0], 1, &y[0], 1, &A[0, 0], lda)


################
# BLAS Level 3 #
################

cdef void _xgemm(char layout, char ta, char tb, int m, int n, int k,
                 floating alpha, floating *A, int lda, floating *B, int ldb,
                 floating beta, floating *C, int ldc) nogil:
    """C := alpha * op(A).op(B) + beta * C"""
    if layout == 'C':
        if floating is float:
            sgemm(&tb, &ta, &n, &m, &k, &alpha, B, &ldb, A, &lda, &beta, C, &ldc)
        else:
            dgemm(&tb, &ta, &n, &m, &k, &alpha, B, &ldb, A, &lda, &beta, C, &ldc)
    elif layout == 'F':
        if floating is float:
            sgemm(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc)
        else:
            dgemm(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc)


cpdef _xgemm_memview(layout, ta, tb, floating alpha, floating[:, :] A,
                     floating[:, :] B, floating beta, floating[:, :] C):
    cdef:
        char layout_ = 'F' if layout == 'F' else 'C'
        char ta_ = 'n' if ta == 'n' else 't'
        char tb_ = 'n' if tb == 'n' else 't'
        int m = A.shape[0] if ta[0] == 'n' else A.shape[1]
        int n = B.shape[1] if tb[0] == 'n' else B.shape[0]
        int k = A.shape[1] if ta[0] == 'n' else A.shape[0]
        int lda, ldb, ldc

    if layout == 'F':
        lda = m if ta == 'n' else k
        ldb = k if tb == 'n' else n
        ldc = m
    else:
        lda = k if ta == 'n' else m
        ldb = n if tb == 'n' else k
        ldc = n

    _xgemm(layout_, ta_, tb_, m, n, k, alpha,
           &A[0, 0], lda, &B[0, 0], ldb, beta, &C[0, 0], ldc)