from cython cimport floating

from scipy.linalg.cython_blas cimport sdot, ddot
from scipy.linalg.cython_blas cimport sasum, dasum
from scipy.linalg.cython_blas cimport saxpy, daxpy
from scipy.linalg.cython_blas cimport snrm2, dnrm2
from scipy.linalg.cython_blas cimport scopy, dcopy
from scipy.linalg.cython_blas cimport sscal, dscal
from scipy.linalg.cython_blas cimport srotg, drotg
from scipy.linalg.cython_blas cimport srot, drot
from scipy.linalg.cython_blas cimport sgemv, dgemv
from scipy.linalg.cython_blas cimport sger, dger
from scipy.linalg.cython_blas cimport sgemm, dgemm


################
# BLAS Level 1 #
################

cdef floating _dot(int n, floating *x, int incx,
                   floating *y, int incy) nogil:
    """x.T.y"""
    if floating is float:
        return sdot(&n, x, &incx, y, &incy)
    else:
        return ddot(&n, x, &incx, y, &incy)


cpdef _dot_memview(floating[::1] x, floating[::1] y):
    return _dot(x.shape[0], &x[0], 1, &y[0], 1)


cdef floating _asum(int n, floating *x, int incx) nogil:
    """sum(|x_i|)"""
    if floating is float:
        return sasum(&n, x, &incx)
    else:
        return dasum(&n, x, &incx)


cpdef _asum_memview(floating[::1] x):
    return _asum(x.shape[0], &x[0], 1)


cdef void _axpy(int n, floating alpha, floating *x, int incx,
                floating *y, int incy) nogil:
    """y := alpha * x + y"""
    if floating is float:
        saxpy(&n, &alpha, x, &incx, y, &incy)
    else:
        daxpy(&n, &alpha, x, &incx, y, &incy)


cpdef _axpy_memview(floating alpha, floating[::1] x, floating[::1] y):
    _axpy(x.shape[0], alpha, &x[0], 1, &y[0], 1)


cdef floating _nrm2(int n, floating *x, int incx) nogil:
    """sqrt(sum((x_i)^2))"""
    if floating is float:
        return snrm2(&n, x, &incx)
    else:
        return dnrm2(&n, x, &incx)


cpdef _nrm2_memview(floating[::1] x):
    return _nrm2(x.shape[0], &x[0], 1)


cdef void _copy(int n, floating *x, int incx, floating *y, int incy) nogil:
    """y := x"""
    if floating is float:
        scopy(&n, x, &incx, y, &incy)
    else:
        dcopy(&n, x, &incx, y, &incy)


cpdef _copy_memview(floating[::1] x, floating[::1] y):
    _copy(x.shape[0], &x[0], 1, &y[0], 1)


cdef void _scal(int n, floating alpha, floating *x, int incx) nogil:
    """x := alpha * x"""
    if floating is float:
        sscal(&n, &alpha, x, &incx)
    else:
        dscal(&n, &alpha, x, &incx)


cpdef _scal_memview(floating alpha, floating[::1] x):
    _scal(x.shape[0], alpha, &x[0], 1)


cdef void _rotg(floating *a, floating *b, floating *c, floating *s) nogil:
    """Generate plane rotation"""
    if floating is float:
        srotg(a, b, c, s)
    else:
        drotg(a, b, c, s)


cpdef _rotg_memview(floating a, floating b, floating c, floating s):
    _rotg(&a, &b, &c, &s)
    return a, b, c, s


cdef void _rot(int n, floating *x, int incx, floating *y, int incy,
               floating c, floating s) nogil:
    """Apply plane rotation"""
    if floating is float:
        srot(&n, x, &incx, y, &incy, &c, &s)
    else:
        drot(&n, x, &incx, y, &incy, &c, &s)


cpdef _rot_memview(floating[::1] x, floating[::1] y, floating c, floating s):
    _rot(x.shape[0], &x[0], 1, &y[0], 1, c, s)


################
# BLAS Level 2 #
################

cdef void _gemv(BLAS_Order order, BLAS_Trans ta, int m, int n, floating alpha,
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
    else:
        if floating is float:
            sgemv(&ta_, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)
        else:
            dgemv(&ta_, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)


cpdef _gemv_memview(BLAS_Trans ta, floating alpha, floating[:, :] A,
                    floating[::1] x, floating beta, floating[::1] y):
    cdef:
        int m = A.shape[0]
        int n = A.shape[1]
        BLAS_Order order = ColMajor if A.strides[0] == A.itemsize else RowMajor
        int lda = m if order == ColMajor else n

    _gemv(order, ta, m, n, alpha, &A[0, 0], lda, &x[0], 1, beta, &y[0], 1)


cdef void _ger(BLAS_Order order, int m, int n, floating alpha, floating *x,
               int incx, floating *y, int incy, floating *A, int lda) nogil:
    """A := alpha * x.y.T + A"""
    if order == RowMajor:
        if floating is float:
            sger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
        else:
            dger(&n, &m, &alpha, y, &incy, x, &incx, A, &lda)
    else:
        if floating is float:
            sger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)
        else:
            dger(&m, &n, &alpha, x, &incx, y, &incy, A, &lda)


cpdef _ger_memview(floating alpha, floating[::1] x, floating[::] y,
                   floating[:, :] A):
    cdef:
        int m = A.shape[0]
        int n = A.shape[1]
        BLAS_Order order = ColMajor if A.strides[0] == A.itemsize else RowMajor
        int lda = m if order == ColMajor else n

    _ger(order, m, n, alpha, &x[0], 1, &y[0], 1, &A[0, 0], lda)


################
# BLAS Level 3 #
################

cdef void _gemm(BLAS_Order order, BLAS_Trans ta, BLAS_Trans tb, int m, int n,
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
    else:
        if floating is float:
            sgemm(&ta_, &tb_, &m, &n, &k, &alpha, A,
                  &lda, B, &ldb, &beta, C, &ldc)
        else:
            dgemm(&ta_, &tb_, &m, &n, &k, &alpha, A,
                  &lda, B, &ldb, &beta, C, &ldc)


cpdef _gemm_memview(BLAS_Trans ta, BLAS_Trans tb, floating alpha,
                    floating[:, :] A, floating[:, :] B, floating beta,
                    floating[:, :] C):
    cdef:
        int m = A.shape[0] if ta == NoTrans else A.shape[1]
        int n = B.shape[1] if tb == NoTrans else B.shape[0]
        int k = A.shape[1] if ta == NoTrans else A.shape[0]
        int lda, ldb, ldc
        BLAS_Order order = ColMajor if A.strides[0] == A.itemsize else RowMajor

    if order == RowMajor:
        lda = k if ta == NoTrans else m
        ldb = n if tb == NoTrans else k
        ldc = n
    else:
        lda = m if ta == NoTrans else k
        ldb = k if tb == NoTrans else n
        ldc = m

    _gemm(order, ta, tb, m, n, k, alpha, &A[0, 0],
          lda, &B[0, 0], ldb, beta, &C[0, 0], ldc)
