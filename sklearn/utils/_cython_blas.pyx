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
include "sklearn/utils/_blas_int.pxi"


################
# BLAS Level 1 #
################

cdef floating _dot(int n, const floating *x, int incx,
                   const floating *y, int incy) noexcept nogil:
    """x.T.y"""
    cdef blas_int n_ = n, incx_ = incx, incy_ = incy
    if floating is float:
        return sdot(&n_, <float *> x, &incx_, <float *> y, &incy_)
    else:
        return ddot(&n_, <double *> x, &incx_, <double *> y, &incy_)


cpdef _dot_memview(const floating[::1] x, const floating[::1] y):
    return _dot(x.shape[0], &x[0], 1, &y[0], 1)


cdef floating _asum(int n, const floating *x, int incx) noexcept nogil:
    """sum(|x_i|)"""
    cdef blas_int n_ = n, incx_ = incx
    if floating is float:
        return sasum(&n_, <float *> x, &incx_)
    else:
        return dasum(&n_, <double *> x, &incx_)


cpdef _asum_memview(const floating[::1] x):
    return _asum(x.shape[0], &x[0], 1)


cdef void _axpy(int n, floating alpha, const floating *x, int incx,
                floating *y, int incy) noexcept nogil:
    """y := alpha * x + y"""
    cdef blas_int n_ = n, incx_ = incx, incy_ = incy
    if floating is float:
        saxpy(&n_, &alpha, <float *> x, &incx_, y, &incy_)
    else:
        daxpy(&n_, &alpha, <double *> x, &incx_, y, &incy_)


cpdef _axpy_memview(floating alpha, const floating[::1] x, floating[::1] y):
    _axpy(x.shape[0], alpha, &x[0], 1, &y[0], 1)


cdef floating _nrm2(int n, const floating *x, int incx) noexcept nogil:
    """sqrt(sum((x_i)^2))"""
    cdef blas_int n_ = n, incx_ = incx
    if floating is float:
        return snrm2(&n_, <float *> x, &incx_)
    else:
        return dnrm2(&n_, <double *> x, &incx_)


cpdef _nrm2_memview(const floating[::1] x):
    return _nrm2(x.shape[0], &x[0], 1)


cdef void _copy(int n, const floating *x, int incx, const floating *y, int incy) noexcept nogil:
    """y := x"""
    cdef blas_int n_ = n, incx_ = incx, incy_ = incy
    if floating is float:
        scopy(&n_, <float *> x, &incx_, <float *> y, &incy_)
    else:
        dcopy(&n_, <double *> x, &incx_, <double *> y, &incy_)


cpdef _copy_memview(const floating[::1] x, const floating[::1] y):
    _copy(x.shape[0], &x[0], 1, &y[0], 1)


cdef void _scal(int n, floating alpha, const floating *x, int incx) noexcept nogil:
    """x := alpha * x"""
    cdef blas_int n_ = n, incx_ = incx
    if floating is float:
        sscal(&n_, &alpha, <float *> x, &incx_)
    else:
        dscal(&n_, &alpha, <double *> x, &incx_)


cpdef _scal_memview(floating alpha, const floating[::1] x):
    _scal(x.shape[0], alpha, &x[0], 1)


cdef void _rotg(floating *a, floating *b, floating *c, floating *s) noexcept nogil:
    """Generate plane rotation"""
    if floating is float:
        srotg(a, b, c, s)
    else:
        drotg(a, b, c, s)


cpdef _rotg_memview(floating a, floating b, floating c, floating s):
    _rotg(&a, &b, &c, &s)
    return a, b, c, s


cdef void _rot(int n, floating *x, int incx, floating *y, int incy,
               floating c, floating s) noexcept nogil:
    """Apply plane rotation"""
    cdef blas_int n_ = n, incx_ = incx, incy_ = incy
    if floating is float:
        srot(&n_, x, &incx_, y, &incy_, &c, &s)
    else:
        drot(&n_, x, &incx_, y, &incy_, &c, &s)


cpdef _rot_memview(floating[::1] x, floating[::1] y, floating c, floating s):
    _rot(x.shape[0], &x[0], 1, &y[0], 1, c, s)


################
# BLAS Level 2 #
################

cdef void _gemv(BLAS_Order order, BLAS_Trans ta, int m, int n, floating alpha,
                const floating *A, int lda, const floating *x, int incx,
                floating beta, floating *y, int incy) noexcept nogil:
    """y := alpha * op(A).x + beta * y"""
    cdef:
        char ta_ = ta
        blas_int m_ = m, n_ = n, lda_ = lda, incx_ = incx, incy_ = incy
    if order == BLAS_Order.RowMajor:
        ta_ = BLAS_Trans.NoTrans if ta == BLAS_Trans.Trans else BLAS_Trans.Trans
        if floating is float:
            sgemv(&ta_, &n_, &m_, &alpha, <float *> A, &lda_, <float *> x,
                  &incx_, &beta, y, &incy_)
        else:
            dgemv(&ta_, &n_, &m_, &alpha, <double *> A, &lda_, <double *> x,
                  &incx_, &beta, y, &incy_)
    else:
        if floating is float:
            sgemv(&ta_, &m_, &n_, &alpha, <float *> A, &lda_, <float *> x,
                  &incx_, &beta, y, &incy_)
        else:
            dgemv(&ta_, &m_, &n_, &alpha, <double *> A, &lda_, <double *> x,
                  &incx_, &beta, y, &incy_)


cpdef _gemv_memview(BLAS_Trans ta, floating alpha, const floating[:, :] A,
                    const floating[::1] x, floating beta, floating[::1] y):
    cdef:
        int m = A.shape[0]
        int n = A.shape[1]
        BLAS_Order order = (
            BLAS_Order.ColMajor if A.strides[0] == A.itemsize else BLAS_Order.RowMajor
        )
        int lda = m if order == BLAS_Order.ColMajor else n

    _gemv(order, ta, m, n, alpha, &A[0, 0], lda, &x[0], 1, beta, &y[0], 1)


cdef void _ger(BLAS_Order order, int m, int n, floating alpha,
               const floating *x, int incx, const floating *y,
               int incy, floating *A, int lda) noexcept nogil:
    """A := alpha * x.y.T + A"""
    cdef blas_int m_ = m, n_ = n, incx_ = incx, incy_ = incy, lda_ = lda
    if order == BLAS_Order.RowMajor:
        if floating is float:
            sger(&n_, &m_, &alpha, <float *> y, &incy_, <float *> x, &incx_, A, &lda_)
        else:
            dger(&n_, &m_, &alpha, <double *> y, &incy_, <double *> x, &incx_, A, &lda_)
    else:
        if floating is float:
            sger(&m_, &n_, &alpha, <float *> x, &incx_, <float *> y, &incy_, A, &lda_)
        else:
            dger(&m_, &n_, &alpha, <double *> x, &incx_, <double *> y, &incy_, A, &lda_)


cpdef _ger_memview(floating alpha, const floating[::1] x,
                   const floating[::1] y, floating[:, :] A):
    cdef:
        int m = A.shape[0]
        int n = A.shape[1]
        BLAS_Order order = (
            BLAS_Order.ColMajor if A.strides[0] == A.itemsize else BLAS_Order.RowMajor
        )
        int lda = m if order == BLAS_Order.ColMajor else n

    _ger(order, m, n, alpha, &x[0], 1, &y[0], 1, &A[0, 0], lda)


################
# BLAS Level 3 #
################

cdef void _gemm(BLAS_Order order, BLAS_Trans ta, BLAS_Trans tb, int m, int n,
                int k, floating alpha, const floating *A, int lda, const floating *B,
                int ldb, floating beta, floating *C, int ldc) noexcept nogil:
    """C := alpha * op(A).op(B) + beta * C"""
    # TODO: Remove the pointer casts below once SciPy uses const-qualification.
    # See: https://github.com/scipy/scipy/issues/14262
    cdef:
        char ta_ = ta
        char tb_ = tb
        blas_int m_ = m, n_ = n, k_ = k, lda_ = lda, ldb_ = ldb, ldc_ = ldc
    if order == BLAS_Order.RowMajor:
        if floating is float:
            sgemm(&tb_, &ta_, &n_, &m_, &k_, &alpha, <float*>B,
                  &ldb_, <float*>A, &lda_, &beta, C, &ldc_)
        else:
            dgemm(&tb_, &ta_, &n_, &m_, &k_, &alpha, <double*>B,
                  &ldb_, <double*>A, &lda_, &beta, C, &ldc_)
    else:
        if floating is float:
            sgemm(&ta_, &tb_, &m_, &n_, &k_, &alpha, <float*>A,
                  &lda_, <float*>B, &ldb_, &beta, C, &ldc_)
        else:
            dgemm(&ta_, &tb_, &m_, &n_, &k_, &alpha, <double*>A,
                  &lda_, <double*>B, &ldb_, &beta, C, &ldc_)


cpdef _gemm_memview(BLAS_Trans ta, BLAS_Trans tb, floating alpha,
                    const floating[:, :] A, const floating[:, :] B, floating beta,
                    floating[:, :] C):
    cdef:
        int m = A.shape[0] if ta == BLAS_Trans.NoTrans else A.shape[1]
        int n = B.shape[1] if tb == BLAS_Trans.NoTrans else B.shape[0]
        int k = A.shape[1] if ta == BLAS_Trans.NoTrans else A.shape[0]
        int lda, ldb, ldc
        BLAS_Order order = (
            BLAS_Order.ColMajor if A.strides[0] == A.itemsize else BLAS_Order.RowMajor
        )

    if order == BLAS_Order.RowMajor:
        lda = k if ta == BLAS_Trans.NoTrans else m
        ldb = n if tb == BLAS_Trans.NoTrans else k
        ldc = n
    else:
        lda = m if ta == BLAS_Trans.NoTrans else k
        ldb = k if tb == BLAS_Trans.NoTrans else n
        ldc = m

    _gemm(order, ta, tb, m, n, k, alpha, &A[0, 0],
          lda, &B[0, 0], ldb, beta, &C[0, 0], ldc)
