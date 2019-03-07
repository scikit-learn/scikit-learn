# cython: language_level=3

from cython cimport floating
cimport numpy as np


cpdef enum BLAS_Order:
    RowMajor  # C contiguous
    ColMajor  # Fortran contiguous


cpdef enum BLAS_Trans:
    NoTrans = 110  # correspond to 'n'
    Trans = 116    # correspond to 't'


cpdef enum BLAS_UPLO:
    Upper = 117  # correspond to 'u'
    Lower = 108  # correspond to 'l'


# BLAS Level 1 ################################################################
cdef floating _dot(int, floating*, int, floating*, int) nogil

cdef floating _asum(int, floating*, int) nogil

cdef void _axpy(int, floating, floating*, int, floating*, int) nogil

cdef floating _nrm2(int, floating*, int) nogil

cdef void _copy(int, floating*, int, floating*, int) nogil

cdef void _scal(int, floating, floating*, int) nogil

cdef void _rotg(floating*, floating*, floating*, floating*) nogil

cdef void _rot(int, floating*, int, floating*, int, floating, floating) nogil

# BLAS Level 2 ################################################################
cdef void _gemv(BLAS_Order, BLAS_Trans, int, int, floating, floating*, int,
                floating*, int, floating, floating*, int) nogil

cdef void _ger(BLAS_Order, int, int, floating, floating*, int, floating*, int,
               floating*, int) nogil

# BLASLevel 3 ################################################################
cdef void _gemm(BLAS_Order, BLAS_Trans, BLAS_Trans, int, int, int, floating,
                floating*, int, floating*, int, floating, floating*,
                int) nogil

cdef void _syrk(BLAS_Order, BLAS_UPLO, BLAS_Trans, int, int, floating,
                floating*, int, floating, floating*, int) nogil

cdef _syrk_helper(BLAS_UPLO, BLAS_Trans, floating,
                  np.ndarray[floating, ndim=2], floating, floating[:, :])