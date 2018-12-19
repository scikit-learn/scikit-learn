# cython: language_level=3

from cython cimport floating


cpdef enum BLAS_Order:
    RowMajor
    ColMajor


cpdef enum BLAS_Trans:
    Trans = 116    # correspond to 'n'
    NoTrans = 110  # correspond to 't'


# BLAS Level 1 ################################################################
cdef floating _xdot(int, floating*, int, floating*, int) nogil
cdef floating _xasum(int, floating*, int) nogil
cdef void _xaxpy(int, floating, floating*, int, floating*, int) nogil
cdef floating _xnrm2(int, floating*, int) nogil
cdef void _xcopy(int, floating*, int, floating*, int) nogil
cdef void _xscal(int, floating, floating*, int) nogil

# BLAS Level 2 ################################################################
cdef void _xgemv(BLAS_Order, BLAS_Trans, int, int, floating, floating*, int,
                 floating*, int, floating, floating*, int) nogil
cdef void _xger(BLAS_Order, int, int, floating, floating*, int, floating*, int,
                floating*, int) nogil

# BLASLevel 3 ################################################################
cdef void _xgemm(BLAS_Order, BLAS_Trans, BLAS_Trans, int, int, int, floating,
                 floating*, int, floating*, int, floating, floating*,
                 int) nogil
