# cython: language_level=3

from cython cimport floating


# BLAS Level 1 ################################################################
cdef floating _xdot(int, floating*, int, floating*, int) nogil
cdef floating _xasum(int, floating*, int) nogil
cdef void _xaxpy(int, floating, floating*, int, floating*, int) nogil
cdef floating _xnrm2(int, floating*, int) nogil
cdef void _xcopy(int, floating*, int, floating*, int) nogil
cdef void _xscal(int, floating, floating*, int) nogil

# BLAS Level 2 ################################################################
cdef void _xgemv(char, char, int, int, floating, floating*, int, floating*,
                 int, floating, floating*, int) nogil
cdef void _xger(char, int, int, floating, floating*, int, floating*, int,
                floating*, int) nogil

# BLASLevel 3 ################################################################
cdef void _xgemm(char, char, char, int, int, int, floating, floating*, int,
                 floating*, int, floating, floating*, int) nogil
