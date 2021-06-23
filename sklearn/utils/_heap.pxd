#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

cimport cython
cimport numpy as np
from libc.math cimport fabs, sqrt, exp, cos, pow
from cython cimport floating, integral

from ._typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t
from ._typedefs import DTYPE, ITYPE

cdef inline void dual_swap(floating* darr, ITYPE_t* iarr,
                           ITYPE_t i1, ITYPE_t i2) nogil:
    """swap the values at inex i1 and i2 of both darr and iarr"""
    cdef floating dtmp = darr[i1]
    darr[i1] = darr[i2]
    darr[i2] = dtmp

    cdef ITYPE_t itmp = iarr[i1]
    iarr[i1] = iarr[i2]
    iarr[i2] = itmp

cdef int _simultaneous_sort(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size
) nogil except -1

cdef int _push(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size,
    floating val,
    ITYPE_t i_val,
) nogil except -1
