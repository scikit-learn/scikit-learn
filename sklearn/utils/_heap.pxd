# cython: language_level=3

from cython cimport floating

from ._typedefs cimport ITYPE_t

cdef int _simultaneous_sort(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size
) nogil except -1

cdef inline int _push(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size,
    floating val,
    ITYPE_t i_val,
) nogil except -1
