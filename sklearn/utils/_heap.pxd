# cython: language_level=3
# Heap routines, used in various Cython implementation.

from cython cimport floating

from ._typedefs cimport ITYPE_t

cdef int simultaneous_sort(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size
) nogil except -1

cdef int heap_push(
    floating* values,
    ITYPE_t* indices,
    ITYPE_t size,
    floating val,
    ITYPE_t val_idx,
) nogil except -1
