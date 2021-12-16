# Heap routines, used in various Cython implementations.

from cython cimport floating

from ._typedefs cimport ITYPE_t

cdef int simultaneous_sort(
    floating* dist,
    ITYPE_t* idx,
    ITYPE_t size
) nogil

cdef int heap_push(
    floating* values,
    ITYPE_t* indices,
    ITYPE_t size,
    floating val,
    ITYPE_t val_idx,
) nogil
