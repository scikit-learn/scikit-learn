from ._typedefs cimport DTYPE_t, ITYPE_t

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil
