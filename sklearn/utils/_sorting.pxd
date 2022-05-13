from cython cimport floating

from ._typedefs cimport DTYPE_t, ITYPE_t


cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil
