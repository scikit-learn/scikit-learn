from ._typedefs cimport DTYPE_t, ITYPE_t

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil

cdef int simultaneous_radix_sort(
    ITYPE_t *values,
    ITYPE_t *indices,
    ITYPE_t size,
    ITYPE_t* value_copies,
    ITYPE_t* index_copies,
) nogil
