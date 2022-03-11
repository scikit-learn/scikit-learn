from ._typedefs cimport DTYPE_t, ITYPE_t

from cython cimport floating

cdef int partition_node_indices(
    DTYPE_t *data,
    ITYPE_t *node_indices,
    ITYPE_t split_dim,
    ITYPE_t split_index,
    ITYPE_t n_features,
    ITYPE_t n_points
) nogil except -1


cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil
