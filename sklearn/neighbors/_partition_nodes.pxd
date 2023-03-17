from cython cimport floating
from ..utils._typedefs cimport DTYPE_t, ITYPE_t

cdef int partition_node_indices(
        floating *data,
        ITYPE_t *node_indices,
        ITYPE_t split_dim,
        ITYPE_t split_index,
        ITYPE_t n_features,
        ITYPE_t n_points) except -1
