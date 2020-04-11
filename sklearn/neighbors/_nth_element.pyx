cdef extern from "_nth_element_inner.h":
    void partition_node_indices_inner[D, I](
            D *data,
            I *node_indices,
            I split_dim,
            I split_index,
            I n_features,
            I n_points) except +


cdef int partition_node_indices(
        DTYPE_t *data,
        ITYPE_t *node_indices,
        ITYPE_t split_dim,
        ITYPE_t split_index,
        ITYPE_t n_features,
        ITYPE_t n_points) except -1:
    partition_node_indices_inner(
        data,
        node_indices,
        split_dim,
        split_index,
        n_features,
        n_points)
    return 0
