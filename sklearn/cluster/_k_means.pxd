from cython cimport floating
cimport numpy as np


cdef void _relocate_empty_clusters_dense(
    np.ndarray[floating, ndim=2, mode='c'],
    floating[::1],
    floating[:, ::1],
    floating[::1],
    int[::1]
)


cdef void _relocate_empty_clusters_sparse(
    floating[::1],
    int[::1],
    int[::1],
    floating[::1],
    floating[:, ::1],
    floating[::1],
    int[::1]
)


cdef void _mean_and_center_shift(
    floating[:, ::1],
    floating[:, ::1],
    floating[::1],
    floating[::1]
)