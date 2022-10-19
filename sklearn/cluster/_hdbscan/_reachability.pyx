# mutual reachability distance computations
# Authors: Leland McInnes <leland.mcinnes@gmail.com>
#          Meekail Zain <zainmeekail@gmail.com>
# License: 3-clause BSD

import numpy as np
from scipy.sparse import issparse

cimport cython
cimport numpy as cnp
from cython.parallel cimport prange
from libc.math cimport isfinite, INFINITY


ctypedef fused integral:
    int
    long long


def mutual_reachability_graph(
    distance_matrix, min_samples=5, max_distance=0.0, copy=False
):
    """Compute the weighted adjacency matrix of the mutual reachability graph.

    The mutual reachability distance used to build the graph is defined as::

        max(d_core(x_p), d_core(x_q), d(x_p, x_q))

    and the core distance `d_core` is defined as the distance between a point
    `x_p` and its k-th nearest neighbor.

    Parameters
    ----------
    distance_matrix : {ndarray, sparse matrix} of shape (n_samples, n_samples)
        Array of distances between samples. If sparse, the array must be in
        `LIL` format.

    min_samples : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    max_distance : float, default=0.0
        The distance which `np.inf` is replaced with. When the true mutual-
        reachability distance is measured to be infinite, it is instead
        truncated to `max_dist`. Only used when `distance_matrix` is a sparse
        matrix.

    copy : bool, default=False
        Whether or not to compute the mutual reachinbility graph in-place, i.e.
        modifying directly `distance_matrix`.

    Returns
    -------
    mututal_reachability_graph: {ndarray, sparse matrix} of shape \
            (n_samples, n_samples)
        Weighted adjacency matrix of the mutual reachability graph.

    References
    ----------
    .. [1] Campello, R. J., Moulavi, D., & Sander, J. (2013, April).
       Density-based clustering based on hierarchical density estimates.
       In Pacific-Asia Conference on Knowledge Discovery and Data Mining
       (pp. 160-172). Springer Berlin Heidelberg.
    """
    if copy:
        distance_matrix = distance_matrix.copy()

    further_neighbor_idx = min_samples - 1
    if issparse(distance_matrix):
        # FIXME: since we convert to a CSR matrix then we do not make the operation
        # in-place.
        distance_matrix = distance_matrix.tocsc()
        _sparse_mutual_reachability_graph(
            distance_matrix.data,
            distance_matrix.indices,
            distance_matrix.indptr,
            distance_matrix.shape,
            further_neighbor_idx=further_neighbor_idx,
            max_distance=max_distance,
        )
    else:
        _dense_mutual_reachability_graph(
            distance_matrix, further_neighbor_idx=further_neighbor_idx
        )
    return distance_matrix



cdef _dense_mutual_reachability_graph(
    cnp.ndarray[dtype=cnp.float64_t, ndim=2] distance_matrix,
    cnp.intp_t further_neighbor_idx=5
):
    """Dense implementation of mutual reachability graph.

    The computation is done in-place, i.e. the distance matrix is modified
    directly.

    Parameters
    ----------
    distance_matrix : ndarray of shape (n_samples, n_samples)
        Array of distances between samples.

    min_samples : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    Returns
    -------
    mututal_reachability_graph : ndarray of shape (n_samples, n_samples)
        Weighted adjacency matrix of the mutual reachability graph. This object
        is the same as `distance_matrix` since the operation is done in-place.
    """
    cdef:
        cnp.intp_t i, j, n_samples = distance_matrix.shape[0]
        cnp.float64_t mutual_reachibility_distance
        cnp.float64_t[:] core_distances

    core_distances = np.partition(
        distance_matrix, further_neighbor_idx, axis=0
    )[further_neighbor_idx]

    with nogil:
        for i in range(n_samples):
            for j in prange(n_samples):
                mutual_reachibility_distance = max(
                    core_distances[i],
                    core_distances[j],
                    distance_matrix[i, j],
                )
                distance_matrix[i, j] = mutual_reachibility_distance


# TODO: Rewrite for CSR.
cdef _sparse_mutual_reachability_graph(
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] data,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] indices,
    cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] indptr,
    cnp.intp_t n_samples,
    cnp.intp_t further_neighbor_idx=5,
    cnp.float64_t max_distance=0.0,
):
    """Sparse implementation of mutual reachability graph.

    The computation is done in-place, i.e. the distance matrix is modified
    directly. This implementation only accepts `LIL` format sparse matrices.

    Parameters
    ----------
    distance_matrix : sparse matrix of shape (n_samples, n_samples)
        Sparse matrix of distances between samples. The sparse format should
        be `LIL`.

    min_samples : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    Returns
    -------
    mututal_reachability_graph : sparse matrix of shape (n_samples, n_samples)
        Weighted adjacency matrix of the mutual reachability graph. This object
        is the same as `distance_matrix` since the operation is done in-place.
    """
    cdef:
        cnp.intp_t i, col_ind, row_ind
        cnp.float64_t mutual_reachibility_distance
        cnp.float64_t[:] core_distances
        cnp.float64_t[:] col_data
        cnp.int32_t[:] row_indices

    core_distances = np.empty(n_samples, dtype=np.float64)

    for i in range(n_samples):
        col_data = data[indptr[i]:indptr[i + 1]]
        if further_neighbor_idx < col_data.size:
            core_distances[i] = np.partition(
                col_data, further_neighbor_idx
            )[further_neighbor_idx]
        else:
            core_distances[i] = INFINITY

    for col_ind in range(n_samples):
        for i in range(indptr[col_ind], indptr[col_ind + 1]):
            row_ind = indices[i]
            mutual_reachibility_distance = max(
                core_distances[col_ind], core_distances[row_ind], data[i]
            )
            if isfinite(mutual_reachibility_distance):
                data[i] = mutual_reachibility_distance
            elif max_distance > 0:
                data[i] = max_distance
