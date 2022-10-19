# mutual reachability distance computations
# Authors: Leland McInnes <leland.mcinnes@gmail.com>
#          Meekail Zain <zainmeekail@gmail.com>
# License: 3-clause BSD

import numpy as np
from scipy.sparse import issparse

cimport numpy as cnp
from cython.parallel cimport prange
from libc.math cimport isfinite, INFINITY


def mutual_reachability_graph(
    distance_matrix, n_neighbors=5, max_distance=0.0, copy=False
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

    n_neighbors : int, default=5
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

    if issparse(distance_matrix):
        # FIXME: since we convert to a CSR matrix then we do not make the operation
        # in-place.
        return _sparse_mutual_reachability_graph(
            distance_matrix, n_neighbors=n_neighbors, max_distance=max_distance
        ).tocsr()

    return _dense_mutual_reachability_graph(distance_matrix, n_neighbors=n_neighbors)


cdef _dense_mutual_reachability_graph(
    cnp.ndarray[dtype=cnp.float64_t, ndim=2] distance_matrix,
    cnp.intp_t n_neighbors=5
):
    """Dense implementation of mutual reachability graph.

    The computation is done in-place, i.e. the distance matrix is modified
    directly.

    Parameters
    ----------
    distance_matrix : ndarray of shape (n_samples, n_samples)
        Array of distances between samples.

    n_neighbors : int, default=5
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
        cnp.intp_t farther_neighbor_idx = n_neighbors - 1
        cnp.float64_t mutual_reachibility_distance
        cnp.float64_t[:] core_distances

    core_distances = np.partition(
        distance_matrix, farther_neighbor_idx, axis=0
    )[farther_neighbor_idx]

    with nogil:
        for i in range(n_samples):
            for j in prange(n_samples):
                mutual_reachibility_distance = max(
                    core_distances[i],
                    core_distances[j],
                    distance_matrix[i, j],
                )
                distance_matrix[i, j] = mutual_reachibility_distance
    return distance_matrix


# TODO: Rewrite for CSR.
cdef _sparse_mutual_reachability_graph(
    object distance_matrix,
    cnp.intp_t n_neighbors=5,
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

    n_neighbors : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    Returns
    -------
    mututal_reachability_graph : sparse matrix of shape (n_samples, n_samples)
        Weighted adjacency matrix of the mutual reachability graph. This object
        is the same as `distance_matrix` since the operation is done in-place.
    """
    cdef:
        cnp.intp_t i, j, sample_idx, n_samples = distance_matrix.shape[0]
        list row_distances
        cnp.intp_t farther_neighbor_idx = n_neighbors - 1
        cnp.float64_t mutual_reachibility_distance
        cnp.float64_t[:] core_distances
        cnp.int32_t[:] nz_row_data, nz_col_data

    core_distances = np.empty(n_samples, dtype=np.float64)

    for i in range(n_samples):
        row_distances = distance_matrix.data[i]
        if farther_neighbor_idx < len(row_distances):
            core_distances[i] = np.partition(
                row_distances, farther_neighbor_idx
            )[farther_neighbor_idx]
        else:
            core_distances[i] = INFINITY

    nz_row_data, nz_col_data = distance_matrix.nonzero()
    for sample_idx in range(nz_row_data.shape[0]):
        i, j = nz_row_data[sample_idx], nz_col_data[sample_idx]
        mutual_reachibility_distance = max(
            core_distances[i], core_distances[j], distance_matrix[i, j]
        )
        if isfinite(mutual_reachibility_distance):
            distance_matrix[i, j] = mutual_reachibility_distance
        elif max_distance > 0:
            # TODO: it seems that we assume that distance_matrix is initialized
            # with zeros.
            distance_matrix[i, j] = max_distance
    return distance_matrix
