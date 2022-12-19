# mutual reachability distance computations
# Authors: Leland McInnes <leland.mcinnes@gmail.com>
#          Meekail Zain <zainmeekail@gmail.com>
#          Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: 3-clause BSD

cimport cython
cimport numpy as cnp

import numpy as np
from scipy.sparse import issparse
from cython cimport floating, integral
from cython.parallel cimport prange
from libc.math cimport isfinite, INFINITY

cnp.import_array()

def mutual_reachability_graph(
    distance_matrix, min_samples=5, max_distance=0.0
):
    """Compute the weighted adjacency matrix of the mutual reachability graph.

    The mutual reachability distance used to build the graph is defined as::

        max(d_core(x_p), d_core(x_q), d(x_p, x_q))

    and the core distance `d_core` is defined as the distance between a point
    `x_p` and its k-th nearest neighbor.

    Note that all computations are done in-place.

    Parameters
    ----------
    distance_matrix : {ndarray, sparse matrix} of shape (n_samples, n_samples)
        Array of distances between samples. If sparse, the array must be in
        `CSR` format.

    min_samples : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    max_distance : float, default=0.0
        The distance which `np.inf` is replaced with. When the true mutual-
        reachability distance is measured to be infinite, it is instead
        truncated to `max_dist`. Only used when `distance_matrix` is a sparse
        matrix.

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
    further_neighbor_idx = min_samples - 1
    if issparse(distance_matrix):
        if distance_matrix.format != "csr":
            raise ValueError(
                "Only sparse CSR matrices are supported for `distance_matrix`."
            )
        _sparse_mutual_reachability_graph(
            distance_matrix.data,
            distance_matrix.indices,
            distance_matrix.indptr,
            distance_matrix.shape[0],
            further_neighbor_idx=further_neighbor_idx,
            max_distance=max_distance,
        )
    else:
        _dense_mutual_reachability_graph(
            distance_matrix, further_neighbor_idx=further_neighbor_idx
        )
    return distance_matrix


def _dense_mutual_reachability_graph(
    cnp.ndarray[dtype=floating, ndim=2] distance_matrix,
    cnp.intp_t further_neighbor_idx,
):
    """Dense implementation of mutual reachability graph.

    The computation is done in-place, i.e. the distance matrix is modified
    directly.

    Parameters
    ----------
    distance_matrix : ndarray of shape (n_samples, n_samples)
        Array of distances between samples.

    further_neighbor_idx : int
        The index of the furthest neighbor to use to define the core distances.
    """
    cdef:
        cnp.intp_t i, j, n_samples = distance_matrix.shape[0]
        floating mutual_reachibility_distance
        floating[::1] core_distances

    # We assume that the distance matrix is symmetric. We choose to sort every
    # row to have the same implementation than the sparse case that requires
    # CSR matrix.
    core_distances = np.ascontiguousarray(
        np.partition(
            distance_matrix, further_neighbor_idx, axis=1
        )[:, further_neighbor_idx]
    )

    with nogil:
        for i in prange(n_samples):
            for j in range(n_samples):
                mutual_reachibility_distance = max(
                    core_distances[i],
                    core_distances[j],
                    distance_matrix[i, j],
                )
                distance_matrix[i, j] = mutual_reachibility_distance

def _sparse_mutual_reachability_graph(
    cnp.ndarray[floating, ndim=1, mode="c"] data,
    cnp.ndarray[integral, ndim=1, mode="c"] indices,
    cnp.ndarray[integral, ndim=1, mode="c"] indptr,
    cnp.intp_t n_samples,
    cnp.intp_t further_neighbor_idx,
    floating max_distance,
):
    """Sparse implementation of mutual reachability graph.

    The computation is done in-place, i.e. the distance matrix is modified
    directly. This implementation only accepts `CSR` format sparse matrices.

    Parameters
    ----------
    distance_matrix : sparse matrix of shape (n_samples, n_samples)
        Sparse matrix of distances between samples. The sparse format should
        be `CSR`.

    further_neighbor_idx : int
        The index of the furthest neighbor to use to define the core distances.

    max_distance : float
        The distance which `np.inf` is replaced with. When the true mutual-
        reachability distance is measured to be infinite, it is instead
        truncated to `max_dist`. Only used when `distance_matrix` is a sparse
        matrix.
    """
    cdef:
        integral i, col_ind, row_ind
        floating mutual_reachibility_distance
        floating[:] core_distances
        floating[:] row_data

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    core_distances = np.empty(n_samples, dtype=dtype)

    for i in range(n_samples):
        row_data = data[indptr[i]:indptr[i + 1]]
        if further_neighbor_idx < row_data.size:
            core_distances[i] = np.partition(
                row_data, further_neighbor_idx
            )[further_neighbor_idx]
        else:
            core_distances[i] = INFINITY

    with nogil:
        for row_ind in range(n_samples):
            for i in range(indptr[row_ind], indptr[row_ind + 1]):
                col_ind = indices[i]
                mutual_reachibility_distance = max(
                    core_distances[row_ind], core_distances[col_ind], data[i]
                )
                if isfinite(mutual_reachibility_distance):
                    data[i] = mutual_reachibility_distance
                elif max_distance > 0:
                    data[i] = max_distance
