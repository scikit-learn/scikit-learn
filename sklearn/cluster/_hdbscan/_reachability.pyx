# mutual reachability distance compiutations
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np
cimport numpy as cnp

import gc

from scipy.sparse import issparse
from scipy.spatial.distance import pdist, squareform

from ...neighbors import BallTree, KDTree

def mutual_reachability(distance_matrix, min_points=5, max_dist=0.0):
    """Compute the weighted adjacency matrix of the mutual reachability
    graph of a distance matrix.

    Parameters
    ----------
    distance_matrix : {ndarray or sparse matrix} of shape (n_samples, n_samples)
        Array of distances between samples.

    min_points : int, default=5
        The number of points in a neighbourhood for a point to be considered
        a core point.

    max_dist : float, default=0.0
        The distance which `np.inf` is replaced with. When the true mutual-
        reachability distance is measured to be infinite, it is instead
        truncated to `max_dist`.

    Returns
    -------
    mututal_reachability: ndarray, shape (n_samples, n_samples)
        Weighted adjacency matrix of the mutual reachability graph.

    References
    ----------
    .. [1] Campello, R. J., Moulavi, D., & Sander, J. (2013, April).
       Density-based clustering based on hierarchical density estimates.
       In Pacific-Asia Conference on Knowledge Discovery and Data Mining
       (pp. 160-172). Springer Berlin Heidelberg.
    """
    if issparse(distance_matrix):
        return _sparse_mutual_reachability(
            distance_matrix,
            min_points=min_points,
            max_dist=max_dist
        )
    return _dense_mutual_reachability(distance_matrix, min_points=min_points)

cdef _dense_mutual_reachability(cnp.ndarray distance_matrix, min_points=5):
    cdef cnp.intp_t i, j, n_samples = distance_matrix.shape[0]

    # Account for index offset
    min_points -= 1

    # Compute the core distances for all samples `x_p` corresponding
    # to the distance of the k-th farthest neighbours (including
    # `x_p`).
    core_distances = np.partition(
        distance_matrix,
        min_points,
        axis=0,
    )[min_points]

    for i in range(n_samples):
        for j in range(n_samples):
            mr_dist = max(core_distances[i], core_distances[j], distance_matrix[i, j])
            distance_matrix[i, j] = mr_dist
    return distance_matrix

# Assumes LIL format.
# TODO: Rewrite for CSR.
cdef _sparse_mutual_reachability(
    object distance_matrix,
    cnp.intp_t min_points=5,
    cnp.float64_t max_dist=0.
):
    cdef cnp.intp_t i, j, n
    cdef cnp.float64_t mr_dist
    cdef cnp.ndarray[dtype=cnp.float64_t, ndim=1] core_distances
    cdef cnp.ndarray[dtype=cnp.int32_t, ndim=1] nz_row_data
    cdef cnp.ndarray[dtype=cnp.int32_t, ndim=1] nz_col_data
    core_distances = np.empty(distance_matrix.shape[0], dtype=np.float64)

    # Account for index offset
    min_points -= 1
    for i in range(distance_matrix.shape[0]):
        if min_points < len(distance_matrix.data[i]):
            core_distances[i] = np.partition(distance_matrix.data[i], min_points)[min_points]
        else:
            core_distances[i] = np.infty

    nz_row_data, nz_col_data = distance_matrix.nonzero()

    for n in range(nz_row_data.shape[0]):
        i = nz_row_data[n]
        j = nz_col_data[n]
        mr_dist = max(core_distances[i], core_distances[j], distance_matrix[i, j])
        if np.isfinite(mr_dist):
            distance_matrix[i, j] = mr_dist
        elif max_dist > 0:
            distance_matrix[i, j] = max_dist

    return distance_matrix.tocsr()
