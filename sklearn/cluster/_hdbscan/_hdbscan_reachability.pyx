# cython: boundscheck=False
# cython: nonecheck=False
# cython: initializedcheck=False
# mutual reachability distance compiutations
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np
cimport numpy as np

from scipy.spatial.distance import pdist, squareform
from scipy.sparse import lil_matrix as sparse_matrix
from sklearn.neighbors import KDTree, BallTree
import gc


def mutual_reachability(distance_matrix, min_points=5, alpha=1.0):
    """Compute the weighted adjacency matrix of the mutual reachability
    graph of a distance matrix.

    Parameters
    ----------
    distance_matrix : ndarray, shape (n_samples, n_samples)
        Array of distances between samples.

    min_points : int, optional (default=5)
        The number of points in a neighbourhood for a point to be considered
        a core point.

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
    size = distance_matrix.shape[0]
    min_points = min(size - 1, min_points)
    try:
        core_distances = np.partition(distance_matrix,
                                      min_points,
                                      axis=0)[min_points]
    except AttributeError:
        core_distances = np.sort(distance_matrix,
                                 axis=0)[min_points]

    if alpha != 1.0:
        distance_matrix = distance_matrix / alpha

    stage1 = np.where(core_distances > distance_matrix,
                      core_distances, distance_matrix)
    result = np.where(core_distances > stage1.T,
                      core_distances.T, stage1.T).T
    return result


cpdef sparse_mutual_reachability(object lil_matrix, np.intp_t min_points=5,
                                 float alpha=1.0, float max_dist=0.):

    cdef np.intp_t i
    cdef np.intp_t j
    cdef np.intp_t n
    cdef np.double_t mr_dist
    cdef list sorted_row_data
    cdef np.ndarray[dtype=np.double_t, ndim=1] core_distance
    cdef np.ndarray[dtype=np.int32_t, ndim=1] nz_row_data
    cdef np.ndarray[dtype=np.int32_t, ndim=1] nz_col_data

    result = sparse_matrix(lil_matrix.shape)
    core_distance = np.empty(lil_matrix.shape[0], dtype=np.double)

    for i in range(lil_matrix.shape[0]):
        sorted_row_data = sorted(lil_matrix.data[i])
        if min_points - 1 < len(sorted_row_data):
            core_distance[i] = sorted_row_data[min_points - 1]
        else:
            core_distance[i] = np.infty

    if alpha != 1.0:
        lil_matrix = lil_matrix / alpha

    nz_row_data, nz_col_data = lil_matrix.nonzero()

    for n in range(nz_row_data.shape[0]):
        i = nz_row_data[n]
        j = nz_col_data[n]

        mr_dist = max(core_distance[i], core_distance[j], lil_matrix[i, j])
        if np.isfinite(mr_dist):
            result[i, j] = mr_dist
        elif max_dist > 0:
            result[i, j] = max_dist

    return result.tocsr()


def kdtree_mutual_reachability(X, distance_matrix, metric, p=2, min_points=5,
                               alpha=1.0, **kwargs):
    dim = distance_matrix.shape[0]
    min_points = min(dim - 1, min_points)

    if metric == 'minkowski':
        tree = KDTree(X, metric=metric, p=p)
    else:
        tree = KDTree(X, metric=metric, **kwargs)

    core_distances = tree.query(X, k=min_points)[0][:, -1]

    if alpha != 1.0:
        distance_matrix = distance_matrix / alpha

    stage1 = np.where(core_distances > distance_matrix,
                      core_distances, distance_matrix)
    result = np.where(core_distances > stage1.T,
                      core_distances.T, stage1.T).T
    return result


def balltree_mutual_reachability(X, distance_matrix, metric, p=2, min_points=5,
                                 alpha=1.0, **kwargs):
    dim = distance_matrix.shape[0]
    min_points = min(dim - 1, min_points)

    tree = BallTree(X, metric=metric, **kwargs)

    core_distances = tree.query(X, k=min_points)[0][:, -1]

    if alpha != 1.0:
        distance_matrix = distance_matrix / alpha

    stage1 = np.where(core_distances > distance_matrix,
                      core_distances, distance_matrix)
    result = np.where(core_distances > stage1.T,
                      core_distances.T, stage1.T).T
    return result


cdef np.ndarray[np.double_t, ndim=1] mutual_reachability_from_pdist(
        np.ndarray[np.double_t, ndim=1] core_distances,
        np.ndarray[np.double_t, ndim=1] dists, np.intp_t dim):

    cdef np.intp_t i
    cdef np.intp_t j
    cdef np.intp_t result_pos

    result_pos = 0
    for i in range(dim):
        for j in range(i + 1, dim):
            if core_distances[i] > core_distances[j]:
                if core_distances[i] > dists[result_pos]:
                    dists[result_pos] = core_distances[i]

            else:
                if core_distances[j] > dists[result_pos]:
                    dists[result_pos] = core_distances[j]

            result_pos += 1

    return dists


def kdtree_pdist_mutual_reachability(X,  metric, p=2, min_points=5, alpha=1.0,
                                     **kwargs):

    dim = X.shape[0]
    min_points = min(dim - 1, min_points)

    if metric == 'minkowski':
        tree = KDTree(X, metric=metric, p=p)
    else:
        tree = KDTree(X, metric=metric, **kwargs)

    core_distances = tree.query(X, k=min_points)[0][:, -1]

    del tree
    gc.collect()

    dists = pdist(X, metric=metric, p=p, **kwargs)

    if alpha != 1.0:
        dists /= alpha

    dists = mutual_reachability_from_pdist(core_distances, dists, dim)

    return dists


def balltree_pdist_mutual_reachability(X, metric, p=2, min_points=5, alpha=1.0,
                                       **kwargs):

    dim = X.shape[0]
    min_points = min(dim - 1, min_points)

    tree = BallTree(X, metric=metric, **kwargs)

    core_distances = tree.query(X, k=min_points)[0][:, -1]

    del tree
    gc.collect()

    dists = pdist(X, metric=metric, p=p, **kwargs)

    if alpha != 1.0:
        dists /= alpha

    dists = mutual_reachability_from_pdist(core_distances, dists, dim)

    return dists
