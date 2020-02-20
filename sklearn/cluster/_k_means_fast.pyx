# cython: profile=True, boundscheck=False, wraparound=False, cdivision=True
# Profiling is enabled by default as the overhead does not seem to be
# measurable on this specific use case.

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Lars Buitinck
#
# License: BSD 3 clause

# TODO: We still need to use ndarrays instead of typed memoryviews when using
# fused types and when the array may be read-only (for instance when it's
# provided by the user). This is fixed in cython > 0.3.

import numpy as np
cimport numpy as np
cimport cython
from cython cimport floating
from libc.math cimport sqrt

from ..utils.extmath import row_norms


np.import_array()


ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INT


cdef floating _euclidean_dense_dense(
        floating* a,  # IN
        floating* b,  # IN
        int n_features,
        bint squared) nogil:
    """Euclidean distance between a dense and b dense"""
    cdef:
        int i
        int n = n_features // 4
        int rem = n_features % 4
        floating result = 0

    # We manually unroll the loop for better cache optimization.
    for i in range(n):
        result += ((a[0] - b[0]) * (a[0] - b[0])
                  +(a[1] - b[1]) * (a[1] - b[1])
                  +(a[2] - b[2]) * (a[2] - b[2])
                  +(a[3] - b[3]) * (a[3] - b[3]))
        a += 4; b += 4

    for i in range(rem):
        result += (a[i] - b[i]) * (a[i] - b[i])

    return result if squared else sqrt(result)


def _euclidean_dense_dense_wrapper(floating[::1] a, floating[::1] b,
                                   bint squared):
    """Wrapper of _euclidean_dense_dense for testing purpose"""
    return _euclidean_dense_dense(&a[0], &b[0], a.shape[0], squared)


cdef floating _euclidean_sparse_dense(
        floating[::1] a_data,  # IN
        int[::1] a_indices,    # IN
        floating[::1] b,       # IN
        floating b_squared_norm,
        bint squared) nogil:
    """Euclidean distance between a sparse and b dense"""
    cdef:
        int nnz = a_indices.shape[0]
        int i
        floating tmp, bi
        floating result = 0.0

    for i in range(nnz):
        bi = b[a_indices[i]]
        tmp = a_data[i] - bi
        result += tmp * tmp - bi * bi

    result += b_squared_norm

    if result < 0: result = 0.0

    return result if squared else sqrt(result)


def _euclidean_sparse_dense_wrapper(
        floating[::1] a_data,
        int[::1] a_indices,
        floating[::1] b,
        floating b_squared_norm,
        bint squared):
    """Wrapper of _euclidean_sparse_dense for testing purpose"""
    return _euclidean_sparse_dense(
        a_data, a_indices, b, b_squared_norm, squared)


cpdef floating _inertia_dense(
        np.ndarray[floating, ndim=2, mode='c'] X,  # IN
        floating[::1] sample_weight,               # IN
        floating[:, ::1] centers,                  # IN
        int[::1] labels):                          # IN
    """Compute inertia for dense input data

    Sum of squared distance between each sample and its assigned center.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

    for i in range(n_samples):
        j = labels[i]
        sq_dist = _euclidean_dense_dense(&X[i, 0], &centers[j, 0],
                                         n_features, True)
        inertia += sq_dist * sample_weight[i]

    return inertia


cpdef floating _inertia_sparse(
        X,                            # IN
        floating[::1] sample_weight,  # IN
        floating[:, ::1] centers,     # IN
        int[::1] labels):             # IN
    """Compute inertia for sparse input data

    Sum of squared distance between each sample and its assigned center.
    """
    cdef:
        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

        floating[::1] centers_squared_norms = row_norms(centers, squared=True)

    for i in range(n_samples):
        j = labels[i]
        sq_dist = _euclidean_sparse_dense(
            X_data[X_indptr[i]: X_indptr[i + 1]],
            X_indices[X_indptr[i]: X_indptr[i + 1]],
            centers[j], centers_squared_norms[j], True)
        inertia += sq_dist * sample_weight[i]

    return inertia


cpdef void _relocate_empty_clusters_dense(
        np.ndarray[floating, ndim=2, mode='c'] X,  # IN
        floating[::1] sample_weight,               # IN
        floating[:, ::1] centers_old,              # IN
        floating[:, ::1] centers_new,              # INOUT
        floating[::1] weight_in_clusters,          # INOUT
        int[::1] labels):                          # IN
    """Relocate centers which have no sample assigned to them."""
    cdef:
        int[::1] empty_clusters = np.where(np.equal(weight_in_clusters, 0))[0].astype(np.int32)
        int n_empty = empty_clusters.shape[0]

    if n_empty == 0:
        return

    cdef:
        int n_features = X.shape[1]

        floating[::1] distances = ((np.asarray(X) - np.asarray(centers_old)[labels])**2).sum(axis=1)
        int[::1] far_from_centers = np.argpartition(distances, -n_empty)[:-n_empty-1:-1].astype(np.int32)

        int new_cluster_id, old_cluster_id, far_idx, idx, k
        floating weight

    for idx in range(n_empty):

        new_cluster_id = empty_clusters[idx]

        far_idx = far_from_centers[idx]
        weight = sample_weight[far_idx]

        old_cluster_id = labels[far_idx]

        for k in range(n_features):
            centers_new[old_cluster_id, k] -= X[far_idx, k] * weight
            centers_new[new_cluster_id, k] = X[far_idx, k] * weight

        weight_in_clusters[new_cluster_id] = weight
        weight_in_clusters[old_cluster_id] -= weight


cpdef void _relocate_empty_clusters_sparse(
        floating[::1] X_data,              # IN
        int[::1] X_indices,                # IN
        int[::1] X_indptr,                 # IN
        floating[::1] sample_weight,       # IN
        floating[:, ::1] centers_old,      # IN
        floating[:, ::1] centers_new,      # INOUT
        floating[::1] weight_in_clusters,  # INOUT
        int[::1] labels):                  # IN
    """Relocate centers which have no sample assigned to them."""
    cdef:
        int[::1] empty_clusters = np.where(np.equal(weight_in_clusters, 0))[0].astype(np.int32)
        int n_empty = empty_clusters.shape[0]

    if n_empty == 0:
        return

    cdef:
        int n_samples = X_indptr.shape[0] - 1
        int n_features = centers_old.shape[1]
        floating x
        int i, j, k

        floating[::1] distances = np.zeros(n_samples, dtype=X_data.base.dtype)
        floating[::1] centers_squared_norms = row_norms(centers_old, squared=True)

    for i in range(n_samples):
        j = labels[i]
        distances[i] = _euclidean_sparse_dense(
            X_data[X_indptr[i]: X_indptr[i + 1]],
            X_indices[X_indptr[i]: X_indptr[i + 1]],
            centers_old[j], centers_squared_norms[j], True)

    cdef:
        int[::1] far_from_centers = np.argpartition(distances, -n_empty)[:-n_empty-1:-1].astype(np.int32)

        int new_cluster_id, old_cluster_id, far_idx, idx
        floating weight

    for idx in range(n_empty):

        new_cluster_id = empty_clusters[idx]

        far_idx = far_from_centers[idx]
        weight = sample_weight[far_idx]

        old_cluster_id = labels[far_idx]

        for k in range(X_indptr[far_idx], X_indptr[far_idx + 1]):
            centers_new[old_cluster_id, X_indices[k]] -= X_data[k] * weight
            centers_new[new_cluster_id, X_indices[k]] = X_data[k] * weight

        weight_in_clusters[new_cluster_id] = weight
        weight_in_clusters[old_cluster_id] -= weight


cdef void _average_centers(
        floating[:, ::1] centers,           # INOUT
        floating[::1] weight_in_clusters):  # IN
    """Average new centers wrt weights."""
    cdef:
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]
        int j, k
        floating alpha

    for j in range(n_clusters):
        if weight_in_clusters[j] > 0:
            alpha = 1.0 / weight_in_clusters[j]
            for k in range(n_features):
                centers[j, k] *= alpha


cdef void _center_shift(
        floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,  # IN
        floating[::1] center_shift):   # OUT
    """Compute shift between old and new centers."""
    cdef:
        int n_clusters = centers_old.shape[0]
        int n_features = centers_old.shape[1]
        int j

    for j in range(n_clusters):
        center_shift[j] = _euclidean_dense_dense(
            &centers_new[j, 0], &centers_old[j, 0], n_features, False)


def _mini_batch_update_csr(X, np.ndarray[floating, ndim=1] sample_weight,
                           np.ndarray[floating, ndim=1] x_squared_norms,
                           np.ndarray[floating, ndim=2] centers,
                           np.ndarray[floating, ndim=1] weight_sums,
                           np.ndarray[INT, ndim=1] nearest_center,
                           np.ndarray[floating, ndim=1] old_center,
                           int compute_squared_diff):
    """Incremental update of the centers for sparse MiniBatchKMeans.

    Parameters
    ----------

    X : CSR matrix, dtype float
        The complete (pre allocated) training set as a CSR matrix.

    centers : array, shape (n_clusters, n_features)
        The cluster centers

    counts : array, shape (n_clusters,)
         The vector in which we keep track of the numbers of elements in a
         cluster

    Returns
    -------
    inertia : float
        The inertia of the batch prior to centers update, i.e. the sum
        of squared distances to the closest center for each sample. This
        is the objective function being minimized by the k-means algorithm.

    squared_diff : float
        The sum of squared update (squared norm of the centers position
        change). If compute_squared_diff is 0, this computation is skipped and
        0.0 is returned instead.

    Both squared diff and inertia are commonly used to monitor the convergence
    of the algorithm.
    """
    cdef:
        np.ndarray[floating, ndim=1] X_data = X.data
        np.ndarray[int, ndim=1] X_indices = X.indices
        np.ndarray[int, ndim=1] X_indptr = X.indptr
        unsigned int n_samples = X.shape[0]
        unsigned int n_clusters = centers.shape[0]
        unsigned int n_features = centers.shape[1]

        unsigned int sample_idx, center_idx, feature_idx
        unsigned int k
        DOUBLE old_weight_sum, new_weight_sum
        DOUBLE center_diff
        DOUBLE squared_diff = 0.0

    # move centers to the mean of both old and newly assigned samples
    for center_idx in range(n_clusters):
        old_weight_sum = weight_sums[center_idx]
        new_weight_sum = old_weight_sum

        # count the number of samples assigned to this center
        for sample_idx in range(n_samples):
            if nearest_center[sample_idx] == center_idx:
                new_weight_sum += sample_weight[sample_idx]

        if new_weight_sum == old_weight_sum:
            # no new sample: leave this center as it stands
            continue

        # rescale the old center to reflect it previous accumulated weight
        # with regards to the new data that will be incrementally contributed
        if compute_squared_diff:
            old_center[:] = centers[center_idx]
        centers[center_idx] *= old_weight_sum

        # iterate of over samples assigned to this cluster to move the center
        # location by inplace summation
        for sample_idx in range(n_samples):
            if nearest_center[sample_idx] != center_idx:
                continue

            # inplace sum with new samples that are members of this cluster
            # and update of the incremental squared difference update of the
            # center position
            for k in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                centers[center_idx, X_indices[k]] += X_data[k]

        # inplace rescale center with updated count
        if new_weight_sum > old_weight_sum:
            # update the count statistics for this center
            weight_sums[center_idx] = new_weight_sum

            # re-scale the updated center with the total new counts
            centers[center_idx] /= new_weight_sum

            # update the incremental computation of the squared total
            # centers position change
            if compute_squared_diff:
                for feature_idx in range(n_features):
                    squared_diff += (old_center[feature_idx]
                                     - centers[center_idx, feature_idx]) ** 2

    return squared_diff
