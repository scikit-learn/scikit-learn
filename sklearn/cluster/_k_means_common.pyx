# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from cython cimport floating
from cython.parallel cimport prange
from libc.math cimport sqrt

from sklearn.utils.extmath import row_norms


# Number of samples per data chunk defined as a global constant.
CHUNK_SIZE = 256


cdef floating _euclidean_dense_dense(
        const floating* a,  # IN
        const floating* b,  # IN
        int n_features,
        bint squared
) noexcept nogil:
    """Euclidean distance between a dense and b dense"""
    cdef:
        int i
        int n = n_features // 4
        int rem = n_features % 4
        floating result = 0

    # We manually unroll the loop for better cache optimization.
    for i in range(n):
        result += (
            (a[0] - b[0]) * (a[0] - b[0]) +
            (a[1] - b[1]) * (a[1] - b[1]) +
            (a[2] - b[2]) * (a[2] - b[2]) +
            (a[3] - b[3]) * (a[3] - b[3])
        )
        a += 4
        b += 4

    for i in range(rem):
        result += (a[i] - b[i]) * (a[i] - b[i])

    return result if squared else sqrt(result)


def _euclidean_dense_dense_wrapper(
    const floating[::1] a,
    const floating[::1] b,
    bint squared
):
    """Wrapper of _euclidean_dense_dense for testing purpose"""
    return _euclidean_dense_dense(&a[0], &b[0], a.shape[0], squared)


cdef floating _euclidean_sparse_dense(
        const floating[::1] a_data,  # IN
        const int[::1] a_indices,    # IN
        const floating[::1] b,       # IN
        floating b_squared_norm,
        bint squared
) noexcept nogil:
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

    if result < 0:
        result = 0.0

    return result if squared else sqrt(result)


def _euclidean_sparse_dense_wrapper(
        const floating[::1] a_data,
        const int[::1] a_indices,
        const floating[::1] b,
        floating b_squared_norm,
        bint squared
):
    """Wrapper of _euclidean_sparse_dense for testing purpose"""
    return _euclidean_sparse_dense(
        a_data, a_indices, b, b_squared_norm, squared)


cpdef floating _inertia_dense(
        const floating[:, ::1] X,           # IN
        const floating[::1] sample_weight,  # IN
        const floating[:, ::1] centers,     # IN
        const int[::1] labels,              # IN
        int n_threads,
        int single_label=-1,
):
    """Compute inertia for dense input data

    Sum of squared distance between each sample and its assigned center.

    If single_label is >= 0, the inertia is computed only for that label.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

    for i in prange(n_samples, nogil=True, num_threads=n_threads,
                    schedule='static'):
        j = labels[i]
        if single_label < 0 or single_label == j:
            sq_dist = _euclidean_dense_dense(&X[i, 0], &centers[j, 0],
                                             n_features, True)
            inertia += sq_dist * sample_weight[i]

    return inertia


cpdef floating _inertia_sparse(
        X,                                  # IN
        const floating[::1] sample_weight,  # IN
        const floating[:, ::1] centers,     # IN
        const int[::1] labels,              # IN
        int n_threads,
        int single_label=-1,
):
    """Compute inertia for sparse input data

    Sum of squared distance between each sample and its assigned center.

    If single_label is >= 0, the inertia is computed only for that label.
    """
    cdef:
        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        int n_samples = X.shape[0]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

        floating[::1] centers_squared_norms = row_norms(centers, squared=True)

    for i in prange(n_samples, nogil=True, num_threads=n_threads,
                    schedule='static'):
        j = labels[i]
        if single_label < 0 or single_label == j:
            sq_dist = _euclidean_sparse_dense(
                X_data[X_indptr[i]: X_indptr[i + 1]],
                X_indices[X_indptr[i]: X_indptr[i + 1]],
                centers[j], centers_squared_norms[j], True)
            inertia += sq_dist * sample_weight[i]

    return inertia


cpdef void _relocate_empty_clusters_dense(
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # INOUT
        floating[::1] weight_in_clusters,    # INOUT
        const int[::1] labels                # IN
):
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

    if np.max(distances) == 0:
        # Happens when there are more clusters than non-duplicate samples. Relocating
        # is pointless in this case.
        return

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
        const floating[::1] X_data,          # IN
        const int[::1] X_indices,            # IN
        const int[::1] X_indptr,             # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # INOUT
        floating[::1] weight_in_clusters,    # INOUT
        const int[::1] labels                # IN
):
    """Relocate centers which have no sample assigned to them."""
    cdef:
        int[::1] empty_clusters = np.where(np.equal(weight_in_clusters, 0))[0].astype(np.int32)
        int n_empty = empty_clusters.shape[0]

    if n_empty == 0:
        return

    cdef:
        int n_samples = X_indptr.shape[0] - 1
        int i, j, k

        floating[::1] distances = np.zeros(n_samples, dtype=X_data.base.dtype)
        floating[::1] centers_squared_norms = row_norms(centers_old, squared=True)

    for i in range(n_samples):
        j = labels[i]
        distances[i] = _euclidean_sparse_dense(
            X_data[X_indptr[i]: X_indptr[i + 1]],
            X_indices[X_indptr[i]: X_indptr[i + 1]],
            centers_old[j], centers_squared_norms[j], True)

    if np.max(distances) == 0:
        # Happens when there are more clusters than non-duplicate samples. Relocating
        # is pointless in this case.
        return

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
        floating[:, ::1] centers,               # INOUT
        const floating[::1] weight_in_clusters  # IN
):
    """Average new centers wrt weights."""
    cdef:
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]
        int j, k
        floating alpha
        int argmax_weight = np.argmax(weight_in_clusters)

    for j in range(n_clusters):
        if weight_in_clusters[j] > 0:
            alpha = 1.0 / weight_in_clusters[j]
            for k in range(n_features):
                centers[j, k] *= alpha
        else:
            # For convenience, we avoid setting empty clusters at the origin but place
            # them at the location of the biggest cluster.
            for k in range(n_features):
                centers[j, k] = centers[argmax_weight, k]


cdef void _center_shift(
        const floating[:, ::1] centers_old,  # IN
        const floating[:, ::1] centers_new,  # IN
        floating[::1] center_shift           # OUT
):
    """Compute shift between old and new centers."""
    cdef:
        int n_clusters = centers_old.shape[0]
        int n_features = centers_old.shape[1]
        int j

    for j in range(n_clusters):
        center_shift[j] = _euclidean_dense_dense(
            &centers_new[j, 0], &centers_old[j, 0], n_features, False)


def _is_same_clustering(
    const int[::1] labels1,
    const int[::1] labels2,
    n_clusters
):
    """Check if two arrays of labels are the same up to a permutation of the labels"""
    cdef int[::1] mapping = np.full(fill_value=-1, shape=(n_clusters,), dtype=np.int32)
    cdef int i

    for i in range(labels1.shape[0]):
        if mapping[labels1[i]] == -1:
            mapping[labels1[i]] = labels2[i]
        elif mapping[labels1[i]] != labels2[i]:
            return False
    return True
