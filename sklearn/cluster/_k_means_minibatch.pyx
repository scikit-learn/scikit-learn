from cython cimport floating
from cython.parallel cimport parallel, prange
from libc.stdlib cimport malloc, free


def _minibatch_update_dense(
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        int n_threads):
    """Update of the centers for dense MiniBatchKMeans.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features), dtype=floating
        The observations to cluster.

    sample_weight : ndarray of shape (n_samples,), dtype=floating
        The weights for each observation in X.

    centers_old : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.

    weight_sums : ndarray of shape (n_clusters,), dtype=floating
        Current sums of the accumulated weights for each center.

    labels : ndarray of shape (n_samples,), dtype=int
        labels assignment.

    n_threads : int
        The number of threads to be used by openmp.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers_old.shape[0]
        int cluster_idx

        int *indices

    with nogil, parallel(num_threads=n_threads):
        indices = <int*> malloc(n_samples * sizeof(int))

        for cluster_idx in prange(n_clusters, schedule="static"):
            update_center_dense(cluster_idx, X, sample_weight,
                                centers_old, centers_new, weight_sums, labels,
                                indices)

        free(indices)


cdef void update_center_dense(
        int cluster_idx,
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        int *indices) noexcept nogil:        # TMP
    """Update of a single center for dense MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        floating alpha
        int n_indices
        int k, sample_idx, feature_idx

        floating wsum = 0

    # indices = np.where(labels == cluster_idx)[0]
    k = 0
    for sample_idx in range(n_samples):
        if labels[sample_idx] == cluster_idx:
            indices[k] = sample_idx
            wsum += sample_weight[sample_idx]
            k += 1
    n_indices = k

    if wsum > 0:
        # Undo the previous count-based scaling for this cluster center
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx] * weight_sums[cluster_idx]

        # Update cluster with new point members
        for k in range(n_indices):
            sample_idx = indices[k]
            for feature_idx in range(n_features):
                centers_new[cluster_idx, feature_idx] += X[sample_idx, feature_idx] * sample_weight[sample_idx]

        # Update the count statistics for this center
        weight_sums[cluster_idx] += wsum

        # Rescale to compute mean of all points (old and new)
        alpha = 1 / weight_sums[cluster_idx]
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] *= alpha
    else:
        # No sample was assigned to this cluster in this batch of data
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx]


def _minibatch_update_sparse(
        X,                                   # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        int n_threads):
    """Update of the centers for sparse MiniBatchKMeans.

    Parameters
    ----------
    X : sparse matrix of shape (n_samples, n_features), dtype=floating
        The observations to cluster. Must be in CSR format.

    sample_weight : ndarray of shape (n_samples,), dtype=floating
        The weights for each observation in X.

    centers_old : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.

    weight_sums : ndarray of shape (n_clusters,), dtype=floating
        Current sums of the accumulated weights for each center.

    labels : ndarray of shape (n_samples,), dtype=int
        labels assignment.

    n_threads : int
        The number of threads to be used by openmp.
    """
    cdef:
        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr
        int n_samples = X.shape[0]
        int n_clusters = centers_old.shape[0]
        int cluster_idx

        int *indices

    with nogil, parallel(num_threads=n_threads):
        indices = <int*> malloc(n_samples * sizeof(int))

        for cluster_idx in prange(n_clusters, schedule="static"):
            update_center_sparse(cluster_idx, X_data, X_indices, X_indptr,
                                 sample_weight, centers_old, centers_new,
                                 weight_sums, labels, indices)

        free(indices)


cdef void update_center_sparse(
        int cluster_idx,
        const floating[::1] X_data,          # IN
        const int[::1] X_indices,            # IN
        const int[::1] X_indptr,             # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        int *indices) noexcept nogil:        # TMP
    """Update of a single center for sparse MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        floating alpha
        int n_indices
        int k, sample_idx, feature_idx

        floating wsum = 0

    # indices = np.where(labels == cluster_idx)[0]
    k = 0
    for sample_idx in range(n_samples):
        if labels[sample_idx] == cluster_idx:
            indices[k] = sample_idx
            wsum += sample_weight[sample_idx]
            k += 1
    n_indices = k

    if wsum > 0:
        # Undo the previous count-based scaling for this cluster center:
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx] * weight_sums[cluster_idx]

        # Update cluster with new point members
        for k in range(n_indices):
            sample_idx = indices[k]
            for feature_idx in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                centers_new[cluster_idx, X_indices[feature_idx]] += X_data[feature_idx] * sample_weight[sample_idx]

        # Update the count statistics for this center
        weight_sums[cluster_idx] += wsum

        # Rescale to compute mean of all points (old and new)
        alpha = 1 / weight_sums[cluster_idx]
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] *= alpha
    else:
        # No sample was assigned to this cluster in this batch of data
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx]
