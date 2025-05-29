from cython cimport floating
from cython.parallel cimport parallel, prange
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt


def _minibatch_update_dense(
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        bint adaptive_lr,                    # IN
        int n_threads):                      # IN
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

    adaptive_lr : bool, default=False
        Whether to use the adaptive learning rate or not.

    n_threads : int
        The number of threads to be used by openmp.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers_old.shape[0]
        int cluster_idx
        floating wsum_batch=0.0
        int *indices

    if adaptive_lr:
        for sample_idx in range(n_samples):
            wsum_batch+= sample_weight[sample_idx]
    with nogil, parallel(num_threads=n_threads):
        indices = <int*> malloc(n_samples * sizeof(int))
        for cluster_idx in prange(n_clusters, schedule="static"):
            update_center_dense(cluster_idx, X, sample_weight,
                                centers_old, centers_new, weight_sums, labels,
                                indices, wsum_batch, adaptive_lr)

        free(indices)


cdef void update_center_dense(
        int cluster_idx,
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # OUT
        floating[::1] weight_sums,           # INOUT
        const int[::1] labels,               # IN
        int *indices,                        # TMP
        floating wsum_batch,                 # IN
        bint adaptive_lr) noexcept nogil:    # IN
    """Update of a single center for dense MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        int n_indices = 0
        int k, sample_idx, feature_idx
        floating wsum_cluster = 0
        floating old_weight, new_weight
        floating alpha, old_scaling_factor, new_scaling_factor

    k = 0
    for sample_idx in range(n_samples):
        if labels[sample_idx] == cluster_idx:
            indices[k] = sample_idx
            wsum_cluster += sample_weight[sample_idx]
            k += 1
    n_indices = k

    if wsum_cluster > 0:
        old_weight = weight_sums[cluster_idx]
        new_weight = old_weight + wsum_cluster
        weight_sums[cluster_idx] = new_weight

        # We want to compute the new center with the update formula
        # C_j^{i+1} = C_j^i*(1-alpha) + alpha*CM(B_j^i).

        # where:
        # - C_j^i is the center representing the j-th cluster at the i-th
        #   iteration
        # - B_j^i is the batch of samples assigned to the j-th cluster at
        #   the i-th iteration
        # - CM(B_j^i) is the (weighted) mean of the samples assigned to
        #   cluster j in iteration i
        # - alpha is the learning rate

        # In the non-adaptive case, alpha = wsum_cluster/(wsum_cluster+old_weight)
        # where:
        # - wsum_cluster is the weight of the points assigned to the cluster in the
        #   current batch
        # - old_weight is the weight of all points assigned to cluster j
        #   in previous iterations.
        # This is equivalent to computing a weighted average of everything
        # assigned to cluster j so far.

        # In the adaptive case (see https://arxiv.org/abs/2304.00419),
        # alpha = sqrt(wsum_cluster/wsum_batch) where wsum_batch is the weight of
        # the batch. This is similar to an exponential moving average but with
        # an adaptive decay rate.

        # For the sake of efficiency, we don't compute the update explicitly.
        # Instead, we skip computing the mean of the batch and instead
        # compute the update by scaling the old center, adding the weighted
        # sum of the batch, and then scaling again.

        # Let Sigma(B_j^i) be the weighted sum of the points assigned to
        # cluster j in the current batch.
        # Therefore (Sigma(B_j^i) = wsum_cluster * CM(B_j^i)).

        # We can rewrite the update formula as:
        # C_{i+1} = C^{i}_j*(1-alpha) + (alpha/wsum_cluster)*Sigma(B_j^i)
        #         = (alpha/wsum_cluster)[C^{i}_j*(1-alpha)(wsum_cluster/alpha) + Sigma(B_j^i)]

        # In the adaptive case, nothing simplifies so we just use the formula
        # as is.
        # In the non-adaptive case, things simplify and we have
        # - (1-alpha)*(wsum_cluster/alpha)
        #   = (old_weight/(w_sum+old_weight))*(wsum_cluster+old_weight) = old_weight
        # - (alpha/wsum_cluster) = 1/(wsum_cluster+old_weight)

        if adaptive_lr:
            alpha = sqrt(wsum_cluster/wsum_batch)
            old_scaling_factor = (1.0-alpha) * (wsum_cluster/alpha)
            new_scaling_factor = alpha/wsum_cluster
        else:
            old_scaling_factor = old_weight
            new_scaling_factor = 1.0/new_weight

        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx]* old_scaling_factor
        for k in range(n_indices):
            sample_idx = indices[k]
            for feature_idx in range(n_features):
                centers_new[cluster_idx, feature_idx] += X[sample_idx, feature_idx] * sample_weight[sample_idx]
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] *= new_scaling_factor

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
        bint adaptive_lr,                    # IN
        int n_threads):                      # IN
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

    adaptive_lr : bool (default: False)
        Whether to use the adaptive learning rate or not.

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
        floating wsum_batch=0.0
        int *indices

    if adaptive_lr:
        for sample_idx in range(n_samples):
            wsum_batch+= sample_weight[sample_idx]
    with nogil, parallel(num_threads=n_threads):
        indices = <int*> malloc(n_samples * sizeof(int))

        for cluster_idx in prange(n_clusters, schedule="static"):
            update_center_sparse(cluster_idx, X_data, X_indices, X_indptr,
                                 sample_weight, centers_old, centers_new,
                                 weight_sums, labels, indices, wsum_batch, adaptive_lr)
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
        int *indices,                        # TMP
        floating wsum_batch,                 # IN
        bint adaptive_lr) noexcept nogil:    # IN
    """Update of a single center for sparse MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        int n_indices = 0
        int k, sample_idx, feature_idx, ptr
        floating wsum_cluster = 0.0
        floating old_weight, new_weight
        floating alpha, old_scaling_factor, new_scaling_factor

    k = 0
    for sample_idx in range(n_samples):
        if labels[sample_idx] == cluster_idx:
            indices[k] = sample_idx
            wsum_cluster += sample_weight[sample_idx]
            k += 1
    n_indices = k

    if wsum_cluster > 0:
        # See update_center_dense for details

        old_weight = weight_sums[cluster_idx]
        new_weight = old_weight + wsum_cluster
        weight_sums[cluster_idx] = new_weight

        if adaptive_lr:
            alpha = sqrt(wsum_cluster / wsum_batch)
            old_scaling_factor = (1.0 - alpha) * (wsum_cluster / alpha)
            new_scaling_factor = alpha / wsum_cluster
        else:
            old_scaling_factor = old_weight
            new_scaling_factor = 1.0 / new_weight

        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = (
                centers_old[cluster_idx, feature_idx] * old_scaling_factor
            )
        for i in range(n_indices):
            sample_idx = indices[i]
            for ptr in range(X_indptr[sample_idx], X_indptr[sample_idx+1]):
                feature_idx = X_indices[ptr]
                centers_new[cluster_idx, feature_idx] += X_data[ptr] * sample_weight[sample_idx]
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] *= new_scaling_factor
    else:
        # No sample was assigned to this cluster in this batch of data
        for feature_idx in range(n_features):
            centers_new[cluster_idx, feature_idx] = centers_old[cluster_idx, feature_idx]
