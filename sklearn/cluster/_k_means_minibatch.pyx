# cython: profile=True, boundscheck=False, wraparound=False, cdivision=True

# TODO: We still need to use ndarrays instead of typed memoryviews when using
# fused types and when the array may be read-only (for instance when it's
# provided by the user). This is fixed in cython > 0.3.

cimport numpy as np
from cython cimport floating
from cython.parallel cimport parallel, prange
from libc.math cimport sqrt
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy


np.import_array()


def _minibatch_update_dense(
        np.ndarray[floating, ndim=2, mode='c'] X,  # IN
        floating[::1] sample_weight,               # IN
        floating[:, ::1] centers_old,              # IN
        floating[:, ::1] centers_new,              # OUT
        floating[::1] weight_sums,                 # INOUT
        int[::1] labels):                          # IN
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
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers_old.shape[0]
        int i

        int *indices
    
    with nogil, parallel():
        indices = <int*> malloc(n_samples * sizeof(int))

        for i in prange(n_clusters):
            update_center_dense(i, &X[0, 0], sample_weight, centers_old,
                                centers_new, weight_sums, labels, indices)
        
        free(indices)


cdef void update_center_dense(
        int i,
        floating *X,                   # IN
        floating[::1] sample_weight,   # IN
        floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,  # OUT
        floating[::1] weight_sums,     # INOUT
        int[::1] labels,               # IN
        int *indices) nogil:           # OUT
    """Update of a single center for dense MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        floating alpha, tmp
        int n_indices
        int j, k, idx

        floating wsum = 0

    # indices = np.where(labels == i)[0]
    k = 0
    for j in range(n_samples):
        if labels[j] == i:
            indices[k] = j
            k += 1
    n_indices = k

    for j in range(n_indices):
        idx = indices[j]
        wsum += sample_weight[idx]

    if wsum > 0:
        # Remove previous count scaling
        for k in range(n_features):
            centers_new[i, k] = centers_old[i, k] * weight_sums[i]

        # Update cluster with new point members
        for j in range(n_indices):
            idx = indices[j]
            for k in range(n_features):
                centers_new[i, k] += X[idx * n_features + k] * sample_weight[idx]

        # Update the count statistics for this center
        weight_sums[i] += wsum

        # Rescale to compute mean of all points (old and new)
        alpha = 1 / weight_sums[i]
        for k in range(n_features):
            centers_new[i, k] *= alpha
    else:
        for k in range(n_features):
            centers_new[i, k] = centers_old[i, k]


def _minibatch_update_sparse(
        X,                             # IN
        floating[::1] sample_weight,   # IN
        floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,  # OUT
        floating[::1] weight_sums,     # INOUT
        int[::1] labels):              # IN
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
    """
    cdef:
        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr
        int n_samples = X.shape[0]
        int n_clusters = centers_old.shape[0]
        int i

        int *indices
    
    with nogil, parallel():
        indices = <int*> malloc(n_samples * sizeof(int))

        for i in prange(n_clusters):
            update_center_sparse(i, X_data, X_indices, X_indptr, sample_weight,
                                 centers_old, centers_new, weight_sums, labels,
                                 indices)
        
        free(indices)


cdef void update_center_sparse(
        int i,
        floating[::1] X_data,          # IN
        int[::1] X_indices,            # IN
        int[::1] X_indptr,             # IN
        floating[::1] sample_weight,   # IN
        floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,  # OUT
        floating[::1] weight_sums,     # INOUT
        int[::1] labels,               # IN
        int *indices) nogil:           # OUT
    """Update of a single center for sparse MinibatchKMeans"""
    cdef:
        int n_samples = sample_weight.shape[0]
        int n_features = centers_old.shape[1]
        floating alpha, tmp
        int n_indices
        int j, k, idx

        floating wsum = 0

    # indices = np.where(labels == i)[0]
    k = 0
    for j in range(n_samples):
        if labels[j] == i:
            indices[k] = j
            k += 1
    n_indices = k

    for j in range(n_indices):
        idx = indices[j]
        wsum += sample_weight[idx]

    if wsum > 0:
        # Remove previous count scaling
        for k in range(n_features):
            centers_new[i, k] = centers_old[i, k] * weight_sums[i]

        # Update cluster with new point members
        for j in range(n_indices):
            idx = indices[j]
            for k in range(X_indptr[idx], X_indptr[idx + 1]):
                centers_new[i, X_indices[k]] += X_data[k] * sample_weight[idx]

        # Update the count statistics for this center
        weight_sums[i] += wsum

        # Rescale to compute mean of all points (old and new)
        alpha = 1 / weight_sums[i]
        for k in range(n_features):
            centers_new[i, k] *= alpha
    else:
        for k in range(n_features):
            centers_new[i, k] = centers_old[i, k]


def _minibatch_update_dense4(np.ndarray[floating, ndim=2, mode='c'] X,
                            floating[::1] sample_weight,
                            floating[:, ::1] centers,
                            floating[:, ::1] centers_new,
                            floating[::1] weight_sums,
                            int[::1] labels):
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int i, j, label
        floating weight_sum, tmp, lr

    # for i in prange(n_samples, nogil=True):
    for i in range(n_samples):
        label = labels[i]

        # update center weight
        weight_sum = weight_sums[label] + sample_weight[i]

        # learning rate
        if weight_sum > 0:
            lr = 1 / weight_sum

            for j in range(n_features):
                centers_new[label, j] = centers[label, j] * (1 - lr) + lr * X[i, j]
        else:
            centers_new[label, j] = centers[label, j]

        weight_sums[label] = weight_sum