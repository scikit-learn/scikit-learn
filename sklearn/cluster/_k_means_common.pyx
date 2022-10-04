# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Lars Buitinck
#
# License: BSD 3 clause

# TODO: We still need to use ndarrays instead of typed memoryviews when using
# fused types and when the array may be read-only (for instance when it's
# provided by the user). This is fixed in cython > 0.3.

import numpy as np
from cython cimport floating
from cython.parallel cimport prange
from libc.math cimport sqrt

from ..utils.extmath import row_norms


# Number of samples per data chunk defined as a global constant.
CHUNK_SIZE = 256


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
        int[::1] a_indptr,
        int sample_index,
        int offset,
        floating[:, ::1] b_arr,       # IN
        floating[::1] b_squared_norms,
        int b_index,
        bint squared) nogil:
    """Euclidean distance between a sparse and b dense.
    
    Parameters
    ----------
    a_data : 1D memview of floating
        The sparse data.
    a_indices : 1D memview of int
        Indices of the sparse data.
    a_indptr : 1D memview of int
        Index pointers for sparse data.
    sample_index : int
        The sample index (i.e. the index within the indptr of 'a_data').
    offset : int
        The offset on the 'a_indices'.
    b_arr : 2D memview of floating of shape (n_samples, n_cols)
        The dense data.
    b_squared_norms : 1D memview of floating
        The squared l2 norms of each row in 'b_arr'.
    b_index : int
        The sample index in the dense data.
    squared : bint
        Whether or not to return squared distance, or not.

    Returns
    -------
    result : floating
        The distance between 'a' and 'b'.
    """
    cdef:
        int nnz = a_indices[a_indptr[sample_index] - offset:a_indptr[sample_index + 1] - offset].shape[0]
        int i
        floating tmp, bi
        floating result = 0.0

    for i in range(nnz):
        bi = b_arr[b_index, a_indices[a_indptr[sample_index] - offset + i]]
        tmp = a_data[a_indptr[sample_index] - offset + i] - bi
        result += tmp * tmp - bi * bi

    result += b_squared_norms[b_index]

    if result < 0: result = 0.0

    return result if squared else sqrt(result)


def _euclidean_sparse_dense_wrapper(
        floating[::1] a_data,
        int[::1] a_indices,
        int[::1] a_indptr,
        int sample_index,
        int offset,
        floating[:, ::1] b_arr,
        floating[::1] b_squared_norms,
        int b_index,
        bint squared):
    """Wrapper of _euclidean_sparse_dense for testing purpose"""
    return _euclidean_sparse_dense(
        a_data, a_indices, a_indptr, sample_index, offset,
        b_arr, b_squared_norms, b_index, squared)


cpdef floating _inertia_dense(
        floating[:, ::1] X,           # IN READ-ONLY
        floating[::1] sample_weight,  # IN READ-ONLY
        floating[:, ::1] centers,     # IN
        int[::1] labels,              # IN
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
        X,                            # IN
        floating[::1] sample_weight,  # IN
        floating[:, ::1] centers,     # IN
        int[::1] labels,              # IN
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
                X_data, X_indices, X_indptr, i, 0,
                centers, centers_squared_norms, j, True)
            inertia += sq_dist * sample_weight[i]

    return inertia


cpdef void _relocate_empty_clusters_dense(
        floating[:, ::1] X,                # IN READ-ONLY
        floating[::1] sample_weight,       # IN READ-ONLY
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
        int i, j, k

        floating[::1] distances = np.zeros(n_samples, dtype=X_data.base.dtype)
        floating[::1] centers_squared_norms = row_norms(centers_old, squared=True)

    for i in range(n_samples):
        j = labels[i]
        distances[i] = _euclidean_sparse_dense(
                X_data, X_indices, X_indptr, i, 0,
                centers_old, centers_squared_norms, j, True)

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


def _is_same_clustering(int[::1] labels1, int[::1] labels2, n_clusters):
    """Check if two arrays of labels are the same up to a permutation of the labels"""
    cdef int[::1] mapping = np.full(fill_value=-1, shape=(n_clusters,), dtype=np.int32)
    cdef int i

    for i in range(labels1.shape[0]):
        if mapping[labels1[i]] == -1:
            mapping[labels1[i]] = labels2[i]
        elif mapping[labels1[i]] != labels2[i]:
            return False
    return True
