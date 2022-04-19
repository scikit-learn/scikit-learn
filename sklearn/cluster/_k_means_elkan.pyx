# Author: Andreas Mueller
#
# Licence: BSD 3 clause

# TODO: We still need to use ndarrays instead of typed memoryviews when using
# fused types and when the array may be read-only (for instance when it's
# provided by the user). This is fixed in cython > 0.3.

import numpy as np
cimport numpy as np
cimport openmp
cimport cython
from cython cimport floating
from cython.parallel import prange, parallel
from libc.math cimport sqrt
from libc.stdlib cimport calloc, free
from libc.string cimport memset, memcpy

from ..utils.extmath import row_norms
from ._k_means_common import CHUNK_SIZE
from ._k_means_common cimport _relocate_empty_clusters_dense
from ._k_means_common cimport _relocate_empty_clusters_sparse
from ._k_means_common cimport _euclidean_dense_dense
from ._k_means_common cimport _euclidean_sparse_dense
from ._k_means_common cimport _average_centers
from ._k_means_common cimport _center_shift


np.import_array()


def init_bounds_dense(
        floating[:, ::1] X,                      # IN READ-ONLY
        floating[:, ::1] centers,                # IN
        floating[:, ::1] center_half_distances,  # IN
        int[::1] labels,                         # OUT
        floating[::1] upper_bounds,              # OUT
        floating[:, ::1] lower_bounds,           # OUT
        int n_threads):
    """Initialize upper and lower bounds for each sample for dense input data.

    Given X, centers and the pairwise distances divided by 2.0 between the
    centers this calculates the upper bounds and lower bounds for each sample.
    The upper bound for each sample is set to the distance between the sample
    and the closest center.

    The lower bound for each sample is a one-dimensional array of n_clusters.
    For each sample i assume that the previously assigned cluster is c1 and the
    previous closest distance is dist, for a new cluster c2, the
    lower_bound[i][c2] is set to distance between the sample and this new
    cluster, if and only if dist > center_half_distances[c1][c2]. This prevents
    computation of unnecessary distances for each sample to the clusters that
    it is unlikely to be assigned to.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features), dtype=floating
        The input data.

    centers : ndarray of shape (n_clusters, n_features), dtype=floating
        The cluster centers.

    center_half_distances : ndarray of shape (n_clusters, n_clusters), \
            dtype=floating
        The half of the distance between any 2 clusters centers.

    labels : ndarray of shape(n_samples), dtype=int
        The label for each sample. This array is modified in place.

    upper_bounds : ndarray of shape(n_samples,), dtype=floating
        The upper bound on the distance between each sample and its closest
        cluster center. This array is modified in place.

    lower_bounds : ndarray, of shape(n_samples, n_clusters), dtype=floating
        The lower bound on the distance between each sample and each cluster
        center. This array is modified in place.

    n_threads : int
        The number of threads to be used by openmp.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers.shape[0]
        int n_features = X.shape[1]

        floating min_dist, dist
        int best_cluster, i, j

    for i in prange(
        n_samples, num_threads=n_threads, schedule='static', nogil=True
    ):
        best_cluster = 0
        min_dist = _euclidean_dense_dense(&X[i, 0], &centers[0, 0],
                                          n_features, False)
        lower_bounds[i, 0] = min_dist
        for j in range(1, n_clusters):
            if min_dist > center_half_distances[best_cluster, j]:
                dist = _euclidean_dense_dense(&X[i, 0], &centers[j, 0],
                                              n_features, False)
                lower_bounds[i, j] = dist
                if dist < min_dist:
                    min_dist = dist
                    best_cluster = j
        labels[i] = best_cluster
        upper_bounds[i] = min_dist


def init_bounds_sparse(
        X,                                       # IN
        floating[:, ::1] centers,                # IN
        floating[:, ::1] center_half_distances,  # IN
        int[::1] labels,                         # OUT
        floating[::1] upper_bounds,              # OUT
        floating[:, ::1] lower_bounds,           # OUT
        int n_threads):
    """Initialize upper and lower bounds for each sample for sparse input data.

    Given X, centers and the pairwise distances divided by 2.0 between the
    centers this calculates the upper bounds and lower bounds for each sample.
    The upper bound for each sample is set to the distance between the sample
    and the closest center.

    The lower bound for each sample is a one-dimensional array of n_clusters.
    For each sample i assume that the previously assigned cluster is c1 and the
    previous closest distance is dist, for a new cluster c2, the
    lower_bound[i][c2] is set to distance between the sample and this new
    cluster, if and only if dist > center_half_distances[c1][c2]. This prevents
    computation of unnecessary distances for each sample to the clusters that
    it is unlikely to be assigned to.

    Parameters
    ----------
    X : sparse matrix of shape (n_samples, n_features), dtype=floating
        The input data. Must be in CSR format.

    centers : ndarray of shape (n_clusters, n_features), dtype=floating
        The cluster centers.

    center_half_distances : ndarray of shape (n_clusters, n_clusters), \
            dtype=floating
        The half of the distance between any 2 clusters centers.

    labels : ndarray of shape(n_samples), dtype=int
        The label for each sample. This array is modified in place.

    upper_bounds : ndarray of shape(n_samples,), dtype=floating
        The upper bound on the distance between each sample and its closest
        cluster center. This array is modified in place.

    lower_bounds : ndarray of shape(n_samples, n_clusters), dtype=floating
        The lower bound on the distance between each sample and each cluster
        center. This array is modified in place.

    n_threads : int
        The number of threads to be used by openmp.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers.shape[0]
        int n_features = X.shape[1]

        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        floating min_dist, dist
        int best_cluster, i, j

        floating[::1] centers_squared_norms = row_norms(centers, squared=True)

    for i in prange(
        n_samples, num_threads=n_threads, schedule='static', nogil=True
    ):
        best_cluster = 0
        min_dist = _euclidean_sparse_dense(
            X_data[X_indptr[i]: X_indptr[i + 1]],
            X_indices[X_indptr[i]: X_indptr[i + 1]],
            centers[0], centers_squared_norms[0], False)

        lower_bounds[i, 0] = min_dist
        for j in range(1, n_clusters):
            if min_dist > center_half_distances[best_cluster, j]:
                dist = _euclidean_sparse_dense(
                    X_data[X_indptr[i]: X_indptr[i + 1]],
                    X_indices[X_indptr[i]: X_indptr[i + 1]],
                    centers[j], centers_squared_norms[j], False)
                lower_bounds[i, j] = dist
                if dist < min_dist:
                    min_dist = dist
                    best_cluster = j
        labels[i] = best_cluster
        upper_bounds[i] = min_dist


def elkan_iter_chunked_dense(
        floating[:, ::1] X,                      # IN READ-ONLY
        floating[::1] sample_weight,             # IN READ-ONLY
        floating[:, ::1] centers_old,            # IN
        floating[:, ::1] centers_new,            # OUT
        floating[::1] weight_in_clusters,        # OUT
        floating[:, ::1] center_half_distances,  # IN
        floating[::1] distance_next_center,      # IN
        floating[::1] upper_bounds,              # INOUT
        floating[:, ::1] lower_bounds,           # INOUT
        int[::1] labels,                         # INOUT
        floating[::1] center_shift,              # OUT
        int n_threads,
        bint update_centers=True):
    """Single iteration of K-means Elkan algorithm with dense input.

    Update labels and centers (inplace), for one iteration, distributed
    over data chunks.

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

    weight_in_clusters : ndarray of shape (n_clusters,), dtype=floating
        Placeholder for the sums of the weights of every observation assigned
        to each center.

    center_half_distances : ndarray of shape (n_clusters, n_clusters), \
            dtype=floating
        Half pairwise distances between centers.

    distance_next_center : ndarray of shape (n_clusters,), dtype=floating
        Distance between each center its closest center.

    upper_bounds : ndarray of shape (n_samples,), dtype=floating
        Upper bound for the distance between each sample and its center,
        updated inplace.

    lower_bounds : ndarray of shape (n_samples, n_clusters), dtype=floating
        Lower bound for the distance between each sample and each center,
        updated inplace.

    labels : ndarray of shape (n_samples,), dtype=int
        labels assignment.

    center_shift : ndarray of shape (n_clusters,), dtype=floating
        Distance between old and new centers.

    n_threads : int
        The number of threads to be used by openmp.

    update_centers : bool
        - If True, the labels and the new centers will be computed, i.e. runs
          the E-step and the M-step of the algorithm.
        - If False, only the labels will be computed, i.e runs the E-step of
          the algorithm. This is useful especially when calling predict on a
          fitted model.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_clusters = centers_new.shape[0]

        # hard-coded number of samples per chunk. Splitting in chunks is
        # necessary to get parallelism. Chunk size chosen to be same as lloyd's
        int n_samples_chunk = CHUNK_SIZE if n_samples > CHUNK_SIZE else n_samples
        int n_chunks = n_samples // n_samples_chunk
        int n_samples_rem = n_samples % n_samples_chunk
        int chunk_idx, n_samples_chunk_eff
        int start, end

        int i, j, k

        floating *centers_new_chunk
        floating *weight_in_clusters_chunk

        openmp.omp_lock_t lock

    # count remainder chunk in total number of chunks
    n_chunks += n_samples != n_chunks * n_samples_chunk

    # number of threads should not be bigger than number of chunks
    n_threads = min(n_threads, n_chunks)

    if update_centers:
        memset(&centers_new[0, 0], 0, n_clusters * n_features * sizeof(floating))
        memset(&weight_in_clusters[0], 0, n_clusters * sizeof(floating))
        openmp.omp_init_lock(&lock)

    with nogil, parallel(num_threads=n_threads):
        # thread local buffers
        centers_new_chunk = <floating*> calloc(n_clusters * n_features, sizeof(floating))
        weight_in_clusters_chunk = <floating*> calloc(n_clusters, sizeof(floating))

        for chunk_idx in prange(n_chunks, schedule='static'):
            start = chunk_idx * n_samples_chunk
            if chunk_idx == n_chunks - 1 and n_samples_rem > 0:
                end = start + n_samples_rem
            else:
                end = start + n_samples_chunk

            _update_chunk_dense(
                X[start: end],
                sample_weight[start: end],
                centers_old,
                center_half_distances,
                distance_next_center,
                labels[start: end],
                upper_bounds[start: end],
                lower_bounds[start: end],
                centers_new_chunk,
                weight_in_clusters_chunk,
                update_centers)

        # reduction from local buffers. The lock is necessary for that to avoid
        # race conditions.
        if update_centers:
            openmp.omp_set_lock(&lock)
            for j in range(n_clusters):
                weight_in_clusters[j] += weight_in_clusters_chunk[j]
                for k in range(n_features):
                    centers_new[j, k] += centers_new_chunk[j * n_features + k]
            openmp.omp_unset_lock(&lock)

        free(centers_new_chunk)
        free(weight_in_clusters_chunk)

    if update_centers:
        openmp.omp_destroy_lock(&lock)
        _relocate_empty_clusters_dense(X, sample_weight, centers_old,
                                       centers_new, weight_in_clusters, labels)

        _average_centers(centers_new, weight_in_clusters)
        _center_shift(centers_old, centers_new, center_shift)

        # update lower and upper bounds
        for i in range(n_samples):
            upper_bounds[i] += center_shift[labels[i]]

            for j in range(n_clusters):
                lower_bounds[i, j] -= center_shift[j]
                if lower_bounds[i, j] < 0:
                    lower_bounds[i, j] = 0


cdef void _update_chunk_dense(
        floating[:, ::1] X,                      # IN READ-ONLY
        floating[::1] sample_weight,             # IN READ-ONLY
        floating[:, ::1] centers_old,            # IN
        floating[:, ::1] center_half_distances,  # IN
        floating[::1] distance_next_center,      # IN
        int[::1] labels,                         # INOUT
        floating[::1] upper_bounds,              # INOUT
        floating[:, ::1] lower_bounds,           # INOUT
        floating *centers_new,                   # OUT
        floating *weight_in_clusters,            # OUT
        bint update_centers) nogil:
    """K-means combined EM step for one dense data chunk.

    Compute the partial contribution of a single data chunk to the labels and
    centers.
    """
    cdef:
        int n_samples = labels.shape[0]
        int n_clusters = centers_old.shape[0]
        int n_features = centers_old.shape[1]

        floating upper_bound, distance
        int i, j, k, label

    for i in range(n_samples):
        upper_bound = upper_bounds[i]
        bounds_tight = 0
        label = labels[i]

        # Next center is not far away from the currently assigned center.
        # Sample might need to be assigned to another center.
        if not distance_next_center[label] >= upper_bound:

            for j in range(n_clusters):

                # If this holds, then center_index is a good candidate for the
                # sample to be relabelled, and we need to confirm this by
                # recomputing the upper and lower bounds.
                if (j != label
                    and (upper_bound > lower_bounds[i, j])
                    and (upper_bound > center_half_distances[label, j])):

                    # Recompute upper bound by calculating the actual distance
                    # between the sample and its current assigned center.
                    if not bounds_tight:
                        upper_bound = _euclidean_dense_dense(
                            &X[i, 0], &centers_old[label, 0], n_features, False)
                        lower_bounds[i, label] = upper_bound
                        bounds_tight = 1

                    # If the condition still holds, then compute the actual
                    # distance between the sample and center. If this is less
                    # than the previous distance, reassign label.
                    if (upper_bound > lower_bounds[i, j]
                        or (upper_bound > center_half_distances[label, j])):

                        distance = _euclidean_dense_dense(
                            &X[i, 0], &centers_old[j, 0], n_features, False)
                        lower_bounds[i, j] = distance
                        if distance < upper_bound:
                            label = j
                            upper_bound = distance

            labels[i] = label
            upper_bounds[i] = upper_bound

        if update_centers:
            weight_in_clusters[label] += sample_weight[i]
            for k in range(n_features):
                centers_new[label * n_features + k] += X[i, k] * sample_weight[i]


def elkan_iter_chunked_sparse(
        X,                                       # IN
        floating[::1] sample_weight,             # IN
        floating[:, ::1] centers_old,            # IN
        floating[:, ::1] centers_new,            # OUT
        floating[::1] weight_in_clusters,        # OUT
        floating[:, ::1] center_half_distances,  # IN
        floating[::1] distance_next_center,      # IN
        floating[::1] upper_bounds,              # INOUT
        floating[:, ::1] lower_bounds,           # INOUT
        int[::1] labels,                         # INOUT
        floating[::1] center_shift,              # OUT
        int n_threads,
        bint update_centers=True):
    """Single iteration of K-means Elkan algorithm with sparse input.

    Update labels and centers (inplace), for one iteration, distributed
    over data chunks.

    Parameters
    ----------
    X : sparse matrix of shape (n_samples, n_features)
        The observations to cluster. Must be in CSR format.

    sample_weight : ndarray of shape (n_samples,), dtype=floating
        The weights for each observation in X.

    centers_old : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : ndarray of shape (n_clusters, n_features), dtype=floating
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.

    weight_in_clusters : ndarray of shape (n_clusters,), dtype=floating
        Placeholder for the sums of the weights of every observation assigned
        to each center.

    center_half_distances : ndarray of shape (n_clusters, n_clusters), \
            dtype=floating
        Half pairwise distances between centers.

    distance_next_center : ndarray of shape (n_clusters,), dtype=floating
        Distance between each center its closest center.

    upper_bounds : ndarray of shape (n_samples,), dtype=floating
        Upper bound for the distance between each sample and its center,
        updated inplace.

    lower_bounds : ndarray of shape (n_samples, n_clusters), dtype=floating
        Lower bound for the distance between each sample and each center,
        updated inplace.

    labels : ndarray of shape (n_samples,), dtype=int
        labels assignment.

    center_shift : ndarray of shape (n_clusters,), dtype=floating
        Distance between old and new centers.

    n_threads : int
        The number of threads to be used by openmp.

    update_centers : bool
        - If True, the labels and the new centers will be computed, i.e. runs
          the E-step and the M-step of the algorithm.
        - If False, only the labels will be computed, i.e runs the E-step of
          the algorithm. This is useful especially when calling predict on a
          fitted model.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_clusters = centers_new.shape[0]

        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        # hard-coded number of samples per chunk. Splitting in chunks is
        # necessary to get parallelism. Chunk size chosen to be same as lloyd's
        int n_samples_chunk = CHUNK_SIZE if n_samples > CHUNK_SIZE else n_samples
        int n_chunks = n_samples // n_samples_chunk
        int n_samples_rem = n_samples % n_samples_chunk
        int chunk_idx, n_samples_chunk_eff
        int start, end

        int i, j, k

        floating[::1] centers_squared_norms = row_norms(centers_old, squared=True)

        floating *centers_new_chunk
        floating *weight_in_clusters_chunk

        openmp.omp_lock_t lock

    # count remainder chunk in total number of chunks
    n_chunks += n_samples != n_chunks * n_samples_chunk

    # number of threads should not be bigger than number of chunks
    n_threads = min(n_threads, n_chunks)

    if update_centers:
        memset(&centers_new[0, 0], 0, n_clusters * n_features * sizeof(floating))
        memset(&weight_in_clusters[0], 0, n_clusters * sizeof(floating))
        openmp.omp_init_lock(&lock)

    with nogil, parallel(num_threads=n_threads):
        # thread local buffers
        centers_new_chunk = <floating*> calloc(n_clusters * n_features, sizeof(floating))
        weight_in_clusters_chunk = <floating*> calloc(n_clusters, sizeof(floating))

        for chunk_idx in prange(n_chunks, schedule='static'):
            start = chunk_idx * n_samples_chunk
            if chunk_idx == n_chunks - 1 and n_samples_rem > 0:
                end = start + n_samples_rem
            else:
                end = start + n_samples_chunk

            _update_chunk_sparse(
                X_data[X_indptr[start]: X_indptr[end]],
                X_indices[X_indptr[start]: X_indptr[end]],
                X_indptr[start: end+1],
                sample_weight[start: end],
                centers_old,
                centers_squared_norms,
                center_half_distances,
                distance_next_center,
                labels[start: end],
                upper_bounds[start: end],
                lower_bounds[start: end],
                centers_new_chunk,
                weight_in_clusters_chunk,
                update_centers)

        # reduction from local buffers. The lock is necessary for that to avoid
        # race conditions.
        if update_centers:
            openmp.omp_set_lock(&lock)
            for j in range(n_clusters):
                weight_in_clusters[j] += weight_in_clusters_chunk[j]
                for k in range(n_features):
                    centers_new[j, k] += centers_new_chunk[j * n_features + k]
            openmp.omp_unset_lock(&lock)

        free(centers_new_chunk)
        free(weight_in_clusters_chunk)

    if update_centers:
        openmp.omp_destroy_lock(&lock)
        _relocate_empty_clusters_sparse(
            X_data, X_indices, X_indptr, sample_weight,
            centers_old, centers_new, weight_in_clusters, labels)

        _average_centers(centers_new, weight_in_clusters)
        _center_shift(centers_old, centers_new, center_shift)

        # update lower and upper bounds
        for i in range(n_samples):
            upper_bounds[i] += center_shift[labels[i]]

            for j in range(n_clusters):
                lower_bounds[i, j] -= center_shift[j]
                if lower_bounds[i, j] < 0:
                    lower_bounds[i, j] = 0


cdef void _update_chunk_sparse(
        floating[::1] X_data,                    # IN
        int[::1] X_indices,                      # IN
        int[::1] X_indptr,                       # IN
        floating[::1] sample_weight,             # IN
        floating[:, ::1] centers_old,            # IN
        floating[::1] centers_squared_norms,     # IN
        floating[:, ::1] center_half_distances,  # IN
        floating[::1] distance_next_center,      # IN
        int[::1] labels,                         # INOUT
        floating[::1] upper_bounds,              # INOUT
        floating[:, ::1] lower_bounds,           # INOUT
        floating *centers_new,                   # OUT
        floating *weight_in_clusters,            # OUT
        bint update_centers) nogil:
    """K-means combined EM step for one sparse data chunk.

    Compute the partial contribution of a single data chunk to the labels and
    centers.
    """
    cdef:
        int n_samples = labels.shape[0]
        int n_clusters = centers_old.shape[0]
        int n_features = centers_old.shape[1]

        floating upper_bound, distance
        int i, j, k, label
        int s = X_indptr[0]

    for i in range(n_samples):
        upper_bound = upper_bounds[i]
        bounds_tight = 0
        label = labels[i]

        # Next center is not far away from the currently assigned center.
        # Sample might need to be assigned to another center.
        if not distance_next_center[label] >= upper_bound:

            for j in range(n_clusters):

                # If this holds, then center_index is a good candidate for the
                # sample to be relabelled, and we need to confirm this by
                # recomputing the upper and lower bounds.
                if (j != label
                    and (upper_bound > lower_bounds[i, j])
                    and (upper_bound > center_half_distances[label, j])):

                    # Recompute upper bound by calculating the actual distance
                    # between the sample and its current assigned center.
                    if not bounds_tight:
                        upper_bound = _euclidean_sparse_dense(
                            X_data[X_indptr[i] - s: X_indptr[i + 1] - s],
                            X_indices[X_indptr[i] - s: X_indptr[i + 1] - s],
                            centers_old[label], centers_squared_norms[label], False)
                        lower_bounds[i, label] = upper_bound
                        bounds_tight = 1

                    # If the condition still holds, then compute the actual
                    # distance between the sample and center. If this is less
                    # than the previous distance, reassign label.
                    if (upper_bound > lower_bounds[i, j]
                        or (upper_bound > center_half_distances[label, j])):
                        distance = _euclidean_sparse_dense(
                            X_data[X_indptr[i] - s: X_indptr[i + 1] - s],
                            X_indices[X_indptr[i] - s: X_indptr[i + 1] - s],
                            centers_old[j], centers_squared_norms[j], False)
                        lower_bounds[i, j] = distance
                        if distance < upper_bound:
                            label = j
                            upper_bound = distance

            labels[i] = label
            upper_bounds[i] = upper_bound

        if update_centers:
            weight_in_clusters[label] += sample_weight[i]
            for k in range(X_indptr[i] - s, X_indptr[i + 1] - s):
                centers_new[label * n_features + X_indices[k]] += X_data[k] * sample_weight[i]
