# cython: profile=True, boundscheck=False, wraparound=False, cdivision=True
# cython: language_level=3
#
# Author: Andreas Mueller
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
cimport cython
cimport openmp
from cython cimport floating
from cython.parallel import prange, parallel
from libc.math cimport sqrt
from libc.stdlib cimport calloc, free
from libc.string cimport memset, memcpy

from ..metrics import euclidean_distances
from ._k_means cimport _relocate_empty_clusters_dense, _mean_and_center_shift


np.import_array()


cdef floating euclidean_dist(floating* a, floating* b, int n_features) nogil:
    """Euclidean distance between a and b, optimized for vectorization"""
    cdef:
        int i
        int n = n_features // 4
        int rem = n_features % 4
        floating result = 0

    for i in range(n):
        result += ((a[0] - b[0]) * (a[0] - b[0])
                  +(a[1] - b[1]) * (a[1] - b[1])
                  +(a[2] - b[2]) * (a[2] - b[2])
                  +(a[3] - b[3]) * (a[3] - b[3]))
        a += 4; b += 4

    for i in range(rem):
        result += (a[i] - b[i]) * (a[i] - b[i])

    return sqrt(result)


cpdef _init_bounds(np.ndarray[floating, ndim=2, mode='c'] X,
                   floating[:, ::1] centers,
                   floating[:, ::1] center_half_distances,
                   int[::1] labels,
                   floating[::1] upper_bounds,
                   floating[:, ::1] lower_bounds):
    """Initialize upper and lower bounds for each sample.

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
    X : nd-array, shape (n_samples, n_features)
        The input data.

    centers : nd-array, shape (n_clusters, n_features)
        The cluster centers.

    center_half_distances : nd-array, shape (n_clusters, n_clusters)
        The half of the distance between any 2 clusters centers.

    labels : nd-array, shape(n_samples)
        The label for each sample. This array is modified in place.

    lower_bounds : nd-array, shape(n_samples, n_clusters)
        The lower bound on the distance between a sample and each cluster
        center. It is modified in place.

    upper_bounds : nd-array, shape(n_samples,)
        The distance of each sample from its closest cluster center.  This is
        modified in place by the function.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_clusters = centers.shape[0]
        int n_features = X.shape[1]

        floating min_dist, dist
        int best_cluster, i, j

    for i in range(n_samples):
        best_cluster = 0
        min_dist = euclidean_dist(&X[i, 0], &centers[0, 0], n_features)
        lower_bounds[i, 0] = min_dist
        for j in range(1, n_clusters):
            if min_dist > center_half_distances[best_cluster, j]:
                dist = euclidean_dist(&X[i, 0], &centers[j, 0], n_features)
                lower_bounds[i, j] = dist
                if dist < min_dist:
                    min_dist = dist
                    best_cluster = j
        labels[i] = best_cluster
        upper_bounds[i] = min_dist


cpdef void _elkan_iter_chunked_dense(np.ndarray[floating, ndim=2, mode='c'] X,
                                     floating[::1] sample_weight,
                                     floating[:, ::1] centers_old,
                                     floating[:, ::1] centers_new,
                                     floating[::1] weight_in_clusters,
                                     floating[:, ::1] center_half_distances,
                                     floating[::1] distance_next_center,
                                     floating[::1] upper_bounds,
                                     floating[:, ::1] lower_bounds,
                                     int[::1] labels,
                                     floating[::1] center_shift,
                                     int n_jobs = -1,
                                     bint update_centers = True):
    """Single interation of K-means elkan algorithm

    Update labels and centers (inplace), for one iteration, distributed
    over data chunks.

    Parameters
    ----------
    X : {float32, float64} array-like, shape (n_samples, n_features)
        The observations to cluster.

    sample_weight : {float32, float64} array-like, shape (n_samples,)
        The weights for each observation in X.

    centers_old : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.

    weight_in_clusters : {float32, float64} array-like, shape (n_clusters,)
        Placeholder for the sums of the weights of every observation assigned
        to each center.

    center_half_distances : {float32, float64} array-like, \
shape (n_clusters, n_clusters)
        Half pairwise distances between centers.

    distance_next_center : {float32, float64} array-like, shape (n_clusters,)
        Distance between each center it's closest center.
    
    upper_bounds : {float32, float64} array-like, shape (n_samples,)
        Upper bound for the distance between each sample and it's center,
        updated inplace.

    lower_bounds : {float32, float64} array-like, shape (n_samples, n_clusters)
        Lower bound for the distance between each sample and each center,
        updated inplace.

    labels : int array-like, shape (n_samples,)
        labels assignment.

    center_shift : {float32, float64} array-like, shape (n_clusters,)
        Distance between old and new centers.

    n_jobs : int
        The number of threads to be used by openmp. If -1, openmp will use as
        many as possible.

    update_centers : bool
        - If True, the labels and the new centers will be computed, i.e. runs
          the E-step and the M-step of the algorithm.
        - If False, only the labels will be computed, i.e runs the E-step of
          the algorithm.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_clusters = centers_new.shape[0]

        # hard-coded number of samples per chunk. Splitting in chunks is
        # necessary to get parallelism. Chunk size chosed to be same as lloyd's
        int n_samples_chunk = 256 if n_samples > 256 else n_samples
        int n_chunks = n_samples // n_samples_chunk
        int n_samples_r = n_samples % n_samples_chunk
        int chunk_idx, n_samples_chunk_eff
        int num_threads

        int i, j, k
        int label
        floating alpha, tmp, x

        floating *centers_new_chunk
        floating *weight_in_clusters_chunk

    # count remainder chunk in total number of chunks
    n_chunks += n_samples != n_chunks * n_samples_chunk

    # re-initialize all arrays at each iteration
    if update_centers:
        memcpy(&centers_old[0, 0], &centers_new[0, 0], n_clusters * n_features * sizeof(floating))
        memset(&centers_new[0, 0], 0, n_clusters * n_features * sizeof(floating))
        memset(&weight_in_clusters[0], 0, n_clusters * sizeof(floating))

    # set number of threads to be used by openmp
    num_threads = n_jobs if n_jobs != -1 else openmp.omp_get_max_threads()

    with nogil, parallel(num_threads=num_threads):
        # thread local buffers
        centers_new_chunk = <floating*> calloc(n_clusters * n_features, sizeof(floating))
        weight_in_clusters_chunk = <floating*> calloc(n_clusters, sizeof(floating))
        
        for chunk_idx in prange(n_chunks):
            if n_samples_r > 0 and chunk_idx == n_chunks - 1:
                n_samples_chunk_eff = n_samples_r
            else:
                n_samples_chunk_eff = n_samples_chunk

            _update_chunk(
                &X[chunk_idx * n_samples_chunk, 0],
                &sample_weight[chunk_idx * n_samples_chunk],
                &centers_old[0, 0],
                centers_new_chunk,
                &center_half_distances[0, 0],
                &distance_next_center[0],
                weight_in_clusters_chunk,
                &labels[chunk_idx * n_samples_chunk],
                &upper_bounds[chunk_idx * n_samples_chunk],
                &lower_bounds[chunk_idx * n_samples_chunk, 0],
                n_samples_chunk_eff,
                n_clusters,
                n_features,
                update_centers)

        # reduction from local buffers. The gil is necessary for that to avoid
        # race conditions.
        if update_centers:
            with gil:
                for j in range(n_clusters):
                    weight_in_clusters[j] += weight_in_clusters_chunk[j]
                    for k in range(n_features):
                        centers_new[j, k] += centers_new_chunk[j * n_features + k]

        free(weight_in_clusters_chunk)
        free(centers_new_chunk)

    if update_centers:
        _relocate_empty_clusters_dense(X, sample_weight, centers_new,
                                       weight_in_clusters, labels)

        _mean_and_center_shift(centers_old, centers_new, weight_in_clusters,
                               center_shift)

        # update lower and upper bounds
        for i in range(n_samples):
            upper_bounds[i] += center_shift[labels[i]]

            for j in range(n_clusters):
                lower_bounds[i, j] -= center_shift[j]
                if lower_bounds[i, j] < 0:
                    lower_bounds[i, j] = 0


cdef void _update_chunk(floating *X,
                        floating *sample_weight,
                        floating *centers_old,
                        floating *centers_new,
                        floating *center_half_distances,
                        floating *distance_next_center,
                        floating *weight_in_clusters, 
                        int *labels,
                        floating *upper_bounds,
                        floating *lower_bounds,
                        int n_samples,
                        int n_clusters,
                        int n_features,
                        bint update_centers) nogil:
    """K-means step for one data chunk using elkan algorithm
    
    Compute the partial contribution of a single data chunk to the labels and
    centers.
    """
    cdef:
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
                    and (upper_bound > lower_bounds[i * n_clusters + j])
                    and (upper_bound > center_half_distances[label * n_clusters + j])):

                    # Recompute upper bound by calculating the actual distance
                    # between the sample and it's current assigned center.
                    if not bounds_tight:
                        upper_bound = euclidean_dist(X + i * n_features,
                                                     centers_old + label * n_features,
                                                     n_features)
                        lower_bounds[i * n_clusters + label] = upper_bound
                        bounds_tight = 1

                    # If the condition still holds, then compute the actual
                    # distance between the sample and center. If this is less
                    #than the previous distance, reassign label.
                    if (upper_bound > lower_bounds[i * n_clusters + j]
                        or (upper_bound > center_half_distances[label * n_clusters + j])):

                        distance = euclidean_dist(X + i * n_features,
                                                  centers_old + j * n_features,
                                                  n_features)
                        lower_bounds[i * n_clusters + j] = distance
                        if distance < upper_bound:
                            label = j
                            upper_bound = distance

            labels[i] = label
            upper_bounds[i] = upper_bound

        if update_centers:
            weight_in_clusters[label] += sample_weight[i]
            for k in range(n_features):  
                centers_new[label * n_features + k] += X[i * n_features + k] * sample_weight[i]