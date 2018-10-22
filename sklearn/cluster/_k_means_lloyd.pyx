# cython: profile=True, boundscheck=False, wraparound=False, cdivision=True
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
cimport cython
cimport openmp
from cython cimport floating
from cython.parallel import prange, parallel
from scipy.linalg.cython_blas cimport sgemm, dgemm
from libc.math cimport sqrt
from libc.stdlib cimport malloc, free
from libc.string cimport memset, memcpy

from ._k_means import (_relocate_empty_clusters_dense,
                       _relocate_empty_clusters_sparse,
                       _mean_and_center_shift)


np.import_array()


cdef:
    float MAX_FLT = np.finfo(np.float32).max
    double MAX_DBL = np.finfo(np.float64).max


cdef void xgemm(char *ta, char *tb, int *m, int *n, int *k, floating *alpha,
                floating *A, int *lda, floating *B, int *ldb, floating *beta,
                floating *C, int *ldc) nogil:
    if floating is float:
        sgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    else:
        dgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)


cpdef void _lloyd_iter_chunked_dense(np.ndarray[floating, ndim=2, mode='c'] X,
                                     floating[::1] sample_weight,
                                     floating[::1] x_squared_norms,
                                     floating[:, ::1] centers_old,
                                     floating[:, ::1] centers_new,
                                     floating[::1] centers_squared_norms,
                                     floating[::1] weight_in_clusters, 
                                     int[::1] labels,
                                     floating[::1] center_shift,
                                     int n_jobs = -1,
                                     bint update_centers = True):
    """Single interation of K-means lloyd algorithm

    Update labels and centers (inplace), for one iteration, distributed
    over data chunks.

    Parameters
    ----------
    X : {float32, float64} array-like, shape (n_samples, n_features)
        The observations to cluster.

    sample_weight : {float32, float64} array-like, shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : {float32, float64} array-like, shape (n_samples,)
        Squared L2 norm of X.
    
    centers_old : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.
    
    centers_squared_norms : {float32, float64} array-like, shape (n_clusters,)
        Squared L2 norm of the centers.

    weight_in_clusters : {float32, float64} array-like, shape (n_clusters,)
        Placeholder for the sums of the weights of every observation assigned
        to each center.

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

        # hard-coded number of samples per chunk. Appeared to be close to
        # optimal in all situations.
        int n_samples_chunk = 256 if n_samples > 256 else n_samples
        int n_chunks = n_samples // n_samples_chunk
        int n_samples_r = n_samples % n_samples_chunk
        int chunk_idx, n_samples_chunk_eff
        int num_threads

        int j, k
        floating alpha

        floating *centers_new_chunk
        floating *weight_in_clusters_chunk
        floating *pairwise_distances_chunk

    # count remainder chunk in total number of chunks
    n_chunks += n_samples != n_chunks * n_samples_chunk
    
    # re-initialize all arrays at each iteration
    memset(&centers_squared_norms[0], 0, n_clusters * sizeof(floating))
    for j in xrange(n_clusters):
        for k in xrange(n_features):
            centers_squared_norms[j] += centers_new[j, k] * centers_new[j, k]

    if update_centers:
        memcpy(&centers_old[0, 0], &centers_new[0, 0],
            n_clusters * n_features * sizeof(floating))
        memset(&centers_new[0, 0], 0,
            n_clusters * n_features * sizeof(floating))
        memset(&weight_in_clusters[0], 0, n_clusters * sizeof(floating))

    # set number of threads to be used by openmp
    num_threads = n_jobs if n_jobs != -1 else openmp.omp_get_max_threads()
    with nogil, parallel(num_threads=num_threads):
        centers_new_chunk = \
            <floating*> malloc(n_clusters * n_features * sizeof(floating))

        weight_in_clusters_chunk = \
            <floating*> malloc(n_clusters * sizeof(floating))

        pairwise_distances_chunk = \
            <floating*> malloc(n_samples_chunk * n_clusters * sizeof(floating))

        # initialize local buffers
        memset(centers_new_chunk, 0,
               n_clusters * n_features * sizeof(floating))
        memset(weight_in_clusters_chunk, 0, n_clusters * sizeof(floating))
        
        for chunk_idx in prange(n_chunks):
            if n_samples_r > 0 and chunk_idx == n_chunks - 1:
                n_samples_chunk_eff = n_samples_r
            else:
                n_samples_chunk_eff = n_samples_chunk

            _update_chunk_dense(
                &X[chunk_idx * n_samples_chunk, 0],
                &sample_weight[chunk_idx * n_samples_chunk],
                &x_squared_norms[chunk_idx * n_samples_chunk],
                &centers_old[0, 0],
                centers_new_chunk,
                &centers_squared_norms[0],
                weight_in_clusters_chunk,
                pairwise_distances_chunk,
                &labels[chunk_idx * n_samples_chunk],
                n_samples_chunk_eff,
                n_clusters,
                n_features,
                update_centers)

        # reduction from local buffers. The gil is necessary for that to avoid
        # race conditions.
        if update_centers:
            with gil:
                for j in xrange(n_clusters):
                    weight_in_clusters[j] += weight_in_clusters_chunk[j]
                    for k in xrange(n_features):
                        centers_new[j, k] += \
                            centers_new_chunk[j * n_features + k]

        free(weight_in_clusters_chunk)
        free(centers_new_chunk)
        free(pairwise_distances_chunk)

    if update_centers:
        _relocate_empty_clusters_dense(X, sample_weight, centers_new,
                                       weight_in_clusters, labels)

        _mean_and_center_shift(centers_old, centers_new, weight_in_clusters,
                               center_shift)


cdef void _update_chunk_dense(floating *X,
                              floating *sample_weight,
                              floating *x_squared_norms,
                              floating *centers_old,
                              floating *centers_new,
                              floating *centers_squared_norms,
                              floating *weight_in_clusters,
                              floating *pairwise_distances,
                              int *labels,
                              int n_samples,
                              int n_clusters,
                              int n_features,
                              bint update_centers) nogil:
    """K-means combined EM step for one data chunk
    
    Compute the partial contribution of a single data chunk to the labels and
    centers.
    """
    cdef:
        floating sq_dist, min_sq_dist
        int i, j, k, best_cluster
    
        # parameters for the BLAS gemm
        floating alpha = -2.0
        floating beta = 1.0
        char *trans_data = 'n'
        char *trans_centers = 't'

    # Instead of computing the full pairwise squared distances matrix,
    # ||X - C||² = ||X||² - 2 X.C^T + ||C||², we only need to store
    # the - 2 X.C^T + ||C||² term since the argmin for a given sample only
    # depends on the centers.
    for i in xrange(n_samples):
        for j in xrange(n_clusters):
            pairwise_distances[i * n_clusters + j] = centers_squared_norms[j]
    
    xgemm(trans_centers, trans_data, &n_clusters, &n_samples, &n_features,
          &alpha, centers_old, &n_features, X, &n_features,
          &beta, pairwise_distances, &n_clusters)

    for i in xrange(n_samples):
        min_sq_dist = pairwise_distances[i * n_clusters]
        best_cluster = 0
        for j in xrange(n_clusters):
            sq_dist = pairwise_distances[i * n_clusters + j]
            if sq_dist < min_sq_dist:
                min_sq_dist = sq_dist
                best_cluster = j

        labels[i] = best_cluster

        if update_centers:
            weight_in_clusters[best_cluster] += sample_weight[i]
            for k in xrange(n_features):  
                centers_new[best_cluster * n_features + k] += \
                    X[i * n_features + k] * sample_weight[i]


cpdef void _lloyd_iter_chunked_sparse(X,
                                      floating[::1] sample_weight,
                                      floating[::1] x_squared_norms,
                                      floating[:, ::1] centers_old,
                                      floating[:, ::1] centers_new,
                                      floating[::1] centers_squared_norms,
                                      floating[::1] weight_in_clusters, 
                                      int[::1] labels,
                                      floating[::1] center_shift,
                                      int n_jobs = -1,
                                      bint update_centers = True):
    """Single interation of K-means lloyd algorithm

    Update labels and centers (inplace), for one iteration, distributed
    over data chunks.

    Parameters
    ----------
    X : {float32, float64} CSR matrix, shape (n_samples, n_features)
        The observations to cluster.

    sample_weight : {float32, float64} array-like, shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : {float32, float64} array-like, shape (n_samples,)
        Squared L2 norm of X.
    
    centers_old : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers before previous iteration, placeholder for the centers after
        previous iteration.

    centers_new : {float32, float64} array-like, shape (n_clusters, n_features)
        Centers after previous iteration, placeholder for the new centers
        computed during this iteration.
    
    centers_squared_norms : {float32, float64} array-like, shape (n_clusters,)
        Squared L2 norm of the centers.

    weight_in_clusters : {float32, float64} array-like, shape (n_clusters,)
        Placeholder for the sums of the weights of every observation assigned
        to each center.

    labels : int array-like, shape (n_samples,)
        labels assignment.
    
    center_shift : {float32, float64} array-like, shape (n_clusters,)
        Distance between old and new centers.

    n_jobs : int
        The number of threads to be used by openmp. If -1, openmp will use as
        many as possible.

    update_centers : bool
        - If True, the labels and the new centers will be computed.
        - If False, only the labels will be computed.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_clusters = centers_new.shape[0]

        # Chosed same as for dense. Does not have the same impact since with
        # sparse data the pairwise distances matrix is not precomputed.
        # However, splitting in chunks is necessary to get parallelism.
        int n_samples_chunk = 256 if n_samples > 256 else n_samples
        int n_chunks = n_samples // n_samples_chunk
        int n_samples_r = n_samples % n_samples_chunk
        int chunk_idx, n_samples_chunk_eff
        int num_threads

        int j, k
        floating alpha

        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        floating *centers_new_chunk
        floating *weight_in_clusters_chunk

    # count remainder for total number of chunks
    n_chunks += n_samples != n_chunks * n_samples_chunk
    
    # re-initialize all arrays at each iteration
    memset(&centers_squared_norms[0], 0, n_clusters * sizeof(floating))
    for j in xrange(n_clusters):
        for k in xrange(n_features):
            centers_squared_norms[j] += centers_new[j, k] * centers_new[j, k]

    if update_centers:
        memcpy(&centers_old[0, 0], &centers_new[0, 0],
               n_clusters * n_features * sizeof(floating))
        memset(&centers_new[0, 0], 0, 
               n_clusters * n_features * sizeof(floating))
        memset(&weight_in_clusters[0], 0, n_clusters * sizeof(floating))

    # set number of threads to be used by openmp
    num_threads = n_jobs if n_jobs != -1 else openmp.omp_get_max_threads()
    with nogil, parallel(num_threads=num_threads):
        centers_new_chunk = \
            <floating*> malloc(n_clusters * n_features * sizeof(floating))

        weight_in_clusters_chunk = \
            <floating*> malloc(n_clusters * sizeof(floating))

        # initialize local buffers
        memset(centers_new_chunk, 0,
               n_clusters * n_features * sizeof(floating))
        memset(weight_in_clusters_chunk, 0, n_clusters * sizeof(floating))

        for chunk_idx in prange(n_chunks):
            if n_samples_r > 0 and chunk_idx == n_chunks - 1:
                n_samples_chunk_eff = n_samples_r
            else:
                n_samples_chunk_eff = n_samples_chunk

            _update_chunk_sparse(
                &X_data[X_indptr[chunk_idx * n_samples_chunk]],
                &X_indices[X_indptr[chunk_idx * n_samples_chunk]],
                &X_indptr[chunk_idx * n_samples_chunk],
                &sample_weight[chunk_idx * n_samples_chunk],
                &x_squared_norms[chunk_idx * n_samples_chunk],
                &centers_old[0, 0],
                centers_new_chunk,
                &centers_squared_norms[0],
                weight_in_clusters_chunk,
                &labels[chunk_idx * n_samples_chunk],
                n_samples_chunk_eff,
                n_clusters,
                n_features,
                update_centers)

        # reduction from local buffers. The gil is necessary for that to avoid
        # race conditions.
        if update_centers:
            with gil:
                for j in xrange(n_clusters):
                    weight_in_clusters[j] += weight_in_clusters_chunk[j]
                    for k in xrange(n_features):
                        centers_new[j, k] += \
                            centers_new_chunk[j * n_features + k]

        free(weight_in_clusters_chunk)
        free(centers_new_chunk)

    if update_centers:
        _relocate_empty_clusters_sparse(X_data, X_indices, X_indptr,
                                        sample_weight, centers_new,
                                        weight_in_clusters, labels)

        _mean_and_center_shift(centers_old, centers_new, weight_in_clusters,
                               center_shift)


cdef void _update_chunk_sparse(floating *X_data,
                               int *X_indices,
                               int *X_indptr,
                               floating *sample_weight,
                               floating *x_squared_norms,
                               floating *centers_old,
                               floating *centers_new,
                               floating *centers_squared_norms,
                               floating *weight_in_cluster,
                               int *labels,
                               int n_samples,
                               int n_clusters,
                               int n_features,
                               bint update_centers) nogil:
    """K-means combined EM step for one data chunk
    
    Compute the partial contribution of a single data chunk to the labels and
    centers.
    """
    cdef:    
        floating sq_dist, min_sq_dist
        int i, j, k, best_cluster
        floating max_floating = MAX_FLT if floating is float else MAX_DBL
        int s = X_indptr[0]

    # XXX Precompute the pairwise distances matrix is not worth for sparse
    # currently. Should be tested when BLAS (sparse x dense) matrix
    # multiplication is available.
    for i in xrange(n_samples):
        min_sq_dist = max_floating
        best_cluster = 0

        for j in xrange(n_clusters):
            sq_dist = 0.0
            for k in xrange(X_indptr[i] - s, X_indptr[i + 1] - s):
                sq_dist += \
                    centers_old[j * n_features + X_indices[k]] * X_data[k]
            
            # Instead of computing the full squared distance with each cluster,
            # ||X - C||² = ||X||² - 2 X.C^T + ||C||², we only need to compute
            # the - 2 X.C^T + ||C||² term since the argmin for a given sample
            # only depends on the centers C.
            sq_dist = centers_squared_norms[j] -2 * sq_dist
            if sq_dist < min_sq_dist:
                min_sq_dist = sq_dist
                best_cluster = j
    
        labels[i] = best_cluster
        
        if update_centers:
            weight_in_cluster[best_cluster] += sample_weight[i]
            for k in xrange(X_indptr[i] - s, X_indptr[i + 1] - s):
                centers_new[best_cluster * n_features + X_indices[k]] += \
                    X_data[k] * sample_weight[i]
