# cython: profile=True
# Profiling is enabled by default as the overhead does not seem to be measurable
# on this specific use case.

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

import numpy as np
from ..utils.extmath import norm
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INT

cdef extern from "math.h":
    double sqrt(double f)

cdef extern from "cblas.h":
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _assign_labels(np.ndarray[DOUBLE, ndim=1] X_data,
                   np.ndarray[INT, ndim=1] X_indices,
                   np.ndarray[INT, ndim=1] X_indptr,
                   np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                   batch_slice,
                   np.ndarray[DOUBLE, ndim=2] centers,
                   np.ndarray[INT, ndim=1] labels):
    cdef int n_clusters = centers.shape[0]
    cdef int n_features = centers.shape[1]
    cdef int n_samples = batch_slice.stop - batch_slice.start
    cdef int batch_slice_start = batch_slice.start
    cdef int sample_idx = -1
    cdef int c = 0
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef DOUBLE min_dist = 0.0
    cdef DOUBLE dist = 0.0

    cdef np.ndarray[DOUBLE, ndim=1] center_squared_norms = np.zeros(
        n_clusters, dtype=np.float64)

    for c in range(n_clusters):
        center_squared_norms[c] = ddot(
            n_features,
            <DOUBLE*>(centers.data + c * n_features * sizeof(DOUBLE)), 1,
            <DOUBLE*>(centers.data + c * n_features * sizeof(DOUBLE)), 1)

    for i in range(n_samples):
        sample_idx = batch_slice_start + i
        min_dist = -1
        for j in range(n_clusters):
            dist = 0.0
            # hardcoded: minimize euclidean distance to cluster center:
            # ||a - b||^2 = ||a||^2 + ||b||^2 -2 <a, b>
            for k in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                dist += centers[j, X_indices[k]] * X_data[k]
            dist *= -2
            dist += center_squared_norms[j] + x_squared_norms[sample_idx]
            if min_dist < 0.0 or dist < min_dist:
                min_dist = dist
                labels[i] = j


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _mini_batch_update_sparse(np.ndarray[DOUBLE, ndim=1] X_data,
                              np.ndarray[INT, ndim=1] X_indices,
                              np.ndarray[INT, ndim=1] X_indptr,
                              np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                              batch_slice,
                              np.ndarray[DOUBLE, ndim=2] centers,
                              np.ndarray[INT, ndim=1] counts):
    """Incremental update of the centers for sparse MiniBatchKMeans.

    Parameters
    ----------

    X_data: array, dtype float64
        Data array of CSR matrix.

    X_indices: array, dtype int32
        Indices array of CSR matrix.

    X_indptr: array, dtype int32
        index pointer array of CSR matrix.

    batch_slice: slice
        The row slice of the mini batch.

    centers: array, shape (n_clusters, n_features)
        The cluster centers

    counts: array, shape (n_clusters,)
         The vector in which we keep track of the numbers of elements in a
         cluster

    """
    cdef int n_samples = batch_slice.stop - batch_slice.start
    cdef int n_clusters = centers.shape[0]
    cdef int n_features = centers.shape[1]

    cdef DOUBLE center_squared_norm = 0.0

    # TODO: replace slice by permutation array view to avoid memory copy
    # of X_shuffled in main minibatch function
    cdef int batch_slice_start = batch_slice.start
    cdef int batch_slice_stop = batch_slice.stop
    cdef int sample_idx = 0
    cdef int center_idx = 0
    cdef int feature_idx = 0
    cdef int i = 0
    cdef int k = 0
    cdef int old_count = 0
    cdef int batch_count = 0

    # TODO: reuse a array preallocated outside of the mini batch main loop
    cdef np.ndarray[INT, ndim=1] nearest_center = np.zeros(
        n_samples, dtype=np.int32)

    # step 1: assign minibatch samples to there nearest center
    _assign_labels(X_data, X_indices, X_indptr, x_squared_norms, batch_slice,
                   centers, nearest_center)

    # step 2: move centers to mean of old and newly assigned samples
    for center_idx in range(n_clusters):
        old_count = counts[center_idx]
        if old_count > 0:
            for feature_idx in range(n_features):
                # inplace remove previous count scaling
                centers[center_idx, feature_idx] *= old_count

        # iterate of over samples assigned to this cluster to move the center
        # location by inplace summation
        batch_count = 0
        for i in range(n_samples):
            if nearest_center[i] != center_idx:
                continue
            sample_idx = batch_slice_start + i
            batch_count += 1

            # inplace sum with new samples that are members of this cluster
            for k in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                centers[center_idx, X_indices[k]] += X_data[k]

        # inplace rescale center with updated count
        if old_count + batch_count > 0:
            for feature_idx in range(n_features):
                centers[center_idx, feature_idx] /= old_count + batch_count

            # update the count statistics for this center
            counts[center_idx] += batch_count


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csr_row_norm_l2(X, squared=True):
    """Get L2 norm of each row in X."""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]
    cdef np.ndarray[DOUBLE, ndim=1] norms = np.zeros((n_samples,),
                                                     dtype=np.float64)

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[INT, ndim=1] X_indices = X.indices
    cdef np.ndarray[INT, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_
    cdef int withsqrt = not squared

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        if withsqrt:
            sum_ = sqrt(sum_)

        norms[i] = sum_
    return norms
