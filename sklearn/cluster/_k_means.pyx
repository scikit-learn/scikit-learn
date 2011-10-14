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


cdef inline DOUBLE array_ddot(int n,
                np.ndarray[DOUBLE, ndim=1] a,
                np.ndarray[DOUBLE, ndim=1] b):
    return ddot(n, <DOUBLE*>(a.data), 1, <DOUBLE*>(b.data), 1)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _assign_labels_array(np.ndarray[DOUBLE, ndim=2] X,
                         np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                         slice_, np.ndarray[DOUBLE, ndim=2] centers,
                         np.ndarray[INT, ndim=1] labels):
    """Compute label assignement and inertia for a slice of a dense array"""
    cdef:
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]
        int n_samples = slice_.stop - slice_.start
        int slice_start = slice_.start
        int sample_idx, center_idx, feature_idx
        int i
        DOUBLE inertia = 0.0
        DOUBLE min_dist
        DOUBLE dist
        np.ndarray[DOUBLE, ndim=1] center_squared_norms = np.zeros(
            n_clusters, dtype=np.float64)

    for center_idx in range(n_clusters):
        center_squared_norms[center_idx] = array_ddot(
            n_features, centers[center_idx], centers[center_idx])

    for i in range(n_samples):
        sample_idx = slice_start + i
        min_dist = -1
        for center_idx in range(n_clusters):
            dist = 0.0
            # hardcoded: minimize euclidean distance to cluster center:
            # ||a - b||^2 = ||a||^2 + ||b||^2 -2 <a, b>
            dist += array_ddot(n_features, X[sample_idx], centers[center_idx])
            dist *= -2
            dist += center_squared_norms[center_idx]
            dist += x_squared_norms[sample_idx]
            if min_dist < 0.0 or dist < min_dist:
                min_dist = dist
                labels[i] = center_idx
        inertia += min_dist

    return inertia

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _assign_labels_csr(X, np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                       slice_, np.ndarray[DOUBLE, ndim=2] centers,
                       np.ndarray[INT, ndim=1] labels):
    """Compute label assignement and inertia for a slice of a CSR input"""
    cdef:
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[INT, ndim=1] X_indices = X.indices
        np.ndarray[INT, ndim=1] X_indptr = X.indptr
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]
        int n_samples = slice_.stop - slice_.start
        int slice_start = slice_.start
        int sample_idx, center_idx, feature_idx
        int i, k
        DOUBLE inertia = 0.0
        DOUBLE min_dist
        DOUBLE dist
        np.ndarray[DOUBLE, ndim=1] center_squared_norms = np.zeros(
            n_clusters, dtype=np.float64)

    for center_idx in range(n_clusters):
        center_squared_norms[center_idx] = array_ddot(
            n_features, centers[center_idx], centers[center_idx])

    for i in range(n_samples):
        sample_idx = slice_start + i
        min_dist = -1
        for center_idx in range(n_clusters):
            dist = 0.0
            # hardcoded: minimize euclidean distance to cluster center:
            # ||a - b||^2 = ||a||^2 + ||b||^2 -2 <a, b>
            for k in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                dist += centers[center_idx, X_indices[k]] * X_data[k]
            dist *= -2
            dist += center_squared_norms[center_idx]
            dist += x_squared_norms[sample_idx]
            if min_dist < 0.0 or dist < min_dist:
                min_dist = dist
                labels[i] = center_idx
        inertia += min_dist

    return inertia


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _mini_batch_update_csr(X, np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                           batch_slice, np.ndarray[DOUBLE, ndim=2] centers,
                           np.ndarray[INT, ndim=1] counts):
    """Incremental update of the centers for sparse MiniBatchKMeans.

    Parameters
    ----------

    X: CSR matrix, dtype float64
        The complete (pre allocated) training set as a CSR matrix.

    batch_slice: slice
        The row slice of the mini batch.

    centers: array, shape (n_clusters, n_features)
        The cluster centers

    counts: array, shape (n_clusters,)
         The vector in which we keep track of the numbers of elements in a
         cluster

    """
    cdef:
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[INT, ndim=1] X_indices = X.indices
        np.ndarray[INT, ndim=1] X_indptr = X.indptr
        int n_samples = batch_slice.stop - batch_slice.start
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]

        DOUBLE center_squared_norm = 0.0

        int batch_slice_start = batch_slice.start
        int batch_slice_stop = batch_slice.stop

        int sample_idx, center_idx, feature_idx
        int i, k
        int old_count, batch_count

        # TODO: reuse a array preallocated outside of the mini batch main loop
        np.ndarray[INT, ndim=1] nearest_center = np.zeros(
            n_samples, dtype=np.int32)

    # step 1: assign minibatch samples to there nearest center
    _assign_labels_csr(X, x_squared_norms, batch_slice, centers,
                       nearest_center)

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
    cdef:
        unsigned int n_samples = X.shape[0]
        unsigned int n_features = X.shape[1]
        np.ndarray[DOUBLE, ndim=1] norms = np.zeros((n_samples,),
                                                    dtype=np.float64)
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[INT, ndim=1] X_indices = X.indices
        np.ndarray[INT, ndim=1] X_indptr = X.indptr

        unsigned int i
        unsigned int j
        double sum_
        int withsqrt = not squared

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        if withsqrt:
            sum_ = sqrt(sum_)

        norms[i] = sum_
    return norms
