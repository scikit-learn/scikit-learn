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


@cython.profile(False)
@cython.wraparound(False)
cdef inline DOUBLE array_ddot(int n,
                np.ndarray[DOUBLE, ndim=2] a, int a_idx,
                np.ndarray[DOUBLE, ndim=2] b, int b_idx):
    return ddot(n, <DOUBLE*>(a.data + a_idx * n * sizeof(DOUBLE)), 1,
                <DOUBLE*>(b.data + b_idx * n * sizeof(DOUBLE)), 1)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef DOUBLE _assign_labels_array(np.ndarray[DOUBLE, ndim=2] X,
                                  np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                                  unsigned int slice_start,
                                  unsigned int slice_stop,
                                  np.ndarray[DOUBLE, ndim=2] centers,
                                  np.ndarray[INT, ndim=1] labels):
    """Compute label assignement and inertia for a slice of a dense array

    Return the inertia (sum of squared distances to the centers).
    """
    cdef:
        unsigned int n_clusters = centers.shape[0]
        unsigned int n_features = centers.shape[1]
        unsigned int n_samples = slice_stop - slice_start
        unsigned int sample_idx, center_idx, feature_idx
        unsigned int i
        DOUBLE inertia = 0.0
        DOUBLE min_dist
        DOUBLE dist
        np.ndarray[DOUBLE, ndim=1] center_squared_norms = np.zeros(
            n_clusters, dtype=np.float64)

    for center_idx in range(n_clusters):
        center_squared_norms[center_idx] = array_ddot(
            n_features, centers, center_idx, centers, center_idx)

    for i in range(n_samples):
        sample_idx = slice_start + i
        min_dist = -1
        for center_idx in range(n_clusters):
            dist = 0.0
            # hardcoded: minimize euclidean distance to cluster center:
            # ||a - b||^2 = ||a||^2 + ||b||^2 -2 <a, b>
            dist += array_ddot(n_features, X, sample_idx, centers, center_idx)
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
cpdef DOUBLE _assign_labels_csr(X, np.ndarray[DOUBLE, ndim=1] x_squared_norms,
                                unsigned int slice_start,
                                unsigned int slice_stop,
                                np.ndarray[DOUBLE, ndim=2] centers,
                                np.ndarray[INT, ndim=1] labels):
    """Compute label assignement and inertia for a slice of a CSR input

    Return the inertia (sum of squared distances to the centers).
    """
    cdef:
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[INT, ndim=1] X_indices = X.indices
        np.ndarray[INT, ndim=1] X_indptr = X.indptr
        unsigned int n_clusters = centers.shape[0]
        unsigned int n_features = centers.shape[1]
        unsigned int n_samples = slice_stop - slice_start
        unsigned int sample_idx, center_idx, feature_idx
        unsigned int i, k
        DOUBLE inertia = 0.0
        DOUBLE min_dist
        DOUBLE dist
        np.ndarray[DOUBLE, ndim=1] center_squared_norms = np.zeros(
            n_clusters, dtype=np.float64)

    for center_idx in range(n_clusters):
        center_squared_norms[center_idx] = array_ddot(
            n_features, centers, center_idx, centers, center_idx)

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

    Return
    ------
    inertia: float
        The inertia of the batch prior to centers update, i.e. the sum distances
        to the closest center for each sample in the slice. This is the
        objective function being minimized by the k-means algorithm.

    squared_diff: float
        The sum of squared update (squared norm of the centers position change).

    Both squared diff and inertia are commonly used to monitor the convergence
    of the algorithm.
    """
    cdef:
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[INT, ndim=1] X_indices = X.indices
        np.ndarray[INT, ndim=1] X_indptr = X.indptr
        unsigned int n_samples = batch_slice.stop - batch_slice.start
        unsigned int n_clusters = centers.shape[0]
        unsigned int n_features = centers.shape[1]

        unsigned int batch_slice_start = batch_slice.start
        unsigned int batch_slice_stop = batch_slice.stop

        unsigned int sample_idx, center_idx, feature_idx
        unsigned int i, k
        int old_count, new_count
        DOUBLE inertia
        DOUBLE center_diff
        DOUBLE squared_diff = 0.0

        # TODO: reuse a array preallocated outside of the mini batch main loop
        np.ndarray[DOUBLE, ndim=1] new_center = np.zeros(
            n_features, dtype=np.double)
        np.ndarray[INT, ndim=1] nearest_center = np.zeros(
            n_samples, dtype=np.int32)

    # step 1: assign minibatch samples to there nearest center
    inertia = _assign_labels_csr(
        X, x_squared_norms, batch_slice_start, batch_slice_stop,
        centers, nearest_center)

    # step 2: move centers to the mean of both old and newly assigned samples
    for center_idx in range(n_clusters):
        old_count = counts[center_idx]
        new_count = old_count

        # rescale the old center to reflect it previous accumulated
        # weight w.r.t. the new data that will be incrementally contributed
        new_center[:] = centers[center_idx]
        new_center *= old_count

        # iterate of over samples assigned to this cluster to move the center
        # location by inplace summation
        for i in range(n_samples):
            if nearest_center[i] != center_idx:
                continue
            sample_idx = batch_slice_start + i
            new_count += 1

            # inplace sum with new samples that are members of this cluster
            # and update of the incremental squared difference update of the
            # center position
            for k in range(X_indptr[sample_idx], X_indptr[sample_idx + 1]):
                new_center[X_indices[k]] += X_data[k]

        # inplace rescale center with updated count
        if new_count > old_count:
            # update the count statistics for this center
            counts[center_idx] = new_count

            # re-scale the updated center with the total new counts
            new_center /= new_count

            # update the incremental computation of the squared total
            # centers position change
            for feature_idx in range(n_features):
                squared_diff += (new_center[feature_idx]
                                 - centers[center_idx, feature_idx]) ** 2

            # save the updated center position
            centers[center_idx] = new_center

    return inertia, squared_diff


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def csr_row_norm_l2(X, squared=True):
    """Get L2 norm of each row in CSR matrix X."""
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
