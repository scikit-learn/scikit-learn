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
    cdef DOUBLE *X_data_ptr = <DOUBLE *>X_data.data
    cdef INT *X_indptr_ptr = <INT *>X_indptr.data
    cdef INT *X_indices_ptr = <INT *>X_indices.data

    cdef int c_stride = centers.strides[0] / centers.strides[1]
    cdef DOUBLE *center_data_ptr = <DOUBLE *> centers.data

    cdef int n_samples = batch_slice.stop - batch_slice.start
    cdef int n_clusters = centers.shape[0]
    cdef int n_features = centers.shape[1]

    cdef int i = 0
    cdef DOUBLE center_squared_norm = 0.0
    cdef int batch_slice_start = batch_slice.start
    cdef int sample_idx = -1
    cdef int offset = -1
    cdef int xnnz = -1
    cdef int c = -1
    cdef double eta = 0.0, u = 0.0

    # TODO: reuse a arrays preallocated outside of the mini batch main loop
    cdef np.ndarray[DOUBLE, ndim=1] center_scales = np.ones(n_clusters,
                                                            dtype=np.float64)
    cdef DOUBLE *center_scales_ptr = <DOUBLE *> center_scales.data
    cdef np.ndarray[INT, ndim=1] nearest_center = np.zeros(n_samples,
                                                           dtype=np.int32)

    # step 1: assign minibatch samples to there nearest center
    _assign_labels(X_data, X_indices, X_indptr, x_squared_norms, batch_slice,
                   centers, nearest_center)

    # step 2: move centers to mean of old and newly assigned samples
    for i in range(n_samples):
        sample_idx = batch_slice_start + i  #batch[i]
        offset = X_indptr_ptr[sample_idx]
        xnnz = X_indptr_ptr[sample_idx + 1] - offset
        c = nearest_center[i]
        counts[c] = counts[c] + 1
        eta = 1.0 / counts[c]
        center_scales_ptr[c] *= (1.0 - eta)

        # check for scale underflow
        if center_scales_ptr[c] < 1e-9:
            scale(center_data_ptr + (c_stride * c), n_features,
                  center_scales_ptr[c])
            center_scales_ptr[c] = 1.0

        add(center_data_ptr + (c_stride * c), center_scales_ptr[c],
            X_data_ptr, X_indices_ptr, offset, xnnz, eta)

    # finally scale by scaling factors.
    for c from 0 <= c < n_clusters:
        scale(center_data_ptr + (c_stride * c), n_features,
              center_scales_ptr[c])


cdef void scale(DOUBLE *dense_vec, int n_features, double c):
    """Scale dense vector by constant c."""
    cdef int j
    for j from 0 <= j < n_features:
        dense_vec[j] *= c


@cython.cdivision(True)
cdef double add(DOUBLE *center_data_ptr, DOUBLE center_scale,
                DOUBLE *X_data_ptr, INT *X_indices_ptr, int offset,
                int xnnz, double c):
    """Scales example x by constant c and adds it to the weight vector w"""
    cdef int j
    cdef int idx
    cdef double val
    cdef double innerprod = 0.0
    cdef double xsqnorm = 0.0
    for j from 0 <= j < xnnz:
        idx = X_indices_ptr[offset + j]
        val = X_data_ptr[offset + j]
        innerprod += (center_data_ptr[idx] * val)
        xsqnorm += (val * val)
        center_data_ptr[idx] += val * c / center_scale
    return (xsqnorm * c * c) + (2.0 * innerprod * c * center_scale)


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
