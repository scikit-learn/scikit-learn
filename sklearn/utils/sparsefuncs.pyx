# Author: Mathieu Blondel
#         Olivier Grisel
#
# License: BSD Style.

cimport numpy as np
import numpy as np
cimport cython

cdef extern from "math.h":
    double fabs(double f)
    double sqrt(double f)

ctypedef np.float64_t DOUBLE


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR matrix

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int ind
    cdef double diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[DOUBLE, ndim=1] means = np.asarray(X.mean(axis=0))[0]

    # variances[j] contains the variance of feature j
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    # counts[j] contains the number of samples where feature j is non-zero
    counts = np.zeros_like(means)

    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            ind = X_indices[j]
            diff = X_data[j] - means[ind]
            variances[ind] += diff * diff
            counts[ind] += 1

    nz = n_samples - counts
    variances += nz * means ** 2
    variances /= n_samples

    return means, variances


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_row_normalize_l1(X):
    """Inplace row normalize using the l1 norm"""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    # the column indices for row i are stored in:
    #    indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in:
    #    data[indptr[i]:indptr[i+1]]
    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += fabs(X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_row_normalize_l2(X):
    """Inplace row normalize using the l2 norm"""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        sum_ = sqrt(sum_)

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def inplace_csr_column_scale(X, np.ndarray[DOUBLE, ndim=1] scale):
    """Inplace column scaling of a CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef unsigned int i, j
    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] *= scale[X_indices[j]]
