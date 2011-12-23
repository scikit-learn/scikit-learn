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
cdef _mean_variance_axis0(unsigned int n_samples,
                          unsigned int n_features,
                          np.ndarray[DOUBLE, ndim=1] X_data,
                          np.ndarray[int, ndim=1] X_indices,
                          np.ndarray[int, ndim=1] X_indptr,
                          np.ndarray[DOUBLE, ndim=1] means,
                          np.ndarray[DOUBLE, ndim=1] variances):
    """Compute the mean and variance of a CSR aling axis=0

    This is the private Cython function with pre-allocated meant to be called
    by the public functions of this module.

    means[j] contains the mean of feature j (must be precomputed by the
    caller)

    variances[j] contains the variance of feature j (must be set to zero by
    the caller)
    """

    # the column indices for row i are stored in:
    #    indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in:
    #    data[indptr[i]:indptr[i+1]]
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int ptr
    cdef unsigned int ind
    cdef double diff

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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def mean_variance_axis0(X):
    """Compute mean and variance along axis 0

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    (means, variances):
        means: np.mean(X, axis=0)
        variances: np.var(X, axis=0)
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    # store the means in a 1d-array
    cdef np.ndarray[DOUBLE, ndim=1] means = np.asarray(X.mean(axis=0))[0]
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros_like(means)

    _mean_variance_axis0(n_samples, n_features, X_data, X_indices, X_indptr,
                         means, variances)
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
def inplace_csr_column_normalize_l2(X):
    """Inplace column normalize using the l2 norm

    This amounts to dividing each non-zero component by the variance of the
    feature assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    Return
    ------
    variances: float array with shape (n_features,)
    """
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[int, ndim=1] X_indices = X.indices
    cdef np.ndarray[int, ndim=1] X_indptr = X.indptr

    cdef np.ndarray[DOUBLE, ndim=1] means = np.asarray(X.mean(axis=0))[0]
    cdef np.ndarray[DOUBLE, ndim=1] variances = np.zeros(n_features)
    _mean_variance_axis0(n_samples, n_features, X_data, X_indices, X_indptr,
                         means, variances)

    cdef unsigned int i, j
    for i in xrange(n_samples):
        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            if variances[X_indices[j]] > 0.0:
                X_data[j] /= variances[X_indices[j]]

    return variances
