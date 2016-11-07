# Authors: Mathieu Blondel
#          Olivier Grisel
#          Peter Prettenhofer
#          Lars Buitinck
#          Giorgio Patrini
#
# License: BSD 3 clause

#!python
#cython: boundscheck=False, wraparound=False, cdivision=True

from libc.math cimport fabs, sqrt, pow
cimport numpy as np
import numpy as np
import scipy.sparse as sp
cimport cython
from cython cimport floating

np.import_array()


ctypedef np.float64_t DOUBLE

def csr_row_norms(X):
    """L2 norm of each row in CSR matrix X."""
    if X.dtype != np.float32:
        X = X.astype(np.float64)
    return _csr_row_norms(X.data, X.shape, X.indices, X.indptr)


def _csr_row_norms(np.ndarray[floating, ndim=1, mode="c"] X_data,
                   shape,
                   np.ndarray[int, ndim=1, mode="c"] X_indices,
                   np.ndarray[int, ndim=1, mode="c"] X_indptr):
    cdef:
        unsigned int n_samples = shape[0]
        unsigned int n_features = shape[1]
        np.ndarray[DOUBLE, ndim=1, mode="c"] norms

        np.npy_intp i, j
        double sum_

    norms = np.zeros(n_samples, dtype=np.float64)

    for i in range(n_samples):
        sum_ = 0.0
        for j in range(X_indptr[i], X_indptr[i + 1]):
            sum_ += X_data[j] * X_data[j]
        norms[i] = sum_

    return norms


def csr_mean_variance_axis0(X):
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
    if X.dtype != np.float32:
        X = X.astype(np.float64)
    return _csr_mean_variance_axis0(X.data, X.shape, X.indices)


def _csr_mean_variance_axis0(np.ndarray[floating, ndim=1, mode="c"] X_data,
                             shape,
                             np.ndarray[int, ndim=1] X_indices):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef unsigned int n_samples = shape[0]
    cdef unsigned int n_features = shape[1]

    cdef unsigned int i
    cdef unsigned int non_zero = X_indices.shape[0]
    cdef unsigned int col_ind
    cdef floating diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[floating, ndim=1] means
    # variances[j] contains the variance of feature j
    cdef np.ndarray[floating, ndim=1] variances

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    means = np.zeros(n_features, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    # counts[j] contains the number of samples where feature j is non-zero
    cdef np.ndarray[int, ndim=1] counts = np.zeros(n_features,
                                                   dtype=np.int32)

    for i in xrange(non_zero):
        col_ind = X_indices[i]
        means[col_ind] += X_data[i]

    means /= n_samples

    for i in xrange(non_zero):
        col_ind = X_indices[i]
        diff = X_data[i] - means[col_ind]
        variances[col_ind] += diff * diff
        counts[col_ind] += 1

    for i in xrange(n_features):
        variances[i] += (n_samples - counts[i]) * means[i] ** 2
        variances[i] /= n_samples

    return means, variances


def csc_mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSC matrix

    Parameters
    ----------
    X: CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    if X.dtype != np.float32:
        X = X.astype(np.float64)
    return _csc_mean_variance_axis0(X.data, X.shape, X.indices, X.indptr)


def _csc_mean_variance_axis0(np.ndarray[floating, ndim=1] X_data,
                             shape,
                             np.ndarray[int, ndim=1] X_indices,
                             np.ndarray[int, ndim=1] X_indptr):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef unsigned int n_samples = shape[0]
    cdef unsigned int n_features = shape[1]

    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int counts
    cdef unsigned int startptr
    cdef unsigned int endptr
    cdef floating diff

    # means[j] contains the mean of feature j
    cdef np.ndarray[floating, ndim=1] means
    # variances[j] contains the variance of feature j
    cdef np.ndarray[floating, ndim=1] variances
    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    means = np.zeros(n_features, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    for i in xrange(n_features):

        startptr = X_indptr[i]
        endptr = X_indptr[i + 1]
        counts = endptr - startptr

        for j in xrange(startptr, endptr):
            means[i] += X_data[j]
        means[i] /= n_samples

        for j in xrange(startptr, endptr):
            diff = X_data[j] - means[i]
            variances[i] += diff * diff

        variances[i] += (n_samples - counts) * means[i] * means[i]
        variances[i] /= n_samples

    return means, variances


def incr_mean_variance_axis0(X, last_mean, last_var, unsigned long last_n):
    """Compute mean and variance along axis 0 on a CSR or CSC matrix.

    last_mean, last_var are the statistics computed at the last step by this
    function. Both must be initilized to 0.0. last_n is the
    number of samples encountered until now and is initialized at 0.

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
      Input data.

    last_mean: float array with shape (n_features,)
      Array of feature-wise means to update with the new data X.

    last_var: float array with shape (n_features,)
      Array of feature-wise var to update with the new data X.

    last_n: int
      Number of samples seen so far, before X.

    Returns
    -------

    updated_mean: float array with shape (n_features,)
      Feature-wise means

    updated_variance: float array with shape (n_features,)
      Feature-wise variances

    updated_n : int
      Updated number of samples seen

    References
    ----------

    T. Chan, G. Golub, R. LeVeque. Algorithms for computing the sample
      variance: recommendations, The American Statistician, Vol. 37, No. 3,
      pp. 242-247

    Also, see the non-sparse implementation of this in
    `utils.extmath._batch_mean_variance_update`.

    """
    if X.dtype != np.float32:
        X = X.astype(np.float64)
    return _incr_mean_variance_axis0(X.data, X.shape, X.indices, X.indptr,
                                     X.format, last_mean, last_var, last_n)


def _incr_mean_variance_axis0(np.ndarray[floating, ndim=1] X_data,
                              shape,
                              np.ndarray[int, ndim=1] X_indices,
                              np.ndarray[int, ndim=1] X_indptr,
                              X_format,
                              last_mean,
                              last_var,
                              unsigned long last_n):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef unsigned long n_samples = shape[0]
    cdef unsigned int n_features = shape[1]
    cdef unsigned int i

    # last = stats until now
    # new = the current increment
    # updated = the aggregated stats
    # when arrays, they are indexed by i per-feature
    cdef np.ndarray[floating, ndim=1] new_mean
    cdef np.ndarray[floating, ndim=1] new_var
    cdef np.ndarray[floating, ndim=1] updated_mean
    cdef np.ndarray[floating, ndim=1] updated_var
    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    new_mean = np.zeros(n_features, dtype=dtype)
    new_var = np.zeros_like(new_mean, dtype=dtype)
    updated_mean = np.zeros_like(new_mean, dtype=dtype)
    updated_var = np.zeros_like(new_mean, dtype=dtype)

    cdef unsigned long new_n
    cdef unsigned long updated_n
    cdef floating last_over_new_n

    # Obtain new stats first
    new_n = n_samples

    if X_format == 'csr':
        # X is a CSR matrix
        new_mean, new_var = _csr_mean_variance_axis0(X_data, shape, X_indices)
    else:
        # X is a CSC matrix
        new_mean, new_var = _csc_mean_variance_axis0(X_data, shape, X_indices,
                                                     X_indptr)

    # First pass
    if last_n == 0:
        return new_mean, new_var, new_n
    # Next passes
    else:
        updated_n = last_n + new_n
        last_over_new_n = last_n / new_n

    for i in xrange(n_features):
        # Unnormalized old stats
        last_mean[i] *= last_n
        last_var[i] *= last_n

        # Unnormalized new stats
        new_mean[i] *= new_n
        new_var[i] *= new_n

        # Update stats
        updated_var[i] = (last_var[i] + new_var[i] +
                          last_over_new_n / updated_n *
                          (last_mean[i] / last_over_new_n - new_mean[i]) ** 2)

        updated_mean[i] = (last_mean[i] + new_mean[i]) / updated_n
        updated_var[i] = updated_var[i] / updated_n

    return updated_mean, updated_var, updated_n


def inplace_csr_row_normalize_l1(X):
    """Inplace row normalize using the l1 norm"""
    _inplace_csr_row_normalize_l1(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l1(np.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  np.ndarray[int, ndim=1] X_indices,
                                  np.ndarray[int, ndim=1] X_indptr):
    cdef unsigned int n_samples = shape[0]
    cdef unsigned int n_features = shape[1]

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


def inplace_csr_row_normalize_l2(X):
    """Inplace row normalize using the l2 norm"""
    _inplace_csr_row_normalize_l2(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l2(np.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  np.ndarray[int, ndim=1] X_indices,
                                  np.ndarray[int, ndim=1] X_indptr):
    cdef unsigned int n_samples = shape[0]
    cdef unsigned int n_features = shape[1]

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


def assign_rows_csr(X,
                    np.ndarray[np.npy_intp, ndim=1] X_rows,
                    np.ndarray[np.npy_intp, ndim=1] out_rows,
                    np.ndarray[floating, ndim=2, mode="c"] out):
    """Densify selected rows of a CSR matrix into a preallocated array.

    Like out[out_rows] = X[X_rows].toarray() but without copying.
    No-copy supported for both dtype=np.float32 and dtype=np.float64.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix, shape=(n_samples, n_features)
    X_rows : array, dtype=np.intp, shape=n_rows
    out_rows : array, dtype=np.intp, shape=n_rows
    out : array, shape=(arbitrary, n_features)
    """
    cdef:
        # npy_intp (np.intp in Python) is what np.where returns,
        # but int is what scipy.sparse uses.
        int i, ind, j
        np.npy_intp rX
        np.ndarray[floating, ndim=1] data = X.data
        np.ndarray[int, ndim=1] indices = X.indices, indptr = X.indptr

    if X_rows.shape[0] != out_rows.shape[0]:
        raise ValueError("cannot assign %d rows to %d"
                         % (X_rows.shape[0], out_rows.shape[0]))

    out[out_rows] = 0.
    for i in range(X_rows.shape[0]):
        rX = X_rows[i]
        for ind in range(indptr[rX], indptr[rX + 1]):
            j = indices[ind]
            out[out_rows[i], j] = data[ind]
