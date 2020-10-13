# Authors: Mathieu Blondel
#          Olivier Grisel
#          Peter Prettenhofer
#          Lars Buitinck
#          Giorgio Patrini
#
# License: BSD 3 clause

#!python
# cython: boundscheck=False, wraparound=False, cdivision=True

from libc.math cimport fabs, sqrt, pow
cimport numpy as np
import numpy as np
import scipy.sparse as sp
cimport cython
from cython cimport floating
from numpy.math cimport isnan

np.import_array()

ctypedef fused integral:
    int
    long long

ctypedef np.float64_t DOUBLE


def csr_row_norms(X):
    """L2 norm of each row in CSR matrix X."""
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)
    return _csr_row_norms(X.data, X.shape, X.indices, X.indptr)


def _csr_row_norms(np.ndarray[floating, ndim=1, mode="c"] X_data,
                   shape,
                   np.ndarray[integral, ndim=1, mode="c"] X_indices,
                   np.ndarray[integral, ndim=1, mode="c"] X_indptr):
    cdef:
        unsigned long long n_samples = shape[0]
        unsigned long long i
        integral j
        double sum_

    norms = np.empty(n_samples, dtype=X_data.dtype)
    cdef floating[::1] norms_view = norms

    for i in range(n_samples):
        sum_ = 0.0
        for j in range(X_indptr[i], X_indptr[i + 1]):
            sum_ += X_data[j] * X_data[j]
        norms_view[i] = sum_

    return norms


def csr_mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR matrix

    Parameters
    ----------
    X : CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------
    means : float array with shape (n_features,)
        Feature-wise means

    variances : float array with shape (n_features,)
        Feature-wise variances

    """
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)
    last_mean = np.zeros(X.shape[1], dtype=X.dtype)
    last_var = np.zeros(X.shape[1], dtype=X.dtype)
    last_n = np.zeros(X.shape[1], dtype=X.dtype)
    weights = np.ones(X.shape[0], dtype=X.dtype)

    means, variances, _ = incr_mean_variance_axis0(
        X, last_mean, last_var, last_n, weights)

    return means, variances


def _csr_mean_variance_axis0(np.ndarray[floating, ndim=1, mode="c"] X_data,
                             unsigned long long n_samples,
                             unsigned long long n_features,
                             np.ndarray[integral, ndim=1] X_indices,
                             np.ndarray[integral, ndim=1] X_row_indices,
                             np.ndarray[floating, ndim=1] weights):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        np.npy_intp i
        unsigned long long non_zero = X_indices.shape[0]
        np.npy_intp col_ind
        floating diff
        # means[j] contains the mean of feature j
        np.ndarray[floating, ndim=1] means
        # variances[j] contains the variance of feature j
        np.ndarray[floating, ndim=1] variances

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    means = np.zeros(n_features, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)

    cdef:
        np.ndarray[np.int64_t, ndim=1] counts = np.zeros(n_features,
                                                         dtype=np.int64)
        np.ndarray[np.float64_t, ndim=1] counts_nan = np.zeros(
            n_features, dtype=np.float64)

    for i in range(non_zero):
        col_ind = X_indices[i]
        row_ind = X_row_indices[i]
        if not isnan(X_data[i]):
            means[col_ind] += (X_data[i] * weights[row_ind])
        else:
            counts_nan[col_ind] += weights[row_ind]

    for i in range(n_features):
        means[i] /= (n_samples - counts_nan[i])

    for i in range(non_zero):
        col_ind = X_indices[i]
        row_ind = X_row_indices[i]
        if not isnan(X_data[i]):
            diff = X_data[i] - means[col_ind]
            variances[col_ind] += diff * diff * weights[row_ind]
            counts[col_ind] += weights[row_ind]

    for i in range(n_features):
        variances[i] += (n_samples - counts_nan[i] - counts[i]) * means[i]**2
        variances[i] /= (n_samples - counts_nan[i])

    return means, variances, counts_nan.astype(dtype, copy=False)


def csc_mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSC matrix

    Parameters
    ----------
    X : CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------
    means : float array with shape (n_features,)
        Feature-wise means

    variances : float array with shape (n_features,)
        Feature-wise variances

    """
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)

    last_mean = np.zeros(X.shape[1], dtype=X.dtype)
    last_var = np.zeros(X.shape[1], dtype=X.dtype)
    last_n = np.zeros(X.shape[1], dtype=X.dtype)
    weights = np.ones(X.shape[0], dtype=X.dtype)
    means, variances, _ = incr_mean_variance_axis0(
        X, last_mean, last_var, last_n, weights)
    return means, variances


def incr_mean_variance_axis0(X, last_mean, last_var, last_n, weights=None):
    """Compute mean and variance along axis 0 on a CSR or CSC matrix.

    last_mean, last_var are the statistics computed at the last step by this
    function. Both must be initialized to 0.0. last_n is the
    number of samples encountered until now and is initialized at 0.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape (n_samples, n_features)
      Input data.

    last_mean : float array with shape (n_features,)
      Array of feature-wise means to update with the new data X.

    last_var : float array with shape (n_features,)
      Array of feature-wise var to update with the new data X.

    last_n : float array with shape (n_features,)
      Sum of the weights seen so far (if weights are all set to 1
      this will be the same as number of samples seen so far, before X).

    weights : float array with shape (n_samples,) or None. If it is set
      to None samples will be equally weighted.

    Returns
    -------
    updated_mean : float array with shape (n_features,)
      Feature-wise means

    updated_variance : float array with shape (n_features,)
      Feature-wise variances

    updated_n : int array with shape (n_features,)
      Updated number of samples seen

    Notes
    -----
    NaNs are ignored during the computation.

    References
    ----------
    T. Chan, G. Golub, R. LeVeque. Algorithms for computing the sample
      variance: recommendations, The American Statistician, Vol. 37, No. 3,
      pp. 242-247

    Also, see the non-sparse implementation of this in
    `utils.extmath._batch_mean_variance_update`.

    """
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)
    X_dtype = X.dtype
    if weights is None:
       weights = np.ones(X.shape[0], dtype=X_dtype)
    elif weights.dtype not in [np.float32, np.float64]:
        weights = weights.astype(np.float64, copy=False)
    if last_n.dtype not in [np.float32, np.float64]:
        last_n = last_n.astype(np.float64, copy=False)

    ind_rows, ind_cols, X_data = sp.find(X)

    return _incr_mean_variance_axis0(X_data,
                                     np.sum(weights),
                                     X.shape[1],
                                     ind_rows,
                                     ind_cols,
                                     X.indptr,
                                     X.format,
                                     last_mean.astype(X_dtype, copy=False),
                                     last_var.astype(X_dtype, copy=False),
                                     last_n.astype(X_dtype, copy=False),
                                     weights.astype(X_dtype, copy=False))


def _incr_mean_variance_axis0(np.ndarray[floating, ndim=1] X_data,
                              floating n_samples,
                              unsigned long long n_features,
                              # find() returns int32 so we can set them to long
                              np.ndarray[int, ndim=1] X_row_ind,
                              np.ndarray[int, ndim=1] X_indices,
                              # X_indptr might be either in32 or int64
                              np.ndarray[integral, ndim=1] X_indptr,
                              str X_format,
                              np.ndarray[floating, ndim=1] last_mean,
                              np.ndarray[floating, ndim=1] last_var,
                              np.ndarray[floating, ndim=1] last_n,
                              # previous sum of the weights (ie float)
                              np.ndarray[floating, ndim=1] weights):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        np.npy_intp i

    # last = stats until now
    # new = the current increment
    # updated = the aggregated stats
    # when arrays, they are indexed by i per-feature
    cdef:
        np.ndarray[floating, ndim=1] new_mean
        np.ndarray[floating, ndim=1] new_var
        np.ndarray[floating, ndim=1] updated_mean
        np.ndarray[floating, ndim=1] updated_var

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    new_mean = np.zeros(n_features, dtype=dtype)
    new_var = np.zeros_like(new_mean, dtype=dtype)
    updated_mean = np.zeros_like(new_mean, dtype=dtype)
    updated_var = np.zeros_like(new_mean, dtype=dtype)

    cdef:
        np.ndarray[floating, ndim=1] new_n
        np.ndarray[floating, ndim=1] updated_n
        np.ndarray[floating, ndim=1] last_over_new_n
        np.ndarray[floating, ndim=1] counts_nan

    # Obtain new stats first
    new_n = np.full(n_features, n_samples, dtype=dtype)
    updated_n = np.zeros_like(new_n, dtype=dtype)
    last_over_new_n = np.zeros_like(new_n, dtype=dtype)

    # X can be a CSR or CSC matrix
    new_mean, new_var, counts_nan = _csr_mean_variance_axis0(
        X_data, n_samples, n_features, X_indices, X_row_ind, weights)

    for i in range(n_features):
        new_n[i] -= counts_nan[i]

    # First pass
    cdef bint is_first_pass = True
    for i in range(n_features):
        if last_n[i] > 0:
            is_first_pass = False
            break
    if is_first_pass:
        return new_mean, new_var, new_n

    # Next passes
    for i in range(n_features):
        if new_n[i] > 0:
            updated_n[i] = last_n[i] + new_n[i]
            last_over_new_n[i] = dtype(last_n[i]) / dtype(new_n[i])
            # Unnormalized stats
            last_mean[i] *= last_n[i]
            last_var[i] *= last_n[i]
            new_mean[i] *= new_n[i]
            new_var[i] *= new_n[i]
            # Update stats
            updated_var[i] = (
                last_var[i] + new_var[i] +
                last_over_new_n[i] / updated_n[i] *
                (last_mean[i] / last_over_new_n[i] - new_mean[i])**2
            )
            updated_mean[i] = (last_mean[i] + new_mean[i]) / updated_n[i]
            updated_var[i] /= updated_n[i]
        else:
            updated_var[i] = last_var[i]
            updated_mean[i] = last_mean[i]
            updated_n[i] = last_n[i]

    return updated_mean, updated_var, updated_n


def inplace_csr_row_normalize_l1(X):
    """Inplace row normalize using the l1 norm"""
    _inplace_csr_row_normalize_l1(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l1(np.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  np.ndarray[integral, ndim=1] X_indices,
                                  np.ndarray[integral, ndim=1] X_indptr):
    cdef unsigned long long n_samples = shape[0]
    cdef unsigned long long n_features = shape[1]

    # the column indices for row i are stored in:
    #    indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in:
    #    data[indptr[i]:indptr[i+1]]
    cdef np.npy_intp i, j
    cdef double sum_

    for i in range(n_samples):
        sum_ = 0.0

        for j in range(X_indptr[i], X_indptr[i + 1]):
            sum_ += fabs(X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        for j in range(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_


def inplace_csr_row_normalize_l2(X):
    """Inplace row normalize using the l2 norm"""
    _inplace_csr_row_normalize_l2(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l2(np.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  np.ndarray[integral, ndim=1] X_indices,
                                  np.ndarray[integral, ndim=1] X_indptr):
    cdef integral n_samples = shape[0]
    cdef integral n_features = shape[1]

    cdef np.npy_intp i, j
    cdef double sum_

    for i in range(n_samples):
        sum_ = 0.0

        for j in range(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        if sum_ == 0.0:
            # do not normalize empty rows (can happen if CSR is not pruned
            # correctly)
            continue

        sum_ = sqrt(sum_)

        for j in range(X_indptr[i], X_indptr[i + 1]):
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
