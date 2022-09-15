# Authors: Mathieu Blondel
#          Olivier Grisel
#          Peter Prettenhofer
#          Lars Buitinck
#          Giorgio Patrini
#
# License: BSD 3 clause

#!python

from libc.math cimport fabs, sqrt
cimport numpy as cnp
import numpy as np
from cython cimport floating
from numpy.math cimport isnan

cnp.import_array()

ctypedef fused integral:
    int
    long long

ctypedef cnp.float64_t DOUBLE


def csr_row_norms(X):
    """L2 norm of each row in CSR matrix X."""
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)
    return _csr_row_norms(X.data, X.shape, X.indices, X.indptr)


def _csr_row_norms(cnp.ndarray[floating, ndim=1, mode="c"] X_data,
                   shape,
                   cnp.ndarray[integral, ndim=1, mode="c"] X_indices,
                   cnp.ndarray[integral, ndim=1, mode="c"] X_indptr):
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


def csr_mean_variance_axis0(X, weights=None, return_sum_weights=False):
    """Compute mean and variance along axis 0 on a CSR matrix

    Uses a np.float64 accumulator.

    Parameters
    ----------
    X : CSR sparse matrix, shape (n_samples, n_features)
        Input data.

    weights : ndarray of shape (n_samples,), dtype=floating, default=None
        If it is set to None samples will be equally weighted.

        .. versionadded:: 0.24

    return_sum_weights : bool, default=False
        If True, returns the sum of weights seen for each feature.

        .. versionadded:: 0.24

    Returns
    -------
    means : float array with shape (n_features,)
        Feature-wise means

    variances : float array with shape (n_features,)
        Feature-wise variances

    sum_weights : ndarray of shape (n_features,), dtype=floating
        Returned if return_sum_weights is True.
    """
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)

    if weights is None:
        weights = np.ones(X.shape[0], dtype=X.dtype)

    means, variances, sum_weights = _csr_mean_variance_axis0(
        X.data, X.shape[0], X.shape[1], X.indices, X.indptr, weights)

    if return_sum_weights:
        return means, variances, sum_weights
    return means, variances


def _csr_mean_variance_axis0(cnp.ndarray[floating, ndim=1, mode="c"] X_data,
                             unsigned long long n_samples,
                             unsigned long long n_features,
                             cnp.ndarray[integral, ndim=1] X_indices,
                             cnp.ndarray[integral, ndim=1] X_indptr,
                             cnp.ndarray[floating, ndim=1] weights):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        cnp.npy_intp i
        unsigned long long row_ind
        integral col_ind
        cnp.float64_t diff
        # means[j] contains the mean of feature j
        cnp.ndarray[cnp.float64_t, ndim=1] means = np.zeros(n_features)
        # variances[j] contains the variance of feature j
        cnp.ndarray[cnp.float64_t, ndim=1] variances = np.zeros(n_features)

        cnp.ndarray[cnp.float64_t, ndim=1] sum_weights = np.full(
            fill_value=np.sum(weights, dtype=np.float64), shape=n_features)
        cnp.ndarray[cnp.float64_t, ndim=1] sum_weights_nz = np.zeros(
            shape=n_features)
        cnp.ndarray[cnp.float64_t, ndim=1] correction = np.zeros(
            shape=n_features)

        cnp.ndarray[cnp.uint64_t, ndim=1] counts = np.full(
            fill_value=weights.shape[0], shape=n_features, dtype=np.uint64)
        cnp.ndarray[cnp.uint64_t, ndim=1] counts_nz = np.zeros(
            shape=n_features, dtype=np.uint64)

    for row_ind in range(len(X_indptr) - 1):
        for i in range(X_indptr[row_ind], X_indptr[row_ind + 1]):
            col_ind = X_indices[i]
            if not isnan(X_data[i]):
                means[col_ind] += <cnp.float64_t>(X_data[i]) * weights[row_ind]
                # sum of weights where X[:, col_ind] is non-zero
                sum_weights_nz[col_ind] += weights[row_ind]
                # number of non-zero elements of X[:, col_ind]
                counts_nz[col_ind] += 1
            else:
                # sum of weights where X[:, col_ind] is not nan
                sum_weights[col_ind] -= weights[row_ind]
                # number of non nan elements of X[:, col_ind]
                counts[col_ind] -= 1

    for i in range(n_features):
        means[i] /= sum_weights[i]

    for row_ind in range(len(X_indptr) - 1):
        for i in range(X_indptr[row_ind], X_indptr[row_ind + 1]):
            col_ind = X_indices[i]
            if not isnan(X_data[i]):
                diff = X_data[i] - means[col_ind]
                # correction term of the corrected 2 pass algorithm.
                # See "Algorithms for computing the sample variance: analysis
                # and recommendations", by Chan, Golub, and LeVeque.
                correction[col_ind] += diff * weights[row_ind]
                variances[col_ind] += diff * diff * weights[row_ind]

    for i in range(n_features):
        if counts[i] != counts_nz[i]:
            correction[i] -= (sum_weights[i] - sum_weights_nz[i]) * means[i]
        correction[i] = correction[i]**2 / sum_weights[i]
        if counts[i] != counts_nz[i]:
            # only compute it when it's guaranteed to be non-zero to avoid
            # catastrophic cancellation.
            variances[i] += (sum_weights[i] - sum_weights_nz[i]) * means[i]**2
        variances[i] = (variances[i] - correction[i]) / sum_weights[i]

    if floating is float:
        return (np.array(means, dtype=np.float32),
                np.array(variances, dtype=np.float32),
                np.array(sum_weights, dtype=np.float32))
    else:
        return means, variances, sum_weights


def csc_mean_variance_axis0(X, weights=None, return_sum_weights=False):
    """Compute mean and variance along axis 0 on a CSC matrix

    Uses a np.float64 accumulator.

    Parameters
    ----------
    X : CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    weights : ndarray of shape (n_samples,), dtype=floating, default=None
        If it is set to None samples will be equally weighted.

        .. versionadded:: 0.24

    return_sum_weights : bool, default=False
        If True, returns the sum of weights seen for each feature.

        .. versionadded:: 0.24

    Returns
    -------
    means : float array with shape (n_features,)
        Feature-wise means

    variances : float array with shape (n_features,)
        Feature-wise variances

    sum_weights : ndarray of shape (n_features,), dtype=floating
        Returned if return_sum_weights is True.
    """
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)

    if weights is None:
        weights = np.ones(X.shape[0], dtype=X.dtype)

    means, variances, sum_weights = _csc_mean_variance_axis0(
        X.data, X.shape[0], X.shape[1], X.indices, X.indptr, weights)

    if return_sum_weights:
        return means, variances, sum_weights
    return means, variances


def _csc_mean_variance_axis0(cnp.ndarray[floating, ndim=1, mode="c"] X_data,
                             unsigned long long n_samples,
                             unsigned long long n_features,
                             cnp.ndarray[integral, ndim=1] X_indices,
                             cnp.ndarray[integral, ndim=1] X_indptr,
                             cnp.ndarray[floating, ndim=1] weights):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        cnp.npy_intp i
        unsigned long long col_ind
        integral row_ind
        cnp.float64_t diff
        # means[j] contains the mean of feature j
        cnp.ndarray[cnp.float64_t, ndim=1] means = np.zeros(n_features)
        # variances[j] contains the variance of feature j
        cnp.ndarray[cnp.float64_t, ndim=1] variances = np.zeros(n_features)

        cnp.ndarray[cnp.float64_t, ndim=1] sum_weights = np.full(
            fill_value=np.sum(weights, dtype=np.float64), shape=n_features)
        cnp.ndarray[cnp.float64_t, ndim=1] sum_weights_nz = np.zeros(
            shape=n_features)
        cnp.ndarray[cnp.float64_t, ndim=1] correction = np.zeros(
            shape=n_features)

        cnp.ndarray[cnp.uint64_t, ndim=1] counts = np.full(
            fill_value=weights.shape[0], shape=n_features, dtype=np.uint64)
        cnp.ndarray[cnp.uint64_t, ndim=1] counts_nz = np.zeros(
            shape=n_features, dtype=np.uint64)

    for col_ind in range(n_features):
        for i in range(X_indptr[col_ind], X_indptr[col_ind + 1]):
            row_ind = X_indices[i]
            if not isnan(X_data[i]):
                means[col_ind] += <cnp.float64_t>(X_data[i]) * weights[row_ind]
                # sum of weights where X[:, col_ind] is non-zero
                sum_weights_nz[col_ind] += weights[row_ind]
                # number of non-zero elements of X[:, col_ind]
                counts_nz[col_ind] += 1
            else:
                # sum of weights where X[:, col_ind] is not nan
                sum_weights[col_ind] -= weights[row_ind]
                # number of non nan elements of X[:, col_ind]
                counts[col_ind] -= 1

    for i in range(n_features):
        means[i] /= sum_weights[i]

    for col_ind in range(n_features):
        for i in range(X_indptr[col_ind], X_indptr[col_ind + 1]):
            row_ind = X_indices[i]
            if not isnan(X_data[i]):
                diff = X_data[i] - means[col_ind]
                # correction term of the corrected 2 pass algorithm.
                # See "Algorithms for computing the sample variance: analysis
                # and recommendations", by Chan, Golub, and LeVeque.
                correction[col_ind] += diff * weights[row_ind]
                variances[col_ind] += diff * diff * weights[row_ind]

    for i in range(n_features):
        if counts[i] != counts_nz[i]:
            correction[i] -= (sum_weights[i] - sum_weights_nz[i]) * means[i]
        correction[i] = correction[i]**2 / sum_weights[i]
        if counts[i] != counts_nz[i]:
            # only compute it when it's guaranteed to be non-zero to avoid
            # catastrophic cancellation.
            variances[i] += (sum_weights[i] - sum_weights_nz[i]) * means[i]**2
        variances[i] = (variances[i] - correction[i]) / sum_weights[i]

    if floating is float:
        return (np.array(means, dtype=np.float32),
                np.array(variances, dtype=np.float32),
                np.array(sum_weights, dtype=np.float32))
    else:
        return means, variances, sum_weights


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

    return _incr_mean_variance_axis0(X.data,
                                     np.sum(weights),
                                     X.shape[1],
                                     X.indices,
                                     X.indptr,
                                     X.format,
                                     last_mean.astype(X_dtype, copy=False),
                                     last_var.astype(X_dtype, copy=False),
                                     last_n.astype(X_dtype, copy=False),
                                     weights.astype(X_dtype, copy=False))


def _incr_mean_variance_axis0(cnp.ndarray[floating, ndim=1] X_data,
                              floating n_samples,
                              unsigned long long n_features,
                              cnp.ndarray[int, ndim=1] X_indices,
                              # X_indptr might be either in32 or int64
                              cnp.ndarray[integral, ndim=1] X_indptr,
                              str X_format,
                              cnp.ndarray[floating, ndim=1] last_mean,
                              cnp.ndarray[floating, ndim=1] last_var,
                              cnp.ndarray[floating, ndim=1] last_n,
                              # previous sum of the weights (ie float)
                              cnp.ndarray[floating, ndim=1] weights):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        cnp.npy_intp i

    # last = stats until now
    # new = the current increment
    # updated = the aggregated stats
    # when arrays, they are indexed by i per-feature
    cdef:
        cnp.ndarray[floating, ndim=1] new_mean
        cnp.ndarray[floating, ndim=1] new_var
        cnp.ndarray[floating, ndim=1] updated_mean
        cnp.ndarray[floating, ndim=1] updated_var

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    new_mean = np.zeros(n_features, dtype=dtype)
    new_var = np.zeros_like(new_mean, dtype=dtype)
    updated_mean = np.zeros_like(new_mean, dtype=dtype)
    updated_var = np.zeros_like(new_mean, dtype=dtype)

    cdef:
        cnp.ndarray[floating, ndim=1] new_n
        cnp.ndarray[floating, ndim=1] updated_n
        cnp.ndarray[floating, ndim=1] last_over_new_n

    # Obtain new stats first
    updated_n = np.zeros(shape=n_features, dtype=dtype)
    last_over_new_n = np.zeros_like(updated_n, dtype=dtype)

    # X can be a CSR or CSC matrix
    if X_format == 'csr':
        new_mean, new_var, new_n = _csr_mean_variance_axis0(
            X_data, n_samples, n_features, X_indices, X_indptr, weights)
    else:  # X_format == 'csc'
        new_mean, new_var, new_n = _csc_mean_variance_axis0(
            X_data, n_samples, n_features, X_indices, X_indptr, weights)

    # First pass
    cdef bint is_first_pass = True
    for i in range(n_features):
        if last_n[i] > 0:
            is_first_pass = False
            break

    if is_first_pass:
        return new_mean, new_var, new_n

    for i in range(n_features):
        updated_n[i] = last_n[i] + new_n[i]

    # Next passes
    for i in range(n_features):
        if new_n[i] > 0:
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


def _inplace_csr_row_normalize_l1(cnp.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  cnp.ndarray[integral, ndim=1] X_indices,
                                  cnp.ndarray[integral, ndim=1] X_indptr):
    cdef unsigned long long n_samples = shape[0]
    cdef unsigned long long n_features = shape[1]

    # the column indices for row i are stored in:
    #    indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in:
    #    data[indptr[i]:indptr[i+1]]
    cdef cnp.npy_intp i, j
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


def _inplace_csr_row_normalize_l2(cnp.ndarray[floating, ndim=1] X_data,
                                  shape,
                                  cnp.ndarray[integral, ndim=1] X_indices,
                                  cnp.ndarray[integral, ndim=1] X_indptr):
    cdef integral n_samples = shape[0]
    cdef integral n_features = shape[1]

    cdef cnp.npy_intp i, j
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
                    cnp.ndarray[cnp.npy_intp, ndim=1] X_rows,
                    cnp.ndarray[cnp.npy_intp, ndim=1] out_rows,
                    cnp.ndarray[floating, ndim=2, mode="c"] out):
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
        cnp.npy_intp rX
        cnp.ndarray[floating, ndim=1] data = X.data
        cnp.ndarray[int, ndim=1] indices = X.indices, indptr = X.indptr

    if X_rows.shape[0] != out_rows.shape[0]:
        raise ValueError("cannot assign %d rows to %d"
                         % (X_rows.shape[0], out_rows.shape[0]))

    out[out_rows] = 0.
    for i in range(X_rows.shape[0]):
        rX = X_rows[i]
        for ind in range(indptr[rX], indptr[rX + 1]):
            j = indices[ind]
            out[out_rows[i], j] = data[ind]
