"""Utilities to work with sparse matrices and arrays written in Cython."""

# Authors: Mathieu Blondel
#          Olivier Grisel
#          Peter Prettenhofer
#          Lars Buitinck
#          Giorgio Patrini
#
# License: BSD 3 clause

from libc.math cimport fabs, sqrt, isnan
from libc.stdint cimport intptr_t

import numpy as np
from cython cimport floating
from ..utils._typedefs cimport float64_t, int32_t, int64_t, intp_t, uint64_t


ctypedef fused integral:
    int32_t
    int64_t


def csr_row_norms(X):
    """Squared L2 norm of each row in CSR matrix X."""
    if X.dtype not in [np.float32, np.float64]:
        X = X.astype(np.float64)
    return _sqeuclidean_row_norms_sparse(X.data, X.indptr)


def _sqeuclidean_row_norms_sparse(
    const floating[::1] X_data,
    const integral[::1] X_indptr,
):
    cdef:
        integral n_samples = X_indptr.shape[0] - 1
        integral i, j

    dtype = np.float32 if floating is float else np.float64

    cdef floating[::1] squared_row_norms = np.zeros(n_samples, dtype=dtype)

    with nogil:
        for i in range(n_samples):
            for j in range(X_indptr[i], X_indptr[i + 1]):
                squared_row_norms[i] += X_data[j] * X_data[j]

    return np.asarray(squared_row_norms)


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


def _csr_mean_variance_axis0(
    const floating[::1] X_data,
    uint64_t n_samples,
    uint64_t n_features,
    const integral[:] X_indices,
    const integral[:] X_indptr,
    const floating[:] weights,
):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        intp_t row_ind
        uint64_t feature_idx
        integral i, col_ind
        float64_t diff
        # means[j] contains the mean of feature j
        float64_t[::1] means = np.zeros(n_features)
        # variances[j] contains the variance of feature j
        float64_t[::1] variances = np.zeros(n_features)

        float64_t[::1] sum_weights = np.full(
            fill_value=np.sum(weights, dtype=np.float64), shape=n_features
        )
        float64_t[::1] sum_weights_nz = np.zeros(shape=n_features)
        float64_t[::1] correction = np.zeros(shape=n_features)

        uint64_t[::1] counts = np.full(
            fill_value=weights.shape[0], shape=n_features, dtype=np.uint64
        )
        uint64_t[::1] counts_nz = np.zeros(shape=n_features, dtype=np.uint64)

    for row_ind in range(len(X_indptr) - 1):
        for i in range(X_indptr[row_ind], X_indptr[row_ind + 1]):
            col_ind = X_indices[i]
            if not isnan(X_data[i]):
                means[col_ind] += <float64_t>(X_data[i]) * weights[row_ind]
                # sum of weights where X[:, col_ind] is non-zero
                sum_weights_nz[col_ind] += weights[row_ind]
                # number of non-zero elements of X[:, col_ind]
                counts_nz[col_ind] += 1
            else:
                # sum of weights where X[:, col_ind] is not nan
                sum_weights[col_ind] -= weights[row_ind]
                # number of non nan elements of X[:, col_ind]
                counts[col_ind] -= 1

    for feature_idx in range(n_features):
        means[feature_idx] /= sum_weights[feature_idx]

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

    for feature_idx in range(n_features):
        if counts[feature_idx] != counts_nz[feature_idx]:
            correction[feature_idx] -= (
                sum_weights[feature_idx] - sum_weights_nz[feature_idx]
            ) * means[feature_idx]
        correction[feature_idx] = correction[feature_idx]**2 / sum_weights[feature_idx]
        if counts[feature_idx] != counts_nz[feature_idx]:
            # only compute it when it's guaranteed to be non-zero to avoid
            # catastrophic cancellation.
            variances[feature_idx] += (
                sum_weights[feature_idx] - sum_weights_nz[feature_idx]
            ) * means[feature_idx]**2
        variances[feature_idx] = (
            (variances[feature_idx] - correction[feature_idx]) /
            sum_weights[feature_idx]
        )

    if floating is float:
        return (
            np.array(means, dtype=np.float32),
            np.array(variances, dtype=np.float32),
            np.array(sum_weights, dtype=np.float32),
        )
    else:
        return (
            np.asarray(means), np.asarray(variances), np.asarray(sum_weights)
        )


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


def _csc_mean_variance_axis0(
    const floating[::1] X_data,
    uint64_t n_samples,
    uint64_t n_features,
    const integral[:] X_indices,
    const integral[:] X_indptr,
    const floating[:] weights,
):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        integral i, row_ind
        uint64_t feature_idx, col_ind
        float64_t diff
        # means[j] contains the mean of feature j
        float64_t[::1] means = np.zeros(n_features)
        # variances[j] contains the variance of feature j
        float64_t[::1] variances = np.zeros(n_features)

        float64_t[::1] sum_weights = np.full(
            fill_value=np.sum(weights, dtype=np.float64), shape=n_features
        )
        float64_t[::1] sum_weights_nz = np.zeros(shape=n_features)
        float64_t[::1] correction = np.zeros(shape=n_features)

        uint64_t[::1] counts = np.full(
            fill_value=weights.shape[0], shape=n_features, dtype=np.uint64
        )
        uint64_t[::1] counts_nz = np.zeros(shape=n_features, dtype=np.uint64)

    for col_ind in range(n_features):
        for i in range(X_indptr[col_ind], X_indptr[col_ind + 1]):
            row_ind = X_indices[i]
            if not isnan(X_data[i]):
                means[col_ind] += <float64_t>(X_data[i]) * weights[row_ind]
                # sum of weights where X[:, col_ind] is non-zero
                sum_weights_nz[col_ind] += weights[row_ind]
                # number of non-zero elements of X[:, col_ind]
                counts_nz[col_ind] += 1
            else:
                # sum of weights where X[:, col_ind] is not nan
                sum_weights[col_ind] -= weights[row_ind]
                # number of non nan elements of X[:, col_ind]
                counts[col_ind] -= 1

    for feature_idx in range(n_features):
        means[feature_idx] /= sum_weights[feature_idx]

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

    for feature_idx in range(n_features):
        if counts[feature_idx] != counts_nz[feature_idx]:
            correction[feature_idx] -= (
                sum_weights[feature_idx] - sum_weights_nz[feature_idx]
            ) * means[feature_idx]
        correction[feature_idx] = correction[feature_idx]**2 / sum_weights[feature_idx]
        if counts[feature_idx] != counts_nz[feature_idx]:
            # only compute it when it's guaranteed to be non-zero to avoid
            # catastrophic cancellation.
            variances[feature_idx] += (
                sum_weights[feature_idx] - sum_weights_nz[feature_idx]
            ) * means[feature_idx]**2
        variances[feature_idx] = (
            (variances[feature_idx] - correction[feature_idx])
        ) / sum_weights[feature_idx]

    if floating is float:
        return (np.array(means, dtype=np.float32),
                np.array(variances, dtype=np.float32),
                np.array(sum_weights, dtype=np.float32))
    else:
        return (
            np.asarray(means), np.asarray(variances), np.asarray(sum_weights)
        )


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


def _incr_mean_variance_axis0(
    const floating[:] X_data,
    floating n_samples,
    uint64_t n_features,
    const int[:] X_indices,
    # X_indptr might be either int32 or int64
    const integral[:] X_indptr,
    str X_format,
    floating[:] last_mean,
    floating[:] last_var,
    floating[:] last_n,
    # previous sum of the weights (ie float)
    const floating[:] weights,
):
    # Implement the function here since variables using fused types
    # cannot be declared directly and can only be passed as function arguments
    cdef:
        uint64_t i

        # last = stats until now
        # new = the current increment
        # updated = the aggregated stats
        # when arrays, they are indexed by i per-feature
        floating[::1] new_mean
        floating[::1] new_var
        floating[::1] updated_mean
        floating[::1] updated_var

    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    new_mean = np.zeros(n_features, dtype=dtype)
    new_var = np.zeros_like(new_mean, dtype=dtype)
    updated_mean = np.zeros_like(new_mean, dtype=dtype)
    updated_var = np.zeros_like(new_mean, dtype=dtype)

    cdef:
        floating[::1] new_n
        floating[::1] updated_n
        floating[::1] last_over_new_n

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
        return np.asarray(new_mean), np.asarray(new_var), np.asarray(new_n)

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

    return (
        np.asarray(updated_mean),
        np.asarray(updated_var),
        np.asarray(updated_n),
    )


def inplace_csr_row_normalize_l1(X):
    """Normalize inplace the rows of a CSR matrix or array by their L1 norm.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix and scipy.sparse.csr_array, \
            shape=(n_samples, n_features)
        The input matrix or array to be modified inplace.

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from sklearn.utils.sparsefuncs_fast import inplace_csr_row_normalize_l1
    >>> X = csr_matrix(([1.0, 2.0, 3.0], [0, 2, 3], [0, 3, 4]), shape=(3, 4))
    >>> X.toarray()
    array([[1., 2., 0., 0.],
           [0., 0., 3., 0.],
           [0., 0., 0., 4.]])
    >>> inplace_csr_row_normalize_l1(X)
    >>> X.toarray()
    array([[0.33...   , 0.66...   , 0.        , 0.        ],
           [0.        , 0.        , 1.        , 0.        ],
           [0.        , 0.        , 0.        , 1.        ]])
    """
    _inplace_csr_row_normalize_l1(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l1(
    floating[:] X_data,
    shape,
    const integral[:] X_indices,
    const integral[:] X_indptr,
):
    cdef:
        uint64_t n_samples = shape[0]

        # the column indices for row i are stored in:
        #    indices[indptr[i]:indices[i+1]]
        # and their corresponding values are stored in:
        #    data[indptr[i]:indptr[i+1]]
        uint64_t i
        integral j
        double sum_

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
    """Normalize inplace the rows of a CSR matrix or array by their L2 norm.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix, shape=(n_samples, n_features)
        The input matrix or array to be modified inplace.

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from sklearn.utils.sparsefuncs_fast import inplace_csr_row_normalize_l2
    >>> X = csr_matrix(([1.0, 2.0, 3.0], [0, 2, 3], [0, 3, 4]), shape=(3, 4))
    >>> X.toarray()
    array([[1., 2., 0., 0.],
           [0., 0., 3., 0.],
           [0., 0., 0., 4.]])
    >>> inplace_csr_row_normalize_l2(X)
    >>> X.toarray()
    array([[0.44...   , 0.89...   , 0.        , 0.        ],
           [0.        , 0.        , 1.        , 0.        ],
           [0.        , 0.        , 0.        , 1.        ]])
    """
    _inplace_csr_row_normalize_l2(X.data, X.shape, X.indices, X.indptr)


def _inplace_csr_row_normalize_l2(
    floating[:] X_data,
    shape,
    const integral[:] X_indices,
    const integral[:] X_indptr,
):
    cdef:
        uint64_t n_samples = shape[0]
        uint64_t i
        integral j
        double sum_

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


def assign_rows_csr(
    X,
    const intptr_t[:] X_rows,
    const intptr_t[:] out_rows,
    floating[:, ::1] out,
):
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
        # intptr_t (npy_intp, np.intp in Python) is what np.where returns,
        # but int is what scipy.sparse uses.
        intp_t i, ind, j, k
        intptr_t rX
        const floating[:] data = X.data
        const int32_t[:] indices = X.indices
        const int32_t[:] indptr = X.indptr

    if X_rows.shape[0] != out_rows.shape[0]:
        raise ValueError("cannot assign %d rows to %d"
                         % (X_rows.shape[0], out_rows.shape[0]))

    with nogil:
        for k in range(out_rows.shape[0]):
            out[out_rows[k]] = 0.0

        for i in range(X_rows.shape[0]):
            rX = X_rows[i]
            for ind in range(indptr[rX], indptr[rX + 1]):
                j = indices[ind]
                out[out_rows[i], j] = data[ind]
