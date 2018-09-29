# Authors: Manoj Kumar
#          Thomas Unterthiner
#          Giorgio Patrini
#
# License: BSD 3 clause
import scipy.sparse as sp
import numpy as np

from .sparsefuncs_fast import (
    csr_mean_variance_axis0 as _csr_mean_var_axis0,
    csc_mean_variance_axis0 as _csc_mean_var_axis0,
    incr_mean_variance_axis0 as _incr_mean_var_axis0)


def _raise_typeerror(X):
    """Raises a TypeError if X is not a CSR or CSC matrix"""
    input_type = X.format if sp.issparse(X) else type(X)
    err = "Expected a CSR or CSC sparse matrix, got %s." % input_type
    raise TypeError(err)


def _raise_error_wrong_axis(axis):
    if axis not in (0, 1):
        raise ValueError(
            "Unknown axis value: %d. Use 0 for rows, or 1 for columns" % axis)


def inplace_csr_column_scale(X, scale):
    """Inplace column scaling of a CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X : CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale : float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    assert scale.shape[0] == X.shape[1]
    X.data *= scale.take(X.indices, mode='clip')


def inplace_csr_row_scale(X, scale):
    """ Inplace row scaling of a CSR matrix.

    Scale each sample of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X : CSR sparse matrix, shape (n_samples, n_features)
        Matrix to be scaled.

    scale : float array with shape (n_samples,)
        Array of precomputed sample-wise values to use for scaling.
    """
    assert scale.shape[0] == X.shape[0]
    X.data *= np.repeat(scale, np.diff(X.indptr))


def mean_variance_axis(X, axis):
    """Compute mean and variance along an axix on a CSR or CSC matrix

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    axis : int (either 0 or 1)
        Axis along which the axis should be computed.

    Returns
    -------

    means : float array with shape (n_features,)
        Feature-wise means

    variances : float array with shape (n_features,)
        Feature-wise variances

    """
    _raise_error_wrong_axis(axis)

    if isinstance(X, sp.csr_matrix):
        if axis == 0:
            return _csr_mean_var_axis0(X)
        else:
            return _csc_mean_var_axis0(X.T)
    elif isinstance(X, sp.csc_matrix):
        if axis == 0:
            return _csc_mean_var_axis0(X)
        else:
            return _csr_mean_var_axis0(X.T)
    else:
        _raise_typeerror(X)


def incr_mean_variance_axis(X, axis, last_mean, last_var, last_n):
    """Compute incremental mean and variance along an axix on a CSR or
    CSC matrix.

    last_mean, last_var are the statistics computed at the last step by this
    function. Both must be initialized to 0-arrays of the proper size, i.e.
    the number of features in X. last_n is the number of samples encountered
    until now.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    axis : int (either 0 or 1)
        Axis along which the axis should be computed.

    last_mean : float array with shape (n_features,)
        Array of feature-wise means to update with the new data X.

    last_var : float array with shape (n_features,)
        Array of feature-wise var to update with the new data X.

    last_n : int with shape (n_features,)
        Number of samples seen so far, excluded X.

    Returns
    -------

    means : float array with shape (n_features,)
        Updated feature-wise means.

    variances : float array with shape (n_features,)
        Updated feature-wise variances.

    n : int with shape (n_features,)
        Updated number of seen samples.

    Notes
    -----
    NaNs are ignored in the algorithm.

    """
    _raise_error_wrong_axis(axis)

    if isinstance(X, sp.csr_matrix):
        if axis == 0:
            return _incr_mean_var_axis0(X, last_mean=last_mean,
                                        last_var=last_var, last_n=last_n)
        else:
            return _incr_mean_var_axis0(X.T, last_mean=last_mean,
                                        last_var=last_var, last_n=last_n)
    elif isinstance(X, sp.csc_matrix):
        if axis == 0:
            return _incr_mean_var_axis0(X, last_mean=last_mean,
                                        last_var=last_var, last_n=last_n)
        else:
            return _incr_mean_var_axis0(X.T, last_mean=last_mean,
                                        last_var=last_var, last_n=last_n)
    else:
        _raise_typeerror(X)


def inplace_column_scale(X, scale):
    """Inplace column scaling of a CSC/CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X : CSC or CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale : float array with shape (n_features,)
        Array of precomputed feature-wise values to use for scaling.
    """
    if isinstance(X, sp.csc_matrix):
        inplace_csr_row_scale(X.T, scale)
    elif isinstance(X, sp.csr_matrix):
        inplace_csr_column_scale(X, scale)
    else:
        _raise_typeerror(X)


def inplace_row_scale(X, scale):
    """ Inplace row scaling of a CSR or CSC matrix.

    Scale each row of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape (n_samples, n_features)
        Matrix to be scaled.

    scale : float array with shape (n_features,)
        Array of precomputed sample-wise values to use for scaling.
    """
    if isinstance(X, sp.csc_matrix):
        inplace_csr_column_scale(X.T, scale)
    elif isinstance(X, sp.csr_matrix):
        inplace_csr_row_scale(X, scale)
    else:
        _raise_typeerror(X)


def inplace_swap_row_csc(X, m, n):
    """
    Swaps two rows of a CSC matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m : int
        Index of the row of X to be swapped.

    n : int
        Index of the row of X to be swapped.
    """
    for t in [m, n]:
        if isinstance(t, np.ndarray):
            raise TypeError("m and n should be valid integers")

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    m_mask = X.indices == m
    X.indices[X.indices == n] = m
    X.indices[m_mask] = n


def inplace_swap_row_csr(X, m, n):
    """
    Swaps two rows of a CSR matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csr_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m : int
        Index of the row of X to be swapped.

    n : int
        Index of the row of X to be swapped.
    """
    for t in [m, n]:
        if isinstance(t, np.ndarray):
            raise TypeError("m and n should be valid integers")

    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    # The following swapping makes life easier since m is assumed to be the
    # smaller integer below.
    if m > n:
        m, n = n, m

    indptr = X.indptr
    m_start = indptr[m]
    m_stop = indptr[m + 1]
    n_start = indptr[n]
    n_stop = indptr[n + 1]
    nz_m = m_stop - m_start
    nz_n = n_stop - n_start

    if nz_m != nz_n:
        # Modify indptr first
        X.indptr[m + 2:n] += nz_n - nz_m
        X.indptr[m + 1] = m_start + nz_n
        X.indptr[n] = n_stop - nz_m

    X.indices = np.concatenate([X.indices[:m_start],
                                X.indices[n_start:n_stop],
                                X.indices[m_stop:n_start],
                                X.indices[m_start:m_stop],
                                X.indices[n_stop:]])
    X.data = np.concatenate([X.data[:m_start],
                             X.data[n_start:n_stop],
                             X.data[m_stop:n_start],
                             X.data[m_start:m_stop],
                             X.data[n_stop:]])


def inplace_swap_row(X, m, n):
    """
    Swaps two rows of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m : int
        Index of the row of X to be swapped.

    n : int
        Index of the row of X to be swapped.
    """
    if isinstance(X, sp.csc_matrix):
        inplace_swap_row_csc(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        inplace_swap_row_csr(X, m, n)
    else:
        _raise_typeerror(X)


def inplace_swap_column(X, m, n):
    """
    Swaps two columns of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape=(n_samples, n_features)
        Matrix whose two columns are to be swapped.

    m : int
        Index of the column of X to be swapped.

    n : int
        Index of the column of X to be swapped.
    """
    if m < 0:
        m += X.shape[1]
    if n < 0:
        n += X.shape[1]
    if isinstance(X, sp.csc_matrix):
        inplace_swap_row_csr(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        inplace_swap_row_csc(X, m, n)
    else:
        _raise_typeerror(X)


def _minor_reduce(X, ufunc):
    major_index = np.flatnonzero(np.diff(X.indptr))
    value = ufunc.reduceat(X.data, X.indptr[major_index])
    return major_index, value


def _min_or_max_axis(X, axis, min_or_max):
    N = X.shape[axis]
    if N == 0:
        raise ValueError("zero-size array to reduction operation")
    M = X.shape[1 - axis]
    mat = X.tocsc() if axis == 0 else X.tocsr()
    mat.sum_duplicates()
    major_index, value = _minor_reduce(mat, min_or_max)
    not_full = np.diff(mat.indptr)[major_index] < N
    value[not_full] = min_or_max(value[not_full], 0)
    mask = value != 0
    major_index = np.compress(mask, major_index)
    value = np.compress(mask, value)

    if axis == 0:
        res = sp.coo_matrix((value, (np.zeros(len(value)), major_index)),
                            dtype=X.dtype, shape=(1, M))
    else:
        res = sp.coo_matrix((value, (major_index, np.zeros(len(value)))),
                            dtype=X.dtype, shape=(M, 1))
    return res.A.ravel()


def _sparse_min_or_max(X, axis, min_or_max):
    if axis is None:
        if 0 in X.shape:
            raise ValueError("zero-size array to reduction operation")
        zero = X.dtype.type(0)
        if X.nnz == 0:
            return zero
        m = min_or_max.reduce(X.data.ravel())
        if X.nnz != np.product(X.shape):
            m = min_or_max(zero, m)
        return m
    if axis < 0:
        axis += 2
    if (axis == 0) or (axis == 1):
        return _min_or_max_axis(X, axis, min_or_max)
    else:
        raise ValueError("invalid axis, use 0 for rows, or 1 for columns")


def _sparse_min_max(X, axis):
        return (_sparse_min_or_max(X, axis, np.minimum),
                _sparse_min_or_max(X, axis, np.maximum))


def _sparse_nan_min_max(X, axis):
    return(_sparse_min_or_max(X, axis, np.fmin),
           _sparse_min_or_max(X, axis, np.fmax))


def min_max_axis(X, axis, ignore_nan=False):
    """Compute minimum and maximum along an axis on a CSR or CSC matrix and
    optionally ignore NaN values.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    axis : int (either 0 or 1)
        Axis along which the axis should be computed.

    ignore_nan : bool, default is False
        Ignore or passing through NaN values.

        .. versionadded:: 0.20

    Returns
    -------

    mins : float array with shape (n_features,)
        Feature-wise minima

    maxs : float array with shape (n_features,)
        Feature-wise maxima
    """
    if isinstance(X, sp.csr_matrix) or isinstance(X, sp.csc_matrix):
        if ignore_nan:
            return _sparse_nan_min_max(X, axis=axis)
        else:
            return _sparse_min_max(X, axis=axis)
    else:
        _raise_typeerror(X)


def count_nonzero(X, axis=None, sample_weight=None):
    """A variant of X.getnnz() with extension to weighting on axis 0

    Useful in efficiently calculating multilabel metrics.

    Parameters
    ----------
    X : CSR sparse matrix, shape = (n_samples, n_labels)
        Input data.

    axis : None, 0 or 1
        The axis on which the data is aggregated.

    sample_weight : array, shape = (n_samples,), optional
        Weight for each row of X.
    """
    if axis == -1:
        axis = 1
    elif axis == -2:
        axis = 0
    elif X.format != 'csr':
        raise TypeError('Expected CSR sparse format, got {0}'.format(X.format))

    # We rely here on the fact that np.diff(Y.indptr) for a CSR
    # will return the number of nonzero entries in each row.
    # A bincount over Y.indices will return the number of nonzeros
    # in each column. See ``csr_matrix.getnnz`` in scipy >= 0.14.
    if axis is None:
        if sample_weight is None:
            return X.nnz
        else:
            return np.dot(np.diff(X.indptr), sample_weight)
    elif axis == 1:
        out = np.diff(X.indptr)
        if sample_weight is None:
            return out
        return out * sample_weight
    elif axis == 0:
        if sample_weight is None:
            return np.bincount(X.indices, minlength=X.shape[1])
        else:
            weights = np.repeat(sample_weight, np.diff(X.indptr))
            return np.bincount(X.indices, minlength=X.shape[1],
                            weights=weights)
    else:
        raise ValueError('Unsupported axis: {0}'.format(axis))


def _get_median(data, n_zeros):
    """Compute the median of data with n_zeros additional zeros.

    This function is used to support sparse matrices; it modifies data in-place
    """
    n_elems = len(data) + n_zeros
    if not n_elems:
        return np.nan
    n_negative = np.count_nonzero(data < 0)
    middle, is_odd = divmod(n_elems, 2)
    data.sort()

    if is_odd:
        return _get_elem_at_rank(middle, data, n_negative, n_zeros)

    return (_get_elem_at_rank(middle - 1, data, n_negative, n_zeros) +
            _get_elem_at_rank(middle, data, n_negative, n_zeros)) / 2.


def _get_elem_at_rank(rank, data, n_negative, n_zeros):
    """Find the value in data augmented with n_zeros for the given rank"""
    if rank < n_negative:
        return data[rank]
    if rank - n_negative < n_zeros:
        return 0
    return data[rank - n_zeros]


def csc_median_axis_0(X):
    """Find the median across axis 0 of a CSC matrix.
    It is equivalent to doing np.median(X, axis=0).

    Parameters
    ----------
    X : CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------
    median : ndarray, shape (n_features,)
        Median.

    """
    if not isinstance(X, sp.csc_matrix):
        raise TypeError("Expected matrix of CSC format, got %s" % X.format)

    indptr = X.indptr
    n_samples, n_features = X.shape
    median = np.zeros(n_features)

    for f_ind, (start, end) in enumerate(zip(indptr[:-1], indptr[1:])):

        # Prevent modifying X in place
        data = np.copy(X.data[start: end])
        nz = n_samples - data.size
        median[f_ind] = _get_median(data, nz)

    return median
