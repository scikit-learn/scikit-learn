# Authors: Manoj Kumar
#          Thomas Unterthiner

# License: BSD 3 clause
import scipy.sparse as sp
import numpy as np

from sklearn.utils.random import choice
from .fixes import sparse_min_max
from .sparsefuncs_fast import (csr_mean_variance_axis0,
                               csc_mean_variance_axis0)


def _raise_typeerror(X):
    """Raises a TypeError if X is not a CSR or CSC matrix"""
    input_type = X.format if sp.issparse(X) else type(X)
    err = "Expected a CSR or CSC sparse matrix, got %s." % input_type
    raise TypeError(err)


def inplace_csr_column_scale(X, scale):
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
    assert scale.shape[0] == X.shape[1]
    X.data *= scale.take(X.indices, mode='clip')


def inplace_csr_row_scale(X, scale):
    """ Inplace row scaling of a CSR matrix.

    Scale each sample of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSR sparse matrix, shape (n_samples, n_features)
    matrix to be scaled.

    scale: float array with shape (n_samples,)
    Array of precomputed sample-wise values to use for scaling.
    """
    assert scale.shape[0] == X.shape[0]
    X.data *= np.repeat(scale, np.diff(X.indptr))


def mean_variance_axis0(X):
    """Compute mean and variance along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    means: float array with shape (n_features,)
        Feature-wise means

    variances: float array with shape (n_features,)
        Feature-wise variances

    """
    if isinstance(X, sp.csr_matrix):
        return csr_mean_variance_axis0(X)
    elif isinstance(X, sp.csc_matrix):
        return csc_mean_variance_axis0(X)
    else:
        _raise_typeerror(X)


def inplace_column_scale(X, scale):
    """Inplace column scaling of a CSC/CSR matrix.

    Scale each feature of the data matrix by multiplying with specific scale
    provided by the caller assuming a (n_samples, n_features) shape.

    Parameters
    ----------
    X: CSC or CSR matrix with shape (n_samples, n_features)
        Matrix to normalize using the variance of the features.

    scale: float array with shape (n_features,)
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
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
    matrix to be scaled.

    scale: float array with shape (n_features,)
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
    X: scipy.sparse.csc_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m: int
        Index of the row of X to be swapped.

    n: int
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
    X: scipy.sparse.csr_matrix, shape=(n_samples, n_features)
        Matrix whose two rows are to be swapped.

    m: int
        Index of the row of X to be swapped.

    n: int
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

    m: int
        Index of the row of X to be swapped.

    n: int
        Index of the row of X to be swapped.
    """
    if isinstance(X, sp.csc_matrix):
        return inplace_swap_row_csc(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return inplace_swap_row_csr(X, m, n)
    else:
        _raise_typeerror(X)


def inplace_swap_column(X, m, n):
    """
    Swaps two columns of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : CSR or CSC sparse matrix, shape=(n_samples, n_features)
        Matrix whose two columns are to be swapped.

    m: int
        Index of the column of X to be swapped.

    n : int
        Index of the column of X to be swapped.
    """
    if m < 0:
        m += X.shape[1]
    if n < 0:
        n += X.shape[1]
    if isinstance(X, sp.csc_matrix):
        return inplace_swap_row_csr(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return inplace_swap_row_csc(X, m, n)
    else:
        _raise_typeerror(X)


def min_max_axis(X, axis):
    """Compute minimum and maximum along axis 0 on a CSR or CSC matrix

    Parameters
    ----------
    X: CSR or CSC sparse matrix, shape (n_samples, n_features)
        Input data.

    Returns
    -------

    mins: float array with shape (n_features,)
        Feature-wise minima

    maxs: float array with shape (n_features,)
        Feature-wise maxima
    """
    if isinstance(X, sp.csr_matrix) or isinstance(X, sp.csc_matrix):
        return sparse_min_max(X, axis=axis)
    else:
        _raise_typeerror(X)


def random_choice_csc(a, size, p=None):
    # XXX Corner case: size 0

    if isinstance(a, int):
        a = np.arange(a)
    a = np.array(a, np.int)

    if p is None:
            p = np.empty(shape=a.shape[0])
            p.fill(1.0/a.shape[0])
    p = np.array(p)

    # XXX Validate: p and a have the same length
    # XXX Warning: Sparse efficieny if 0 is not in a, or p for 0 is small

    if 0 not in a:
        a = np.insert(a, 0, 0)
        p = np.insert(p, 0, 0.0)

    if 0 in a and a.shape[0] is 1:
        return sp.csc_matrix((size, 1))

    nnz = size - int(size * p[np.where(a == 0)[0][0]])  # XXX maybe round
    indices = choice(a=range(size), size=nnz, replace=False)

    # Normalize probabilites for the nonzero elements
    p_nz = p[a != 0]
    p_nz_norm = p_nz / np.sum(p_nz)
    data = choice(a=a[a != 0], size=nnz, p=p_nz_norm, replace=True)
    indptr = [0, indices.shape[0]]

    return sp.csc_matrix((data, indices, indptr), shape=(size, 1))
