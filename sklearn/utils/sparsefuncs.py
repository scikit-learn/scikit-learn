# Authors: Manoj Kumar

# License: BSD 3 clause
import scipy.sparse as sp
import numpy as np

from .sparsefuncs_fast import (csr_mean_variance_axis0,
                               csc_mean_variance_axis0,
                               inplace_csr_column_scale,
                               inplace_csc_column_scale)

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
        raise TypeError(
                "Unsupported type; expected a CSR or CSC sparse matrix.")


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
    if isinstance(X, sp.csr_matrix):
        return inplace_csr_column_scale(X, scale)
    elif isinstance(X, sp.csc_matrix):
        return inplace_csc_column_scale(X, scale)
    else:
        raise TypeError(
                "Unsupported type; expected a CSR or CSC sparse matrix.")


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
        raise TypeError(
            "Unsupported type; expected a CSR or CSC sparse matrix.")


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
        raise TypeError(
            "Unsupported type; expected a CSR or CSC sparse matrix.")
