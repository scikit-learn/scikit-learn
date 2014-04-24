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


def swap_row_csc(X, m, n):
    """
    Swaps two rows of a CSC matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first sample
    n : int, index of second sample
    """
    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]

    m_mask = X.indices == m
    X.indices[X.indices == n] = m
    X.indices[m_mask] = n


def swap_row_csr(X, m, n):
    """
    Swaps two rows of a CSR matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first sample
    n : int, index of second sample
    """
    if m < 0:
        m += X.shape[0]
    if n < 0:
        n += X.shape[0]
    if m > n:
        m, n = n, m

    indptr = X.indptr
    m_ptr1 = indptr[m]
    m_ptr2 = indptr[m + 1]
    n_ptr1 = indptr[n]
    n_ptr2 = indptr[n + 1]
    nz_m = m_ptr2 - m_ptr1
    nz_n = n_ptr2 - n_ptr1


    if nz_m != nz_n:
        # Modify indptr first
        X.indptr[m + 2:n] += nz_n - nz_m
        X.indptr[m + 1] = X.indptr[m] + nz_n
        X.indptr[n] = X.indptr[n + 1] - nz_m

    X.indices = np.concatenate([X.indices[:m_ptr1], X.indices[n_ptr1:n_ptr2],
                                X.indices[m_ptr2:n_ptr1],
                                X.indices[m_ptr1:m_ptr2],
                                X.indices[n_ptr2:]])
    X.data = np.concatenate([X.data[:m_ptr1], X.data[n_ptr1:n_ptr2],
                             X.data[m_ptr2:n_ptr1], X.data[m_ptr1:m_ptr2],
                             X.data[n_ptr2:]])


def swap_row(X, m, n):
    """
    Swaps two rows of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first sample
    n : int, index of second sample
    """
    if isinstance(X, sp.csc_matrix):
        return swap_row_csc(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return swap_row_csr(X, m, n)
    else:
        raise TypeError(
            "Unsupported type; expected a CSR or CSC sparse matrix.")


def swap_column(X, m, n):
    """
    Swaps two columns of a CSC/CSR matrix in-place.

    Parameters
    ----------
    X : scipy.sparse.csc_matrix, shape=(n_samples, n_features)
    m : int, index of first sample
    n : int, index of second sample
    """
    if m < 0:
        m += X.shape[1]
    if n < 0:
        n += X.shape[1]
    if isinstance(X, sp.csc_matrix):
        return swap_row_csr(X, m, n)
    elif isinstance(X, sp.csr_matrix):
        return swap_row_csc(X, m, n)
    else:
        raise TypeError(
            "Unsupported type; expected a CSR or CSC sparse matrix.")
