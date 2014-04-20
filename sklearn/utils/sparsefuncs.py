# Authors: Manoj Kumar

# License: BSD 3 clause
import scipy.sparse as sp

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
