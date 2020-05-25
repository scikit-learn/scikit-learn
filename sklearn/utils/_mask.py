import numpy as np
from scipy import sparse as sp

from . import is_scalar_nan
from .fixes import _object_dtype_isnan


def _get_dense_mask(X, value_to_mask):
    if is_scalar_nan(value_to_mask):
        if X.dtype.kind == "f":
            Xt = np.isnan(X)
        elif X.dtype.kind in ("i", "u"):
            # can't have NaNs in integer array.
            Xt = np.zeros(X.shape, dtype=bool)
        else:
            # np.isnan does not work on object dtypes.
            Xt = _object_dtype_isnan(X)
    else:
        Xt = X == value_to_mask

    return Xt


def _get_mask(X, value_to_mask, sparse=False):
    """Compute the boolean mask X == value_to_mask.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Input data, where ``n_samples`` is the number of samples and
        ``n_features`` is the number of features.

    value_to_mask : {int, float, nan}
                    The value which is to be masked in X.

    sparse : bool, default=False
             Whether or not we need to reconstruct sparse matrix.
             If True, X is considered considered sparse and sparse
             mask matrix is created.
             If False, a dense mask matrix is returned of the same
             shape as X.

    """
    if not (sp.issparse(X) and sparse):
        # For all cases apart of a sparse input where we need to reconstruct
        # a sparse output
        return _fit_mask(X, value_to_mask)

    if not sp.issparse(X):
        raise ValueError("Cannot reconstruct sparse mask as passed"
                         " input is not sparse.")
    Xt = _get_dense_mask(X.data, value_to_mask)

    sparse_constructor = (sp.csr_matrix if X.format == 'csr'
                          else sp.csc_matrix)
    Xt_sparse = sparse_constructor(
        (Xt, X.indices.copy(), X.indptr.copy()), shape=X.shape, dtype=bool
    )

    return Xt_sparse
