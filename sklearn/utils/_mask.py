import numpy as np
from scipy import sparse

from . import is_scalar_nan
from .fixes import _object_dtype_isnan


def _get_mask(X, value_to_mask, reconstruct_sparse=False):
    """Compute the boolean mask X == missing_values."""
    X_ = X.copy()
    if reconstruct_sparse:
        # if sparse.issparse(X):
        X_ = X.data

    if is_scalar_nan(value_to_mask):
        if X_.dtype.kind == "f":
            Xt = np.isnan(X_)
        elif X.dtype.kind in ("i", "u"):
            # can't have NaNs in integer array.
            Xt = np.zeros(X_.shape, dtype=bool)
        else:
            # np.isnan does not work on object dtypes.
            Xt = _object_dtype_isnan(X_)
    else:
        # X == value_to_mask with object dtypes does not always perform
        # element-wise for old versions of numpy
        Xt = np.equal(X_, value_to_mask)

    if not reconstruct_sparse:
        return Xt

    else:
        sparse_constructor = (sparse.csr_matrix
                              if X.format == 'csr'
                              else sparse.csc_matrix)
        re_sparse = sparse_constructor(
            (Xt, X.indices.copy(), X.indptr.copy()),
            shape=X.shape, dtype=bool)

        return Xt, re_sparse
