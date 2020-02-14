import numpy as np
from scipy import sparse

from . import is_scalar_nan
from .fixes import _object_dtype_isnan


def _fit_mask(X, value_to_mask):
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
        Xt == value_to_mask
    return Xt


def _get_mask(X, value_to_mask, reconstruct_sparse=False):
    """Compute the boolean mask X == missing_values."""
    # We get entire sparse matrix when reconstruct is True
    # Otherwise we get X.data if X is sparse, and X when X is dense
    if sparse.issparse(X) and reconstruct_sparse:
        Xt = _fit_mask(X.data, value_to_mask)
    # if we need not reconstruct sparse or X is dense we can directly
    # convert to mask
    else:
        return _fit_mask(X, value_to_mask)

    # following code will execute only when we need to reconstruct sparse
    sparse_constructor = (sparse.csr_matrix if X.format == 'csr'
                          else sparse.csc_matrix)
    re_sparse = sparse_constructor(
        (Xt, X.indices.copy(), X.indptr.copy()), shape=X.shape, dtype=bool
    )

    return Xt, re_sparse
