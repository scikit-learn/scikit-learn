"""Helpers for preprocessing"""

import numpy as np
from scipy import sparse

from ..utils import check_array
from ..utils.fixes import _astype_copy_false
from ..utils.validation import FLOAT_DTYPES



    X = check_array(X, accept_sparse='csc', copy=copy, dtype=FLOAT_DTYPES)

    if sparse.issparse(X) and retain_order:
        raise ValueError("The retain_order option can only be set to True "
                         "for dense matrices.")

    if isinstance(selected, str) and selected == "all":
        return transform(X)

    if len(selected) == 0:
        return X

    n_features = X.shape[1]
    ind = np.arange(n_features)
    sel = np.zeros(n_features, dtype=bool)
    sel[np.asarray(selected)] = True
    not_sel = np.logical_not(sel)
    n_selected = np.sum(sel)

    if n_selected == 0:
        # No features selected.
        return X
    elif n_selected == n_features:
        # All features selected.
        return transform(X)
    else:
        X_sel = transform(X[:, ind[sel]])
        # The columns of X which are not transformed need
        # to be casted to the desire dtype before concatenation.
        # Otherwise, the stacking will cast to the higher-precision dtype.
        X_not_sel = X[:, ind[not_sel]].astype(dtype, **_astype_copy_false(X))

    if retain_order:
        if X_sel.shape[1] + X_not_sel.shape[1] != n_features:
            raise ValueError("The retain_order option can only be set to True "
                             "if the dimensions of the input array match the "
                             "dimensions of the transformed array.")

        # Fancy indexing not supported for sparse matrices
        X[:, ind[sel]] = X_sel
        return X

    if sparse.issparse(X_sel) or sparse.issparse(X_not_sel):
        return sparse.hstack((X_sel, X_not_sel))
    else:
        return np.hstack((X_sel, X_not_sel))
