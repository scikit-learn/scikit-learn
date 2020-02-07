import numpy as np


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if X.dtype.kind in ("i", "u"):
        # can't have NaNs in integer array.
        return np.zeros(X.shape, dtype=bool)
    else:
        return np.isnan(X)
