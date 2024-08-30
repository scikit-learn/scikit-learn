# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from scipy.sparse import issparse

from sklearn.utils._array_api import _is_numpy_namespace, get_namespace


def attach_unique(y, xp=None):
    """Attach unique values of y to y and return the result.

    The result is a view of y, and the metadata (unique) is not attached to y.
    """
    try:
        # avoid recalculating unique in nested calls.
        if "unique" in y.dtype.metadata:
            return y
    except AttributeError:
        pass

    xp, _ = get_namespace(y, xp=xp)
    if not _is_numpy_namespace(xp) or issparse(y):
        # we don't support non-numpy arrays here
        return y

    unique = np.unique_values(y)
    try:
        unique_dtype = np.dtype(y.dtype, metadata={"unique": unique})
        return y.view(dtype=unique_dtype)
    except AttributeError:
        # if y is not a numpy array, we can't attach metadata
        return y


def _cached_unique(y, xp=None):
    """Return the unique values of y.

    Use the cached values from dtype.metadata if present.

    This function does NOT cache the values in y, i.e. it doesn't change y.

    Call `attach_unique` to attach the unique values to y.
    """
    try:
        if y.dtype.metadata is not None and "unique" in y.dtype.metadata:
            return y.dtype.metadata["unique"]
    except AttributeError:
        # in case y is not a numpy array
        pass
    xp, _ = get_namespace(y, xp=xp)
    return xp.unique_values(y)


def cached_unique(*ys, xp=None):
    """Return the unique values of ys.

    Use the cached values from dtype.metadata if present.

    This function does NOT cache the values in y, i.e. it doesn't change y.

    Call `attach_unique` to attach the unique values to y.
    """
    return (_cached_unique(y, xp=xp) for y in ys)
