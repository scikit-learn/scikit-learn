# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.utils._array_api import get_namespace


def _attach_unique(y):
    """Attach unique values of y to y and return the result.

    The result is a view of y, and the metadata (unique) is not attached to y.
    """
    if not isinstance(y, np.ndarray):
        return y
    try:
        # avoid recalculating unique in nested calls.
        if "unique" in y.dtype.metadata:
            return y
    except (AttributeError, TypeError):
        pass

    unique = np.unique(y)
    unique_dtype = np.dtype(y.dtype, metadata={"unique": unique})
    return y.view(dtype=unique_dtype)


def attach_unique(*ys, return_tuple=False):
    """Attach unique values of ys to ys and return the results.

    The result is a view of y, and the metadata (unique) is not attached to y.

    IMPORTANT: The output of this function should NEVER be returned in functions.
    This is to avoid this pattern:

    .. code:: python

        y = np.array([1, 2, 3])
        y = attach_unique(y)
        y[1] = -1
        # now np.unique(y) will be different from cached_unique(y)

    Parameters
    ----------
    *ys : sequence of array-like
        Input data arrays.

    return_tuple : bool, default=False
        If True, always return a tuple even if there is only one array.

    Returns
    -------
    ys : tuple of array-like or array-like
        Input data with unique values attached.
    """
    res = tuple(_attach_unique(y) for y in ys)
    if len(res) == 1 and not return_tuple:
        return res[0]
    return res


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

    Parameters
    ----------
    *ys : sequence of array-like
        Input data arrays.

    xp : module, default=None
        Precomputed array namespace module. When passed, typically from a caller
        that has already performed inspection of its own inputs, skips array
        namespace inspection.

    Returns
    -------
    res : tuple of array-like or array-like
        Unique values of ys.
    """
    res = tuple(_cached_unique(y, xp=xp) for y in ys)
    if len(res) == 1:
        return res[0]
    return res
