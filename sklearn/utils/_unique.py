# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.utils._array_api import get_namespace


def _get_metadata(y):
    """Return dtype metadata, or None when dtype metadata is unsupported."""
    if not isinstance(y, np.ndarray) or y.dtype.kind == "T":
        # StringDType cannot carry arbitrary metadata, including on NumPy 2.5.
        return None
    return y.dtype.metadata or {}


def _attach_metadata(y, **metadata):
    """Return a view with metadata, or the input if metadata is unsupported."""
    current_metadata = _get_metadata(y)
    if current_metadata is None:
        return y
    dtype = np.dtype(y.dtype, metadata={**current_metadata, **metadata})
    return y.view(dtype=dtype)


def _attach_unique(y):
    """Attach unique values when the dtype supports metadata, without mutating y."""
    metadata = _get_metadata(y)
    if metadata is None or "unique" in metadata:
        return y
    return _attach_metadata(y, unique=np.unique(y))


def attach_unique(*ys, return_tuple=False):
    """Attach unique values of ys to ys and return the results.

    For NumPy dtypes supporting metadata, the result is a view of y with
    cached unique values. Other inputs are returned unchanged. The input is
    never modified.

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
    metadata = _get_metadata(y)
    if metadata is not None and "unique" in metadata:
        return metadata["unique"]
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
