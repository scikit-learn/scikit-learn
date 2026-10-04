# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from contextlib import contextmanager
from contextvars import ContextVar

import numpy as np

from sklearn.utils._array_api import get_namespace

_METADATA_CACHE = ContextVar("sklearn_array_metadata", default=None)


class _MetadataCache(dict):
    def __init__(self):
        super().__init__()
        self.active = True


def _active_metadata_cache():
    cache = _METADATA_CACHE.get()
    return cache if cache is not None and cache.active else None


@contextmanager
def _metadata_cache():
    """Share metadata during a read-only operation, including nested calls.

    Entries retain their input objects to prevent identity reuse. The outermost
    scope releases all entries on exit, including when validation raises. Inputs
    must not be mutated inside the scope. Nothing is cached between operations.
    """
    if _active_metadata_cache() is not None:
        yield
        return
    cache = _MetadataCache()
    token = _METADATA_CACHE.set(cache)
    try:
        yield
    finally:
        # Copied contexts may outlive the operation. Release retained arrays and
        # prevent those contexts from treating this closed scope as active.
        cache.active = False
        cache.clear()
        _METADATA_CACHE.reset(token)


def _remember_metadata(y, metadata):
    cache = _active_metadata_cache()
    if cache is not None:
        cache[id(y)] = (y, metadata)


def _get_metadata(y):
    """Read operation-local metadata, falling back to NumPy dtype metadata."""
    cache = _active_metadata_cache()
    if cache is not None and id(y) in cache:
        return cache[id(y)][1]
    if isinstance(y, np.ndarray) and y.dtype.kind != "T":
        return y.dtype.metadata or {}
    return {} if cache is not None else None


def _attach_metadata(y, **metadata):
    """Attach metadata through a NumPy view or the active operation-local cache."""
    current_metadata = _get_metadata(y)
    if current_metadata is None:
        return y
    metadata = {**current_metadata, **metadata}
    _remember_metadata(y, metadata)
    if not isinstance(y, np.ndarray) or y.dtype.kind == "T":
        return y
    try:
        dtype = np.dtype(y.dtype, metadata=metadata)
    except TypeError:
        # Other custom NumPy dtypes can also lack metadata support.
        return y
    view = y.view(dtype=dtype)
    _remember_metadata(view, metadata)
    return view


def _transfer_unique(source, target):
    """Reuse unique values after a conversion known to preserve every value."""
    metadata = _get_metadata(source)
    if metadata is not None and "unique" in metadata:
        current = _get_metadata(target) or {}
        _remember_metadata(target, {**current, "unique": metadata["unique"]})


def _attach_unique(y):
    """Cache NumPy unique values using dtype metadata or the active scope."""
    if not isinstance(y, np.ndarray):
        return y
    metadata = _get_metadata(y)
    if metadata is None or "unique" in metadata:
        return y
    return _attach_metadata(y, unique=np.unique(y))


def attach_unique(*ys, return_tuple=False):
    """Attach unique values of ys to ys and return the results.

    For NumPy dtypes supporting metadata, the result is a view of y with
    cached unique values. Within ``_metadata_cache``, StringDType arrays use
    operation-local metadata and are returned unchanged. Other inputs are
    returned unchanged. The input is never modified.

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

    Use operation-local or dtype metadata if present.

    The input is never modified. Within ``_metadata_cache``, a computed result
    is retained until the outermost read-only operation finishes.

    Outside a cache scope, call `attach_unique` to cache on supported NumPy dtypes.
    """
    metadata = _get_metadata(y)
    if metadata is not None and "unique" in metadata:
        return metadata["unique"]
    xp, _ = get_namespace(y, xp=xp)
    unique = xp.unique_values(y)
    _remember_metadata(y, {**(metadata or {}), "unique": unique})
    return unique


def cached_unique(*ys, xp=None):
    """Return the unique values of ys.

    Use operation-local or dtype metadata if present.

    The input is never modified. Within ``_metadata_cache``, a computed result
    is retained until the outermost read-only operation finishes.

    Outside a cache scope, call `attach_unique` to cache on supported NumPy dtypes.

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
