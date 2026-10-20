"""Delegation layer for searching functions."""

from typing import Literal

from . import _agnostic, default_dtype
from ._lib import _compat
from ._lib._typing import Array, ArrayNamespace

__all__ = ["searchsorted"]


def searchsorted(
    x1: Array,
    x2: Array,
    /,
    *,
    side: Literal["left", "right"] = "left",
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Find indices where elements should be inserted to maintain order.

    Find the indices into a sorted array ``x1`` such that if the elements in ``x2``
    were inserted before the indices, the resulting array would remain sorted.

    The behavior of this function is similar to that of :func:`array_api.searchsorted`,
    but it relaxes the requirement that `x1` must be one-dimensional.
    This function is vectorized, treating slices along the last axis
    as elements and preceding axes as batch (or "loop") dimensions.

    Parameters
    ----------
    x1 : Array
        Input array. Should have a real-valued data type. Must be sorted in ascending
        order along the last axis.
    x2 : Array
        Array containing search values. Should have a real-valued data type. Must have
        the same shape as ``x1`` except along the last axis.
    side : {'left', 'right'}, optional
        Argument controlling which index is returned if an element of ``x2`` is equal to
        one or more elements of ``x1``: ``'left'`` returns the index of the first of
        these elements; ``'right'`` returns the next index after the last of these
        elements. Default: ``'left'``.
    xp : array_namespace, optional
        The standard-compatible namespace for the array arguments. Default: infer.

    Returns
    -------
    Array: integer array
        An array of indices with the same shape as ``x2``.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> x = xp.asarray([11, 12, 13, 13, 14, 15])
    >>> xpx.searchsorted(x, xp.asarray([10, 11.5, 14.5, 16]), xp=xp)
    Array([0, 1, 5, 6], dtype=array_api_strict.int64)
    >>> xpx.searchsorted(x, xp.asarray(13), xp=xp)
    Array(2, dtype=array_api_strict.int64)
    >>> xpx.searchsorted(x, xp.asarray(13), side='right', xp=xp)
    Array(4, dtype=array_api_strict.int64)

    `searchsorted` is vectorized along the last axis.

    >>> x1 = xp.asarray([[1., 2., 3., 4.], [5., 6., 7., 8.]])
    >>> x2 = xp.asarray([[1.1, 3.3], [6.6, 8.8]])
    >>> xpx.searchsorted(x1, x2, xp=xp)
    Array([[1, 3],
           [2, 4]], dtype=array_api_strict.int64)
    """
    if xp is None:
        xp = _compat.array_namespace(x1, x2)

    if side not in {"left", "right"}:
        message = "`side` must be either 'left' or 'right'."
        raise ValueError(message)

    xp_default_int = default_dtype(xp, kind="integral")
    x2_0d = x2.ndim == 0
    x1_1d = x1.ndim <= 1

    if x1_1d or _compat.is_torch_namespace(xp):
        x2 = xp.reshape(x2, ()) if (x2_0d and x1_1d) else x2
        out = xp.searchsorted(x1, x2, side=side)
        return xp.astype(out, xp_default_int, copy=False)

    return _agnostic._searching.searchsorted(x1, x2, side=side, xp=xp)
