"""Delegation layer for sorting functions."""

from . import _agnostic
from ._lib import _compat, _helpers
from ._lib._typing import Array, ArrayNamespace

__all__ = ["argpartition", "partition"]


def partition(
    a: Array,
    kth: int,
    /,
    axis: int | None = -1,
    *,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Return a partitioned copy of an array.

    Creates a copy of the array and partially sorts it in such a way that the value
    of the element in k-th position is in the position it would be in a sorted array.
    In the output array, all elements smaller than the k-th element are located to
    the left of this element and all equal or greater are located to its right.
    The ordering of the elements in the two partitions on the either side of
    the k-th element in the output array is undefined.

    Parameters
    ----------
    a : Array
        Input array.
    kth : int
        Element index to partition by.
    axis : int, optional
        Axis along which to partition. The default is ``-1`` (the last axis).
        If ``None``, the flattened array is used.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    partitioned_array
        Array of the same type and shape as `a`.

    Notes
    -----
    If `xp` implements ``partition`` or an equivalent function
    (e.g. ``topk`` for torch), complexity will likely be O(n).
    If not, this function simply calls ``xp.sort`` and complexity is O(n log n).
    """
    # Validate inputs.
    if xp is None:
        xp = _compat.array_namespace(a)
    if a.ndim < 1:
        msg = "`a` must be at least 1-dimensional"
        raise TypeError(msg)
    if axis is None:
        return partition(xp.reshape(a, (-1,)), kth, axis=0, xp=xp)
    (size,) = _helpers.eager_shape(a, axis)
    if not (0 <= kth < size):
        msg = f"kth(={kth}) out of bounds [0 {size})"
        raise ValueError(msg)

    # Delegate where possible.
    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.partition(a, kth, axis=axis)

    # Use top-k when possible:
    if _compat.is_torch_namespace(xp):
        if not (axis == -1 or axis == a.ndim - 1):
            a = xp.transpose(a, axis, -1)

        out = xp.empty_like(a)
        ranks = xp.arange(a.shape[-1]).expand_as(a)

        split_value, indices = xp.kthvalue(a, kth + 1, keepdim=True)
        del indices  # indices won't be used => del ASAP to reduce peak memory usage

        # fill the left-side of the partition
        mask_src = a < split_value
        n_left = mask_src.sum(dim=-1, keepdim=True)
        mask_dest = ranks < n_left
        out[mask_dest] = a[mask_src]

        # fill the middle of the partition
        mask_src = a == split_value
        n_left += mask_src.sum(dim=-1, keepdim=True)
        mask_dest ^= ranks < n_left
        out[mask_dest] = a[mask_src]

        # fill the right-side of the partition
        mask_src = a > split_value
        mask_dest = ranks >= n_left
        out[mask_dest] = a[mask_src]

        if not (axis == -1 or axis == a.ndim - 1):
            out = xp.transpose(out, axis, -1)
        return out

    # Note: dask topk/argtopk sort the return values, so it's
    # not much more efficient than sorting everything when
    # kth is not small compared to x.size

    return _agnostic._sorting.partition(a, kth, axis=axis, xp=xp)


def argpartition(
    a: Array,
    kth: int,
    /,
    axis: int | None = -1,
    *,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Perform an indirect partition along the given axis.

    It returns an array of indices of the same shape as `a` that
    index data along the given axis in partitioned order.

    Parameters
    ----------
    a : Array
        Input array.
    kth : int
        Element index to partition by.
    axis : int, optional
        Axis along which to partition. The default is ``-1`` (the last axis).
        If ``None``, the flattened array is used.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    index_array
        Array of indices that partition `a` along the specified axis.

    Notes
    -----
    If `xp` implements ``argpartition`` or an equivalent function
    e.g. ``topk`` for torch), complexity will likely be O(n).
    If not, this function simply calls ``xp.argsort`` and complexity is O(n log n).
    """
    # Validate inputs.
    if xp is None:
        xp = _compat.array_namespace(a)
    if _compat.is_pydata_sparse_namespace(xp):
        msg = "Not implemented for sparse backend: no argsort"
        raise NotImplementedError(msg)
    if a.ndim < 1:
        msg = "`a` must be at least 1-dimensional"
        raise TypeError(msg)
    if axis is None:
        return argpartition(xp.reshape(a, (-1,)), kth, axis=0, xp=xp)
    (size,) = _helpers.eager_shape(a, axis)
    if not (0 <= kth < size):
        msg = f"kth(={kth}) out of bounds [0 {size})"
        raise ValueError(msg)

    # Delegate where possible.
    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.argpartition(a, kth, axis=axis)

    # Use top-k when possible:
    if _compat.is_torch_namespace(xp):
        # see `partition` above for commented details of those steps:
        if not (axis == -1 or axis == a.ndim - 1):
            a = xp.transpose(a, axis, -1)

        ranks = xp.arange(a.shape[-1]).expand_as(a)
        out = xp.empty_like(ranks)

        split_value, indices = xp.kthvalue(a, kth + 1, keepdim=True)
        del indices  # indices won't be used => del ASAP to reduce peak memory usage

        mask_src = a < split_value
        n_left = mask_src.sum(dim=-1, keepdim=True)
        mask_dest = ranks < n_left
        out[mask_dest] = ranks[mask_src]

        mask_src = a == split_value
        n_left += mask_src.sum(dim=-1, keepdim=True)
        mask_dest ^= ranks < n_left
        out[mask_dest] = ranks[mask_src]

        mask_src = a > split_value
        mask_dest = ranks >= n_left
        out[mask_dest] = ranks[mask_src]

        if not (axis == -1 or axis == a.ndim - 1):
            out = xp.transpose(out, axis, -1)
        return out

    # Note: dask topk/argtopk sort the return values, so it's
    # not much more efficient than sorting everything when
    # kth is not small compared to x.size

    return _agnostic._sorting.argpartition(a, kth, axis=axis, xp=xp)
