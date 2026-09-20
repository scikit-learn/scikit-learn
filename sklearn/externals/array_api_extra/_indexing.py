"""Delegation layer for indexing functions."""

from . import _agnostic
from ._lib import _compat
from ._lib._typing import Array, ArrayNamespace, Device

__all__ = ["diag_indices", "tril_indices", "triu_indices", "unravel_index"]


def diag_indices(
    n: int, /, *, ndim: int = 2, device: Device | None = None, xp: ArrayNamespace
) -> tuple[Array, ...]:
    """
    Return the indices to access the main diagonal of an array.

    Equivalent to :func:`numpy.diag_indices`.

    Parameters
    ----------
    n : int
        The size of each dimension of the (hyper-)cube ``(n, n, ..., n)``
        that the returned indices index into.
    ndim : int, optional
        The number of dimensions. Default: ``2``.
    device : Device, optional
        The device on which to place the returned arrays. Default: current device.
    xp : array_namespace
        The standard-compatible namespace to create the indices in.

    Returns
    -------
    tuple of array
        1-D integer arrays of length ``n`` that together index
        the main diagonal of an array of shape ``(n,) * ndim``.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> rows, cols = xpx.diag_indices(3, xp=xp)
    >>> rows
    Array([0, 1, 2], dtype=array_api_strict.int64)
    >>> cols
    Array([0, 1, 2], dtype=array_api_strict.int64)
    """
    if n < 0:
        msg = f"`n` must be non-negative, got {n}"
        raise ValueError(msg)
    if ndim < 1:
        msg = f"`ndim` must be >= 1, got {ndim}"
        raise ValueError(msg)
    if device is None and (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.diag_indices(n, ndim=ndim)
    return _agnostic._indexing.diag_indices(n, ndim=ndim, device=device, xp=xp)


def tril_indices(
    n: int,
    /,
    *,
    offset: int = 0,
    m: int | None = None,
    device: Device | None = None,
    xp: ArrayNamespace,
) -> tuple[Array, Array]:
    """
    Return the indices of the lower triangle of an ``(n, m)`` array.

    Equivalent to :func:`numpy.tril_indices` with parameter ``k`` renamed to
    ``offset`` to match :func:`array_api.linalg.diagonal`'s naming.

    Parameters
    ----------
    n : int
        The row dimension of the array.
    offset : int, optional
        Diagonal offset; ``0`` (default) is the main diagonal. Corresponds
        to ``k`` in :func:`numpy.tril_indices`.
    m : int, optional
        The column dimension. If ``None`` (default), assumed equal to `n`.
    device : Device, optional
        The device on which to place the returned arrays. Default: current device.
    xp : array_namespace
        The standard-compatible namespace to create the indices in.

    Returns
    -------
    tuple of array
        Row and column indices ``(rows, cols)`` of the lower triangle of
        the ``(n, m)`` matrix, shifted by `offset`.

    Notes
    -----
    The generic fallback uses :func:`array_api.nonzero`, so namespaces without
    ``nonzero`` are not supported on that path.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> rows, cols = xpx.tril_indices(3, xp=xp)
    >>> rows
    Array([0, 1, 1, 2, 2, 2], dtype=array_api_strict.int64)
    >>> cols
    Array([0, 0, 1, 0, 1, 2], dtype=array_api_strict.int64)
    """
    if n < 0:
        msg = f"`n` must be non-negative, got {n}"
        raise ValueError(msg)
    if m is not None and m < 0:
        msg = f"`m` must be non-negative, got {m}"
        raise ValueError(msg)
    if device is None and (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_dask_namespace(xp)
    ):
        return xp.tril_indices(n, k=offset, m=m)
    if _compat.is_torch_namespace(xp):
        # `torch.tril_indices` returns a 2xN tensor, not a tuple, and
        # takes (row, col) rather than (n, *, m=None).
        cols = n if m is None else m
        idx = xp.tril_indices(n, cols, offset=offset, device=device)
        return (idx[0], idx[1])
    return _agnostic._indexing.tril_indices(n, offset=offset, m=m, device=device, xp=xp)


def triu_indices(
    n: int,
    /,
    *,
    offset: int = 0,
    m: int | None = None,
    device: Device | None = None,
    xp: ArrayNamespace,
) -> tuple[Array, Array]:
    """
    Return the indices of the upper triangle of an ``(n, m)`` array.

    Equivalent to :func:`numpy.triu_indices` with parameter ``k`` renamed to
    ``offset`` to match :func:`array_api.linalg.diagonal`'s naming.

    Parameters
    ----------
    n : int
        The row dimension of the array.
    offset : int, optional
        Diagonal offset; ``0`` (default) is the main diagonal. Corresponds
        to ``k`` in :func:`numpy.triu_indices`.
    m : int, optional
        The column dimension. If ``None`` (default), assumed equal to `n`.
    device : Device, optional
        The device on which to place the returned arrays. Default: current device.
    xp : array_namespace
        The standard-compatible namespace to create the indices in.

    Returns
    -------
    tuple of array
        Row and column indices ``(rows, cols)`` of the upper triangle of
        the ``(n, m)`` matrix, shifted by `offset`.

    Notes
    -----
    The generic fallback uses :func:`array_api.nonzero`, so namespaces without
    ``nonzero`` are not supported on that path.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> rows, cols = xpx.triu_indices(3, xp=xp)
    >>> rows
    Array([0, 0, 0, 1, 1, 2], dtype=array_api_strict.int64)
    >>> cols
    Array([0, 1, 2, 1, 2, 2], dtype=array_api_strict.int64)
    """
    if n < 0:
        msg = f"`n` must be non-negative, got {n}"
        raise ValueError(msg)
    if m is not None and m < 0:
        msg = f"`m` must be non-negative, got {m}"
        raise ValueError(msg)
    if device is None and (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_dask_namespace(xp)
    ):
        return xp.triu_indices(n, k=offset, m=m)
    if _compat.is_torch_namespace(xp):
        cols = n if m is None else m
        idx = xp.triu_indices(n, cols, offset=offset, device=device)
        return (idx[0], idx[1])
    return _agnostic._indexing.triu_indices(n, offset=offset, m=m, device=device, xp=xp)


def unravel_index(
    indices: Array,
    shape: tuple[int, ...],
    /,
    *,
    xp: ArrayNamespace | None = None,
) -> tuple[Array, ...]:
    """
    Convert a flat index or array of flat indices into a tuple of coordinate arrays.

    Parameters
    ----------
    indices : array
        An integer array whose elements are indices into the flattened version
        of an array of dimensions `shape`.

    shape : tuple of ints
        The shape to use for unraveling `indices`.

    xp : array_namespace, optional
        The standard-compatible namespace for `indices`. Default: infer.

    Returns
    -------
    tuple of array
        A tuple of unraveled indices. Each array in the tuple has the same shape
        as the `indices` array.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> xs, ys = xpx.unravel_index(xp.asarray([1, 2, 4, 5, 6, 8]), (4, 3))
    >>> xs
    Array([0, 0, 1, 1, 2, 2], dtype=array_api_strict.int64)
    >>> ys
    Array([1, 2, 1, 2, 0, 2], dtype=array_api_strict.int64)
    >>> [(int(x), int(y)) for x, y in zip(xs, ys)]
    [(0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 2)]
    >>> xs, ys = xpx.unravel_index(xp.arange(6), (2, 2))
    >>> [(int(x), int(y)) for x, y in zip(xs, ys)]
    [(0, 0), (0, 1), (1, 0), (1, 1), (0, 0), (0, 1)]
    """
    if xp is None:
        xp = _compat.array_namespace(indices)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return xp.unravel_index(indices, shape)

    return _agnostic._indexing.unravel_index(indices, shape)
