"""Delegation layer for manipulation functions."""

from collections.abc import Sequence
from typing import Literal

from . import _agnostic
from ._lib import _compat, _helpers
from ._lib._typing import Array, ArrayNamespace

__all__ = ["atleast_nd", "broadcast_shapes", "expand_dims", "pad"]


def atleast_nd(x: Array, /, *, ndim: int, xp: ArrayNamespace | None = None) -> Array:
    """
    Recursively expand the dimension of an array to at least `ndim`.

    Parameters
    ----------
    x : array
        Input array.
    ndim : int
        The minimum number of dimensions for the result.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        An array with ``res.ndim`` >= `ndim`.
        If ``x.ndim`` >= `ndim`, `x` is returned.
        If ``x.ndim`` < `ndim`, `x` is expanded by prepending new axes
        until ``res.ndim`` equals `ndim`.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> x = xp.asarray([1])
    >>> xpx.atleast_nd(x, ndim=3, xp=xp)
    Array([[[1]]], dtype=array_api_strict.int64)

    >>> x = xp.asarray([[[1, 2],
    ...                  [3, 4]]])
    >>> xpx.atleast_nd(x, ndim=1, xp=xp) is x
    True
    """
    if xp is None:
        xp = _compat.array_namespace(x)

    if 1 <= ndim <= 2 and (
        _compat.is_numpy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return getattr(xp, f"atleast_{ndim}d")(x)

    return _agnostic._manipulation.atleast_nd(x, ndim=ndim, xp=xp)


@_helpers.deprecated(
    "`xpx.broadcast_shapes` is deprecated and will be removed in v1.0.0. "
    "`xp.broadcast_shapes` exists in the standard as of v2025.12."
)
def broadcast_shapes(
    *shapes: tuple[float | None, ...], xp: ArrayNamespace | None = None
) -> tuple[int | None, ...]:
    """
    Compute the shape of the broadcasted arrays.

    .. deprecated:: 0.11.0
        :func:`broadcast_shapes` is deprecated and will be removed in v1.0.0.
        :func:`array_api.broadcast_shapes` exists in the standard as of v2025.12.

    Duplicates :func:`numpy.broadcast_shapes`, with additional support for
    None and NaN sizes.

    Parameters
    ----------
    *shapes : tuple[int | None, ...]
        Shapes of the arrays to broadcast.
    xp : array_namespace, optional
        The standard-compatible namespace to use for native delegation.
        Default: use the array-agnostic implementation.

    Returns
    -------
    tuple[int | None, ...]
        The shape of the broadcasted arrays.

    See Also
    --------
    numpy.broadcast_shapes : Equivalent NumPy function.
    array_api.broadcast_arrays : Function to broadcast actual arrays.

    Notes
    -----
    This function accepts the Array API's ``None`` for unknown sizes,
    as well as Dask's non-standard ``math.nan``.
    Regardless of input, the output always contains ``None`` for unknown sizes.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> xpx.broadcast_shapes((2, 3), (2, 1))
    (2, 3)
    >>> xpx.broadcast_shapes((4, 2, 3), (2, 1), (1, 3))
    (4, 2, 3)
    """
    if (
        xp is not None
        and all(isinstance(size, int) for shape in shapes for size in shape)
        and (
            _compat.is_numpy_namespace(xp)
            or _compat.is_cupy_namespace(xp)
            or _compat.is_jax_namespace(xp)
            or _compat.is_torch_namespace(xp)
        )
    ):
        return xp.broadcast_shapes(*shapes)

    return _agnostic._manipulation.broadcast_shapes(*shapes)


@_helpers.deprecated(
    "`xpx.expand_dims` is deprecated and will be removed in v1.0.0. "
    "`xp.expand_dims` with support for a tuple of ints in `axis` "
    "exists in the standard as of v2025.12."
)
def expand_dims(
    a: Array, /, *, axis: int | tuple[int, ...] = (0,), xp: ArrayNamespace | None = None
) -> Array:
    """
    Expand the shape of an array.

    .. deprecated:: 0.11.0
        :func:`expand_dims` is deprecated and will be removed in v1.0.0.
        :func:`array_api.expand_dims` with support for a tuple of ints in `axis`
        exists in the standard as of v2025.12.

    Insert (a) new axis/axes that will appear at the position(s) specified by
    `axis` in the expanded array shape.

    Parameters
    ----------
    a : array
        Array to have its shape expanded.
    axis : int or tuple of ints, optional
        Position(s) in the expanded axes where the new axis (or axes) is/are placed.
        If multiple positions are provided, they should be unique (note that a position
        given by a positive index could also be referred to by a negative index -
        that will also result in an error).
        Default: ``(0,)``.
    xp : array_namespace, optional
        The standard-compatible namespace for `a`. Default: infer.

    Returns
    -------
    array
        `a` with an expanded shape.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> x = xp.asarray([1, 2])
    >>> x.shape
    (2,)

    The following is equivalent to ``x[xp.newaxis, :]`` or ``x[xp.newaxis]``:

    >>> y = xpx.expand_dims(x, axis=0, xp=xp)
    >>> y
    Array([[1, 2]], dtype=array_api_strict.int64)
    >>> y.shape
    (1, 2)

    The following is equivalent to ``x[:, xp.newaxis]``:

    >>> y = xpx.expand_dims(x, axis=1, xp=xp)
    >>> y
    Array([[1],
           [2]], dtype=array_api_strict.int64)
    >>> y.shape
    (2, 1)

    ``axis`` may also be a tuple:

    >>> y = xpx.expand_dims(x, axis=(0, 1), xp=xp)
    >>> y
    Array([[[1, 2]]], dtype=array_api_strict.int64)

    >>> y = xpx.expand_dims(x, axis=(2, 0), xp=xp)
    >>> y
    Array([[[1],
            [2]]], dtype=array_api_strict.int64)
    """
    if xp is None:
        xp = _compat.array_namespace(a)

    if not isinstance(axis, tuple):
        axis = (axis,)
    ndim = a.ndim + len(axis)
    if axis != () and (min(axis) < -ndim or max(axis) >= ndim):
        err_msg = (
            f"a provided axis position is out of bounds for array of dimension {a.ndim}"
        )
        raise IndexError(err_msg)
    axis = tuple(dim % ndim for dim in axis)
    if len(set(axis)) != len(axis):
        err_msg = "Duplicate dimensions specified in `axis`."
        raise ValueError(err_msg)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.expand_dims(a, axis=axis)

    return _agnostic._manipulation.expand_dims(a, axis=axis, xp=xp)


def pad(
    x: Array,
    pad_width: int | tuple[int, int] | Sequence[tuple[int, int]],
    mode: Literal["constant"] = "constant",
    *,
    constant_values: complex = 0,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Pad the input array.

    Parameters
    ----------
    x : array
        Input array.
    pad_width : int or tuple of ints or sequence of pairs of ints
        Pad the input array with this many elements from each side.
        If a sequence of tuples, ``[(before_0, after_0), ... (before_N, after_N)]``,
        each pair applies to the corresponding axis of ``x``.
        A single tuple, ``(before, after)``, is equivalent to a list of ``x.ndim``
        copies of this tuple.
    mode : str, optional
        Only "constant" mode is currently supported, which pads with
        the value passed to `constant_values`.
    constant_values : python scalar, optional
        Use this value to pad the input. Default is zero.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        The input array,
        padded with ``pad_width`` elements equal to ``constant_values``.
    """
    xp = _compat.array_namespace(x) if xp is None else xp

    if mode != "constant":
        msg = "Only `'constant'` mode is currently supported"
        raise NotImplementedError(msg)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_pydata_sparse_namespace(xp)
    ):
        return xp.pad(x, pad_width, mode, constant_values=constant_values)

    if _compat.is_torch_namespace(xp):
        # normalize `pad_width` on the host rather than through a tensor as done in
        # `torch/_numpy`'s implementation (avoids device transfers)
        pad_width_seq = _helpers.normalize_pad_width(pad_width, x.ndim)
        # torch.nn.functional.pad counts dimensions from the last one
        flat_pad_width = [w for pair in reversed(pad_width_seq) for w in pair]
        return xp.nn.functional.pad(x, tuple(flat_pad_width), value=constant_values)

    return _agnostic._manipulation.pad(
        x, pad_width, constant_values=constant_values, xp=xp
    )
