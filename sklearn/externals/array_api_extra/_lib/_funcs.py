"""Array-agnostic implementations for the public API."""

import math
import warnings
from collections.abc import Callable, Sequence
from types import ModuleType, NoneType
from typing import Literal, cast, overload

from ._at import at
from ._utils import _compat, _helpers
from ._utils._compat import (
    array_namespace,
    is_dask_namespace,
    is_jax_array,
)
from ._utils._helpers import (
    asarrays,
    capabilities,
    eager_shape,
    meta_namespace,
    ndindex,
    normalize_pad_width,
)
from ._utils._typing import Array, Device, DType

__all__ = [
    "angle",
    "apply_where",
    "argpartition",
    "atleast_nd",
    "broadcast_shapes",
    "cov",
    "create_diagonal",
    "default_dtype",
    "diag_indices",
    "expand_dims",
    "isclose",
    "isin",
    "kron",
    "nan_to_num",
    "nanmax",
    "nanmin",
    "nunique",
    "one_hot",
    "pad",
    "partition",
    "searchsorted",
    "setdiff1d",
    "sinc",
    "tril_indices",
    "triu_indices",
    "union1d",
    "unravel_index",
]


@overload
def apply_where(  # numpydoc ignore=GL08
    cond: Array,
    args: Array | tuple[Array, ...],
    f1: Callable[..., Array],
    f2: Callable[..., Array],
    /,
    *,
    kwargs: dict[str, Array] | None = None,
    xp: ModuleType | None = None,
) -> Array: ...


@overload
def apply_where(  # numpydoc ignore=GL08
    cond: Array,
    args: Array | tuple[Array, ...],
    f1: Callable[..., Array],
    /,
    *,
    fill_value: Array | complex,
    kwargs: dict[str, Array] | None = None,
    xp: ModuleType | None = None,
) -> Array: ...


def apply_where(  # numpydoc ignore=PR01,PR02
    cond: Array,
    args: Array | tuple[Array, ...],
    f1: Callable[..., Array],
    f2: Callable[..., Array] | None = None,
    /,
    *,
    fill_value: Array | complex | None = None,
    kwargs: dict[str, Array] | None = None,
    xp: ModuleType | None = None,
) -> Array:
    """
    Run one of two elementwise functions depending on a condition.

    Equivalent to ``f1(*args) if cond else fill_value`` performed elementwise
    when `fill_value` is defined, otherwise to ``f1(*args) if cond else f2(*args)``.

    Parameters
    ----------
    cond : array
        The condition, expressed as a boolean array.
    args : Array or tuple of Arrays
        Argument(s) to `f1` (and `f2`). Must be broadcastable with `cond`.
    f1 : callable
        Elementwise function of `args`, returning a single array.
        Where `cond` is True, output will be ``f1(arg0[cond], arg1[cond], ...)``.
    f2 : callable, optional
        Elementwise function of `args`, returning a single array.
        Where `cond` is False, output will be ``f2(arg0[cond], arg1[cond], ...)``.
        Mutually exclusive with `fill_value`.
    fill_value : Array or scalar, optional
        If provided, value with which to fill output array where `cond` is False.
        It does not need to be scalar; it needs however to be broadcastable with
        `cond` and `args`.
        Mutually exclusive with `f2`. You must provide one or the other.
    kwargs : dict of str : Array pairs
        Keyword argument(s) to `f1` (and `f2`). Values must be broadcastable with
        `cond`.
    xp : array_namespace, optional
        The standard-compatible namespace for `cond` and `args`. Default: infer.

    Returns
    -------
    Array
        An array with elements from the output of `f1` where `cond` is True and either
        the output of `f2` or `fill_value` where `cond` is False. The returned array has
        data type determined by type promotion rules between the output of `f1` and
        either `fill_value` or the output of `f2`.

    Notes
    -----
    ``xp.where(cond, f1(*args), f2(*args))`` requires explicitly evaluating `f1` even
    when `cond` is False, and `f2` when cond is True. This function evaluates each
    function only for their matching condition, if the backend allows for it.

    On Dask, `f1` and `f2` are applied to the individual chunks and should use functions
    from the namespace of the chunks.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> a = xp.asarray([5, 4, 3])
    >>> b = xp.asarray([0, 2, 2])
    >>> def f(a, b):
    ...     return a // b
    >>> xpx.apply_where(b != 0, (a, b), f, fill_value=xp.nan)
    array([ nan,  2., 1.])
    """
    # Parse and normalize arguments
    if (f2 is None) == (fill_value is None):
        msg = "Exactly one of `fill_value` or `f2` must be given."
        raise TypeError(msg)
    args_ = list(args) if isinstance(args, tuple) else [args]
    del args

    kwargs_ = {} if kwargs is None else kwargs
    kwkeys = list(kwargs_.keys())
    args_ = [*args_, *kwargs_.values()]
    del kwargs

    xp = array_namespace(cond, fill_value, *args_) if xp is None else xp

    if isinstance(fill_value, int | float | complex | NoneType):
        cond, *args_ = xp.broadcast_arrays(cond, *args_)
    else:
        cond, fill_value, *args_ = xp.broadcast_arrays(cond, fill_value, *args_)

    if is_dask_namespace(xp):
        meta_xp = meta_namespace(cond, fill_value, *args_, xp=xp)
        # map_blocks doesn't descend into tuples of Arrays
        return xp.map_blocks(
            _apply_where, cond, f1, f2, fill_value, *args_, kwkeys=kwkeys, xp=meta_xp
        )

    return _apply_where(cond, f1, f2, fill_value, *args_, kwkeys=kwkeys, xp=xp)


def _apply_where(  # numpydoc ignore=PR01,RT01
    cond: Array,
    f1: Callable[..., Array],
    f2: Callable[..., Array] | None,
    fill_value: Array | int | float | complex | bool | None,
    *args: Array,
    kwkeys: list[str],
    xp: ModuleType,
) -> Array:
    """Helper of `apply_where`. On Dask, this runs on a single chunk."""

    nargs = len(args) - len(kwkeys)
    kwargs = dict(zip(kwkeys, args[nargs:], strict=True))
    args = args[:nargs]

    if not capabilities(xp, device=_compat.device(cond))["boolean indexing"]:
        # jax.jit does not support assignment by boolean mask
        return xp.where(
            cond,
            f1(*args, **kwargs),
            f2(*args, **kwargs) if f2 is not None else fill_value,
        )

    temp1 = f1(
        *(arr[cond] for arr in args), **{key: val[cond] for key, val in kwargs.items()}
    )

    if f2 is None:
        dtype = xp.result_type(temp1, fill_value)
        if isinstance(fill_value, int | float | complex):
            out = xp.full_like(cond, dtype=dtype, fill_value=fill_value)
        else:
            out = xp.astype(fill_value, dtype, copy=True)
    else:
        ncond = ~cond
        temp2 = f2(
            *(arr[ncond] for arr in args),
            **{key: val[ncond] for key, val in kwargs.items()},
        )
        dtype = xp.result_type(temp1, temp2)
        out = xp.empty_like(cond, dtype=dtype)
        out = at(out, ncond).set(temp2)

    return at(out, cond).set(temp1)


def atleast_nd(x: Array, /, *, ndim: int, xp: ModuleType) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""

    if x.ndim < ndim:
        x = xp.expand_dims(x, axis=0)
        x = atleast_nd(x, ndim=ndim, xp=xp)
    return x


# `float` in signature to accept `math.nan` for Dask.
# `int`s are still accepted as `float` is a superclass of `int` in typing
def broadcast_shapes(  # numpydoc ignore=PR01,RT01
    *shapes: tuple[float | None, ...],
) -> tuple[int | None, ...]:
    """See docstring in array_api_extra._delegation."""
    if not shapes:
        return ()  # Match NumPy output

    ndim = max(len(shape) for shape in shapes)
    out: list[int | None] = []
    for axis in range(-ndim, 0):
        sizes = {shape[axis] for shape in shapes if axis >= -len(shape)}
        # Dask uses NaN for unknown shape, which predates the Array API spec for None
        none_size = None in sizes or math.nan in sizes  # noqa: PLW0177
        sizes -= {1, None, math.nan}
        if len(sizes) > 1:
            msg = (
                "shape mismatch: objects cannot be broadcast to a single shape: "
                f"{shapes}."
            )
            raise ValueError(msg)
        out.append(None if none_size else cast(int, sizes.pop()) if sizes else 1)

    return tuple(out)


def cov(m: Array, /, *, xp: ModuleType) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    m = xp.asarray(m, copy=True)
    dtype = (
        xp.float64 if xp.isdtype(m.dtype, "integral") else xp.result_type(m, xp.float64)
    )

    m = atleast_nd(m, ndim=2, xp=xp)
    m = xp.astype(m, dtype)

    avg = xp.mean(m, axis=-1, keepdims=True)

    m_shape = eager_shape(m)
    fact = m_shape[-1] - 1

    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice", RuntimeWarning, stacklevel=2)
        fact = 0

    m -= avg
    m_transpose = xp.matrix_transpose(m)
    if xp.isdtype(m_transpose.dtype, "complex floating"):
        m_transpose = xp.conj(m_transpose)
    c = xp.matmul(m, m_transpose)
    c /= fact
    axes = tuple(axis for axis, length in enumerate(c.shape) if length == 1)
    return xp.squeeze(c, axis=axes)


def one_hot(
    x: Array,
    /,
    num_classes: int,
    *,
    xp: ModuleType,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""
    # TODO: Benchmark whether this is faster on the NumPy backend:
    # if is_numpy_array(x):
    #     out = xp.zeros((x.size, num_classes), dtype=dtype)
    #     out[xp.arange(x.size), xp.reshape(x, (-1,))] = 1
    #     return xp.reshape(out, (*x.shape, num_classes))
    range_num_classes = xp.arange(num_classes, dtype=x.dtype, device=_compat.device(x))
    return x[..., xp.newaxis] == range_num_classes


def create_diagonal(
    x: Array, /, *, offset: int = 0, xp: ModuleType
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    x_shape = eager_shape(x)
    batch_dims = x_shape[:-1]
    n = x_shape[-1] + abs(offset)
    diag = xp.zeros((*batch_dims, n**2), dtype=x.dtype, device=_compat.device(x))

    target_slice = slice(
        offset if offset >= 0 else abs(offset) * n,
        min(n * (n - offset), diag.shape[-1]),
        n + 1,
    )
    for index in ndindex(*batch_dims):
        diag = at(diag)[(*index, target_slice)].set(x[(*index, slice(None))])
    return xp.reshape(diag, (*batch_dims, n, n))


def diag_indices(
    n: int, /, *, ndim: int, device: Device | None, xp: ModuleType
) -> tuple[Array, ...]:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    idx = xp.arange(n, device=device)
    return (idx,) * ndim


def _tri_indices(
    n: int,
    *,
    offset: int,
    m: int | None,
    upper: bool,
    device: Device | None,
    xp: ModuleType,
) -> tuple[Array, Array]:  # numpydoc ignore=PR01,RT01
    """Shared implementation for `tril_indices` and `triu_indices`."""
    cols = n if m is None else m
    rows = xp.arange(n, device=device)[:, xp.newaxis]
    cols_a = xp.arange(cols, device=device)[xp.newaxis, :]
    delta = cols_a - rows
    mask = delta >= offset if upper else delta <= offset
    r, c = xp.nonzero(mask)
    return (r, c)


def tril_indices(
    n: int,
    /,
    *,
    offset: int,
    m: int | None,
    device: Device | None,
    xp: ModuleType,
) -> tuple[Array, Array]:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    return _tri_indices(n, offset=offset, m=m, upper=False, device=device, xp=xp)


def triu_indices(
    n: int,
    /,
    *,
    offset: int,
    m: int | None,
    device: Device | None,
    xp: ModuleType,
) -> tuple[Array, Array]:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    return _tri_indices(n, offset=offset, m=m, upper=True, device=device, xp=xp)


def default_dtype(
    xp: ModuleType,
    kind: Literal[
        "real floating", "complex floating", "integral", "indexing"
    ] = "real floating",
    *,
    device: Device | None = None,
) -> DType:
    """
    Return the default dtype for the given namespace and device.

    This is a convenience shorthand for
    ``xp.__array_namespace_info__().default_dtypes(device=device)[kind]``.

    Parameters
    ----------
    xp : array_namespace
        The standard-compatible namespace for which to get the default dtype.
    kind : {'real floating', 'complex floating', 'integral', 'indexing'}, optional
        The kind of dtype to return. Default is 'real floating'.
    device : Device, optional
        The device for which to get the default dtype. Default: current device.

    Returns
    -------
    dtype
        The default dtype for the given namespace, kind, and device.
    """
    dtypes = xp.__array_namespace_info__().default_dtypes(device=device)
    try:
        return dtypes[kind]
    except KeyError as e:
        domain = ("real floating", "complex floating", "integral", "indexing")
        assert set(dtypes) == set(domain), f"Non-compliant namespace: {dtypes}"
        msg = f"Unknown kind '{kind}'. Expected one of {domain}."
        raise ValueError(msg) from e


def expand_dims(a: Array, /, *, axis: tuple[int, ...] = (0,), xp: ModuleType) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    for i in sorted(axis):
        a = xp.expand_dims(a, axis=i)
    return a


def isclose(
    a: Array | complex,
    b: Array | complex,
    *,
    rtol: float = 1e-05,
    atol: float = 1e-08,
    equal_nan: bool = False,
    xp: ModuleType,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""
    a, b = asarrays(a, b, xp=xp)

    a_inexact = xp.isdtype(a.dtype, ("real floating", "complex floating"))
    b_inexact = xp.isdtype(b.dtype, ("real floating", "complex floating"))
    if a_inexact or b_inexact:
        # prevent warnings on NumPy and Dask on inf - inf
        mxp = meta_namespace(a, b, xp=xp)
        out = apply_where(
            xp.isinf(a) | xp.isinf(b),
            (a, b),
            lambda a, b: mxp.isinf(a) & mxp.isinf(b) & (mxp.sign(a) == mxp.sign(b)),  # pyright: ignore[reportUnknownArgumentType]
            # Note: inf <= inf is True!
            lambda a, b: mxp.abs(a - b) <= (atol + rtol * mxp.abs(b)),  # pyright: ignore[reportUnknownArgumentType]
            xp=xp,
        )
        if equal_nan:
            out = xp.where(xp.isnan(a) & xp.isnan(b), True, out)
        return out

    if xp.isdtype(a.dtype, "bool") or xp.isdtype(b.dtype, "bool"):
        if atol >= 1 or rtol >= 1:
            return xp.ones_like(a == b)
        return a == b

    # integer types
    atol = int(atol)
    if rtol == 0:
        return xp.abs(a - b) <= atol

    # Don't rely on OverflowError, as it is not guaranteed by the Array API.
    nrtol = int(1.0 / rtol)
    if nrtol > xp.iinfo(b.dtype).max:
        # rtol * max_int < 1, so it's inconsequential
        return xp.abs(a - b) <= atol
    return xp.abs(a - b) <= (atol + xp.abs(b) // nrtol)


def kron(
    a: Array,
    b: Array,
    /,
    *,
    xp: ModuleType,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in array_api_extra._delegation."""

    singletons = (1,) * (b.ndim - a.ndim)
    a = cast(Array, xp.broadcast_to(a, singletons + a.shape))

    nd_b, nd_a = b.ndim, a.ndim
    nd_max = max(nd_b, nd_a)
    if nd_a == 0 or nd_b == 0:
        return xp.multiply(a, b)

    a_shape = eager_shape(a)
    b_shape = eager_shape(b)

    # Equalise the shapes by prepending smaller one with 1s
    a_shape = (1,) * max(0, nd_b - nd_a) + a_shape
    b_shape = (1,) * max(0, nd_a - nd_b) + b_shape

    # Insert empty dimensions
    a_arr = expand_dims(a, axis=tuple(range(nd_b - nd_a)), xp=xp)
    b_arr = expand_dims(b, axis=tuple(range(nd_a - nd_b)), xp=xp)

    # Compute the product
    a_arr = expand_dims(a_arr, axis=tuple(range(1, nd_max * 2, 2)), xp=xp)
    b_arr = expand_dims(b_arr, axis=tuple(range(0, nd_max * 2, 2)), xp=xp)
    result = xp.multiply(a_arr, b_arr)

    # Reshape back and return
    res_shape = tuple(a_s * b_s for a_s, b_s in zip(a_shape, b_shape, strict=True))
    return xp.reshape(result, res_shape)


def nan_to_num(  # numpydoc ignore=PR01,RT01
    x: Array,
    /,
    fill_value: int | float = 0.0,
    *,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""

    def perform_replacements(  # numpydoc ignore=PR01,RT01
        x: Array,
        fill_value: int | float,
        xp: ModuleType,
    ) -> Array:
        """Internal function to perform the replacements."""
        x = xp.where(xp.isnan(x), fill_value, x)

        # convert infinities to finite values
        finfo = xp.finfo(x.dtype)
        idx_posinf = xp.isinf(x) & ~xp.signbit(x)
        idx_neginf = xp.isinf(x) & xp.signbit(x)
        x = xp.where(idx_posinf, finfo.max, x)
        return xp.where(idx_neginf, finfo.min, x)

    if xp.isdtype(x.dtype, "complex floating"):
        return perform_replacements(
            xp.real(x),
            fill_value,
            xp,
        ) + 1j * perform_replacements(
            xp.imag(x),
            fill_value,
            xp,
        )

    if xp.isdtype(x.dtype, "numeric"):
        return perform_replacements(x, fill_value, xp)

    return x


def nunique(x: Array, /, *, xp: ModuleType | None = None) -> Array:
    """
    Count the number of unique elements in an array.

    Compatible with JAX and Dask, whose laziness would be otherwise
    problematic.

    Parameters
    ----------
    x : Array
        Input array.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array: 0-dimensional integer array
        The number of unique elements in `x`. It can be lazy.
    """
    if xp is None:
        xp = array_namespace(x)

    if is_jax_array(x):
        # size= is JAX-specific
        # https://github.com/data-apis/array-api/issues/883
        _, counts = xp.unique_counts(x, size=_compat.size(x))
        return (counts > 0).sum()

    # There are 3 general use cases:
    # 1. backend has unique_counts and it returns an array with known shape
    # 2. backend has unique_counts and it returns a None-sized array;
    #    e.g. Dask, ndonnx
    # 3. backend does not have unique_counts; e.g. wrapped JAX
    if capabilities(xp, device=_compat.device(x))["data-dependent shapes"]:
        # xp has unique_counts; O(n) complexity
        _, counts = xp.unique_counts(x)
        n = _compat.size(counts)
        if n is None:
            return xp.sum(xp.ones_like(counts))
        return xp.asarray(n, device=_compat.device(x))

    # xp does not have unique_counts; O(n*logn) complexity
    x = xp.reshape(x, (-1,))
    x = xp.sort(x, stable=False)
    mask = x != xp.roll(x, -1)
    default_int = default_dtype(xp, "integral", device=_compat.device(x))
    return xp.maximum(
        # Special cases:
        # - array is size 0
        # - array has all elements equal to each other
        xp.astype(xp.any(~mask), default_int),
        xp.sum(xp.astype(mask, default_int)),
    )


def pad(
    x: Array,
    pad_width: int | tuple[int, int] | Sequence[tuple[int, int]],
    *,
    constant_values: complex = 0,
    xp: ModuleType,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""
    pad_width_seq = normalize_pad_width(pad_width, x.ndim)

    slices: list[slice] = []
    newshape: list[int] = []
    for ax, w_tpl in enumerate(pad_width_seq):
        if len(w_tpl) != 2:
            msg = f"expect a 2-tuple (before, after), got {w_tpl}."
            raise ValueError(msg)

        sh = eager_shape(x)[ax]

        if w_tpl[0] == 0 and w_tpl[1] == 0:
            sl = slice(None, None, None)
        else:
            stop: int | None
            start, stop = w_tpl
            stop = None if stop == 0 else -stop

            sl = slice(start, stop, None)
            sh += w_tpl[0] + w_tpl[1]

        newshape.append(sh)
        slices.append(sl)

    padded = xp.full(
        tuple(newshape),
        fill_value=constant_values,
        dtype=x.dtype,
        device=_compat.device(x),
    )
    return at(padded, tuple(slices)).set(x)


def searchsorted(
    x1: Array,
    x2: Array,
    /,
    *,
    side: Literal["left", "right"] = "left",
    xp: ModuleType,
) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""
    a = xp.full(x2.shape, 0, device=_compat.device(x1))

    if x1.shape[-1] == 0:
        return a

    n = xp.count_nonzero(~xp.isnan(x1), axis=-1, keepdims=True)
    b = xp.broadcast_to(n, x2.shape)

    compare = xp.less_equal if side == "left" else xp.less

    # while xp.any(b - a > 1):
    # refactored to for loop with ~log2(n) iterations for JAX JIT
    for _ in range(int(math.log2(x1.shape[-1])) + 1):  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]
        c = (a + b) // 2
        x0 = xp.take_along_axis(x1, c, axis=-1)
        j = compare(x2, x0)
        b = xp.where(j, c, b)
        a = xp.where(j, a, c)

    out = xp.where(compare(x2, xp.min(x1, axis=-1, keepdims=True)), 0, b)
    out = xp.where(xp.isnan(x2), x1.shape[-1], out) if side == "right" else out
    return xp.astype(out, default_dtype(xp, kind="integral"), copy=False)


def setdiff1d(
    x1: Array | complex,
    x2: Array | complex,
    /,
    *,
    assume_unique: bool = False,
    xp: ModuleType,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""

    # https://github.com/microsoft/pyright/issues/10103
    x1_, x2_ = asarrays(x1, x2, xp=xp)

    if assume_unique:
        x1_ = xp.reshape(x1_, (-1,))
        x2_ = xp.reshape(x2_, (-1,))
    else:
        x1_ = xp.unique_values(x1_)
        x2_ = xp.unique_values(x2_)

    return x1_[_helpers.in1d(x1_, x2_, assume_unique=True, invert=True, xp=xp)]


def sinc(x: Array, /, *, xp: ModuleType) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""

    # no scalars in `where` - array-api#807
    y = xp.pi * xp.where(
        xp.astype(x, xp.bool),
        x,
        xp.asarray(xp.finfo(x.dtype).eps, dtype=x.dtype, device=_compat.device(x)),
    )
    return xp.sin(y) / y


def partition(  # numpydoc ignore=PR01,RT01
    x: Array,
    kth: int,  # noqa: ARG001
    /,
    axis: int = -1,
    *,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""
    return xp.sort(x, axis=axis, stable=False)


def argpartition(  # numpydoc ignore=PR01,RT01
    x: Array,
    kth: int,  # noqa: ARG001
    /,
    axis: int = -1,
    *,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""
    return xp.argsort(x, axis=axis, stable=False)


def isin(  # numpydoc ignore=PR01,RT01
    a: Array,
    b: Array,
    /,
    *,
    assume_unique: bool = False,
    invert: bool = False,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""
    original_a_shape = a.shape
    a = xp.reshape(a, (-1,))
    b = xp.reshape(b, (-1,))
    return xp.reshape(
        _helpers.in1d(a, b, assume_unique=assume_unique, invert=invert, xp=xp),
        original_a_shape,
    )


def union1d(a: Array, b: Array, /, *, xp: ModuleType) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""
    a = xp.reshape(a, (-1,))
    b = xp.reshape(b, (-1,))
    # XXX: `sparse` returns NumPy arrays from `unique_values`
    return xp.asarray(xp.unique_values(xp.concat([a, b])))


def angle(z: Array, /, *, deg: bool = False, xp: ModuleType | None = None) -> Array:
    """
    Return the angle of the complex argument.

    Parameters
    ----------
    z : Array
        Input array.
    deg : bool, optional
        Return angle in degrees if True, radians if False (default).
    xp : array_namespace, optional
        The standard-compatible namespace for `z`. Default: infer.

    Returns
    -------
    array
        The counterclockwise angle from the positive real axis on the complex
        plane in the range ``(-pi, pi]``.

    Notes
    -----
    Real input ``x`` is interpreted as ``x + 0j``.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> xpx.angle(xp.asarray([1.0, 1.0j, 1 + 1j]), xp=xp)
    Array([0.        , 1.57079633, 0.78539816], dtype=array_api_strict.float64)
    >>> xpx.angle(xp.asarray([1.0, 1.0j, 1 + 1j]), deg=True, xp=xp)
    Array([ 0., 90., 45.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = array_namespace(z)
    if xp.isdtype(z.dtype, "complex floating"):
        zimag = xp.imag(z)
        zreal = xp.real(z)
    else:
        if not xp.isdtype(z.dtype, "real floating"):
            z = xp.astype(z, default_dtype(xp, device=_compat.device(z)))
        zimag = xp.zeros_like(z)
        zreal = z
    a = xp.atan2(zimag, zreal)
    if deg:
        a = a * 180 / xp.pi
    return a


def unravel_index(indices: Array, shape: tuple[int, ...], /) -> tuple[Array, ...]:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._delegation.py`."""
    coords: list[Array] = []
    for dim in reversed(shape):
        coords.append(indices % dim)
        indices = indices // dim
    return tuple(reversed(coords))


def nanmin(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""
    mask = xp.isnan(a)
    device_a = _compat.device(a)
    x = xp.min(
        xp.where(mask, xp.asarray(+xp.inf, dtype=a.dtype, device=device_a), a),
        axis=axis,
    )
    # Replace Infs from all NaN slices with NaN again
    mask = xp.all(mask, axis=axis)
    if xp.any(mask):
        x = xp.where(mask, xp.asarray(xp.nan, dtype=x.dtype, device=device_a), x)
    return x


def nanmax(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ModuleType,
) -> Array:
    """See docstring in `array_api_extra._delegation.py`."""
    mask = xp.isnan(a)
    device_a = _compat.device(a)
    x = xp.max(
        xp.where(mask, xp.asarray(-xp.inf, dtype=a.dtype, device=device_a), a),
        axis=axis,
    )
    # Replace Infs from all NaN slices with NaN again
    mask = xp.all(mask, axis=axis)
    if xp.any(mask):
        x = xp.where(mask, xp.asarray(xp.nan, dtype=x.dtype, device=device_a), x)
    return x
