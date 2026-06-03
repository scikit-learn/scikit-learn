"""Delegation to existing implementations for Public API Functions."""

from collections.abc import Sequence
from types import ModuleType
from typing import Literal

from ._lib import _funcs
from ._lib._utils._compat import (
    array_namespace,
    is_cupy_namespace,
    is_dask_namespace,
    is_jax_namespace,
    is_numpy_namespace,
    is_pydata_sparse_namespace,
    is_torch_namespace,
)
from ._lib._utils._compat import device as get_device
from ._lib._utils._helpers import asarrays
from ._lib._utils._typing import Array, DType

__all__ = ["isclose", "nan_to_num", "one_hot", "pad"]


def isclose(
    a: Array | complex,
    b: Array | complex,
    *,
    rtol: float = 1e-05,
    atol: float = 1e-08,
    equal_nan: bool = False,
    xp: ModuleType | None = None,
) -> Array:
    """
    Return a boolean array where two arrays are element-wise equal within a tolerance.

    The tolerance values are positive, typically very small numbers. The relative
    difference ``(rtol * abs(b))`` and the absolute difference `atol` are added together
    to compare against the absolute difference between `a` and `b`.

    NaNs are treated as equal if they are in the same place and if ``equal_nan=True``.
    Infs are treated as equal if they are in the same place and of the same sign in both
    arrays.

    Parameters
    ----------
    a, b : Array | int | float | complex | bool
        Input objects to compare. At least one must be an array.
    rtol : array_like, optional
        The relative tolerance parameter (see Notes).
    atol : array_like, optional
        The absolute tolerance parameter (see Notes).
    equal_nan : bool, optional
        Whether to compare NaN's as equal. If True, NaN's in `a` will be considered
        equal to NaN's in `b` in the output array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a` and `b`. Default: infer.

    Returns
    -------
    Array
        A boolean array of shape broadcasted from `a` and `b`, containing ``True`` where
        `a` is close to `b`, and ``False`` otherwise.

    Warnings
    --------
    The default `atol` is not appropriate for comparing numbers with magnitudes much
    smaller than one (see notes).

    See Also
    --------
    math.isclose : Similar function in stdlib for Python scalars.

    Notes
    -----
    For finite values, `isclose` uses the following equation to test whether two
    floating point values are equivalent::

        absolute(a - b) <= (atol + rtol * absolute(b))

    Unlike the built-in `math.isclose`,
    the above equation is not symmetric in `a` and `b`,
    so that ``isclose(a, b)`` might be different from ``isclose(b, a)`` in some rare
    cases.

    The default value of `atol` is not appropriate when the reference value `b` has
    magnitude smaller than one. For example, it is unlikely that ``a = 1e-9`` and
    ``b = 2e-9`` should be considered "close", yet ``isclose(1e-9, 2e-9)`` is ``True``
    with default settings. Be sure to select `atol` for the use case at hand, especially
    for defining the threshold below which a non-zero value in `a` will be considered
    "close" to a very small or zero value in `b`.

    The comparison of `a` and `b` uses standard broadcasting, which means that `a` and
    `b` need not have the same shape in order for ``isclose(a, b)`` to evaluate to
    ``True``.

    `isclose` is not defined for non-numeric data types.
    ``bool`` is considered a numeric data-type for this purpose.
    """
    xp = array_namespace(a, b) if xp is None else xp

    if (
        is_numpy_namespace(xp)
        or is_cupy_namespace(xp)
        or is_dask_namespace(xp)
        or is_jax_namespace(xp)
    ):
        return xp.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    if is_torch_namespace(xp):
        a, b = asarrays(a, b, xp=xp)  # Array API 2024.12 support
        return xp.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    return _funcs.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan, xp=xp)


def nan_to_num(
    x: Array | float | complex,
    /,
    *,
    fill_value: int | float = 0.0,
    xp: ModuleType | None = None,
) -> Array:
    """
    Replace NaN with zero and infinity with large finite numbers (default behaviour).

    If `x` is inexact, NaN is replaced by zero or by the user defined value in the
    `fill_value` keyword, infinity is replaced by the largest finite floating
    point value representable by ``x.dtype``, and -infinity is replaced by the
    most negative finite floating point value representable by ``x.dtype``.

    For complex dtypes, the above is applied to each of the real and
    imaginary components of `x` separately.

    Parameters
    ----------
    x : array | float | complex
        Input data.
    fill_value : int | float, optional
        Value to be used to fill NaN values. If no value is passed
        then NaN values will be replaced with 0.0.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        `x`, with the non-finite values replaced.

    See Also
    --------
    array_api.isnan : Shows which elements are Not a Number (NaN).

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> xpx.nan_to_num(xp.inf)
    1.7976931348623157e+308
    >>> xpx.nan_to_num(-xp.inf)
    -1.7976931348623157e+308
    >>> xpx.nan_to_num(xp.nan)
    0.0
    >>> x = xp.asarray([xp.inf, -xp.inf, xp.nan, -128, 128])
    >>> xpx.nan_to_num(x)
    array([ 1.79769313e+308, -1.79769313e+308,  0.00000000e+000, # may vary
           -1.28000000e+002,  1.28000000e+002])
    >>> y = xp.asarray([complex(xp.inf, xp.nan), xp.nan, complex(xp.nan, xp.inf)])
    array([  1.79769313e+308,  -1.79769313e+308,   0.00000000e+000, # may vary
         -1.28000000e+002,   1.28000000e+002])
    >>> xpx.nan_to_num(y)
    array([  1.79769313e+308 +0.00000000e+000j, # may vary
             0.00000000e+000 +0.00000000e+000j,
             0.00000000e+000 +1.79769313e+308j])
    """
    if isinstance(fill_value, complex):
        msg = "Complex fill values are not supported."
        raise TypeError(msg)

    xp = array_namespace(x) if xp is None else xp

    # for scalars we want to output an array
    y = xp.asarray(x)

    if (
        is_cupy_namespace(xp)
        or is_jax_namespace(xp)
        or is_numpy_namespace(xp)
        or is_torch_namespace(xp)
    ):
        return xp.nan_to_num(y, nan=fill_value)

    return _funcs.nan_to_num(y, fill_value=fill_value, xp=xp)


def one_hot(
    x: Array,
    /,
    num_classes: int,
    *,
    dtype: DType | None = None,
    axis: int = -1,
    xp: ModuleType | None = None,
) -> Array:
    """
    One-hot encode the given indices.

    Each index in the input `x` is encoded as a vector of zeros of length `num_classes`
    with the element at the given index set to one.

    Parameters
    ----------
    x : array
        An array with integral dtype whose values are between `0` and `num_classes - 1`.
    num_classes : int
        Number of classes in the one-hot dimension.
    dtype : DType, optional
        The dtype of the return value.  Defaults to the default float dtype (usually
        float64).
    axis : int, optional
        Position in the expanded axes where the new axis is placed. Default: -1.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        An array having the same shape as `x` except for a new axis at the position
        given by `axis` having size `num_classes`.  If `axis` is unspecified, it
        defaults to -1, which appends a new axis.

        If ``x < 0`` or ``x >= num_classes``, then the result is undefined, may raise
        an exception, or may even cause a bad state.  `x` is not checked.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> xpx.one_hot(xp.asarray([1, 2, 0]), 3)
    Array([[0., 1., 0.],
          [0., 0., 1.],
          [1., 0., 0.]], dtype=array_api_strict.float64)
    """
    # Validate inputs.
    if xp is None:
        xp = array_namespace(x)
    if not xp.isdtype(x.dtype, "integral"):
        msg = "x must have an integral dtype."
        raise TypeError(msg)
    if dtype is None:
        dtype = _funcs.default_dtype(xp, device=get_device(x))
    # Delegate where possible.
    if is_jax_namespace(xp):
        from jax.nn import one_hot as jax_one_hot

        return jax_one_hot(x, num_classes, dtype=dtype, axis=axis)
    if is_torch_namespace(xp):
        from torch.nn.functional import one_hot as torch_one_hot

        x = xp.astype(x, xp.int64)  # PyTorch only supports int64 here.
        try:
            out = torch_one_hot(x, num_classes)
        except RuntimeError as e:
            raise IndexError from e
    else:
        out = _funcs.one_hot(x, num_classes, xp=xp)
    out = xp.astype(out, dtype, copy=False)
    if axis != -1:
        out = xp.moveaxis(out, -1, axis)
    return out


def pad(
    x: Array,
    pad_width: int | tuple[int, int] | Sequence[tuple[int, int]],
    mode: Literal["constant"] = "constant",
    *,
    constant_values: complex = 0,
    xp: ModuleType | None = None,
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
    xp = array_namespace(x) if xp is None else xp

    if mode != "constant":
        msg = "Only `'constant'` mode is currently supported"
        raise NotImplementedError(msg)

    if (
        is_numpy_namespace(xp)
        or is_cupy_namespace(xp)
        or is_jax_namespace(xp)
        or is_pydata_sparse_namespace(xp)
    ):
        return xp.pad(x, pad_width, mode, constant_values=constant_values)

    # https://github.com/pytorch/pytorch/blob/cf76c05b4dc629ac989d1fb8e789d4fac04a095a/torch/_numpy/_funcs_impl.py#L2045-L2056
    if is_torch_namespace(xp):
        pad_width = xp.asarray(pad_width)
        pad_width = xp.broadcast_to(pad_width, (x.ndim, 2))
        pad_width = xp.flip(pad_width, axis=(0,)).flatten()
        return xp.nn.functional.pad(x, tuple(pad_width), value=constant_values)  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]

    return _funcs.pad(x, pad_width, constant_values=constant_values, xp=xp)
