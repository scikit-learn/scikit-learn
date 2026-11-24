"""Delegation layer for element-wise functions."""

from . import _agnostic
from ._lib import _compat, _helpers
from ._lib._typing import Array, ArrayNamespace

__all__ = ["deg2rad", "isclose", "nan_to_num", "rad2deg", "sinc"]


def deg2rad(x: Array, /, *, xp: ArrayNamespace | None = None) -> Array:
    """
    Convert angles from degrees to radians.

    Parameters
    ----------
    x : array
        Input array in degrees. Must have an integral or floating-point dtype.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        The corresponding angles in radians. Integral inputs are converted to the
        default floating-point dtype.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> xpx.deg2rad(xp.asarray([0, 90, 180]), xp=xp)
    Array([0.        , 1.57079633, 3.14159265], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(x)
    if xp.isdtype(x.dtype, "integral"):
        dtype = _agnostic._inspection.default_dtype(xp, device=_compat.device(x))
        x = xp.astype(x, dtype)
    elif not xp.isdtype(x.dtype, ("real floating", "complex floating")):
        msg = "`x` must have an integral, real floating, or complex floating dtype."
        raise TypeError(msg)

    if _compat.is_jax_namespace(xp) or (
        not xp.isdtype(x.dtype, "complex floating")
        and (
            _compat.is_numpy_namespace(xp)
            or _compat.is_cupy_namespace(xp)
            or _compat.is_torch_namespace(xp)
            or _compat.is_dask_namespace(xp)
        )
    ):
        return xp.deg2rad(x)

    return _agnostic._elementwise.deg2rad(x, xp=xp)


def rad2deg(x: Array, /, *, xp: ArrayNamespace | None = None) -> Array:
    """
    Convert angles from radians to degrees.

    Parameters
    ----------
    x : array
        Input array in radians. Must have an integral or floating-point dtype.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        The corresponding angles in degrees. Integral inputs are converted to the
        default floating-point dtype.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> xpx.rad2deg(xp.asarray([0.0, xp.pi / 2, xp.pi]), xp=xp)
    Array([  0.,  90., 180.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(x)
    if xp.isdtype(x.dtype, "integral"):
        dtype = _agnostic._inspection.default_dtype(xp, device=_compat.device(x))
        x = xp.astype(x, dtype)
    elif not xp.isdtype(x.dtype, ("real floating", "complex floating")):
        msg = "`x` must have an integral, real floating, or complex floating dtype."
        raise TypeError(msg)

    if _compat.is_jax_namespace(xp) or (
        not xp.isdtype(x.dtype, "complex floating")
        and (
            _compat.is_numpy_namespace(xp)
            or _compat.is_cupy_namespace(xp)
            or _compat.is_torch_namespace(xp)
            or _compat.is_dask_namespace(xp)
        )
    ):
        return xp.rad2deg(x)

    return _agnostic._elementwise.rad2deg(x, xp=xp)


def isclose(
    a: Array | complex,
    b: Array | complex,
    *,
    rtol: float = 1e-05,
    atol: float = 1e-08,
    equal_nan: bool = False,
    xp: ArrayNamespace | None = None,
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
    xp = _compat.array_namespace(a, b) if xp is None else xp

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    if _compat.is_torch_namespace(xp):
        a, b = _helpers.asarrays(a, b, xp=xp)  # Array API 2024.12 support
        return xp.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    return _agnostic._elementwise.isclose(
        a, b, rtol=rtol, atol=atol, equal_nan=equal_nan, xp=xp
    )


def nan_to_num(
    x: Array | complex,
    /,
    *,
    fill_value: float = 0.0,
    xp: ArrayNamespace | None = None,
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
    >>> xpx.nan_to_num(xp.inf, xp=xp)
    Array(1.79769313e+308, dtype=array_api_strict.float64)
    >>> xpx.nan_to_num(-xp.inf, xp=xp)
    Array(-1.79769313e+308, dtype=array_api_strict.float64)
    >>> xpx.nan_to_num(xp.nan, xp=xp)
    Array(0., dtype=array_api_strict.float64)
    >>> x = xp.asarray([xp.inf, -xp.inf, xp.nan, -128, 128])
    >>> xpx.nan_to_num(x)
    Array([ 1.79769313e+308, -1.79769313e+308,  0.00000000e+000,
           -1.28000000e+002,  1.28000000e+002],
          dtype=array_api_strict.float64)
    >>> y = xp.asarray([complex(xp.inf, xp.nan), xp.nan, complex(xp.nan, xp.inf)])
    >>> xpx.nan_to_num(y)
    Array([1.79769313e+308+0.00000000e+000j,
           0.00000000e+000+0.00000000e+000j,
           0.00000000e+000+1.79769313e+308j],
          dtype=array_api_strict.complex128)
    """
    if isinstance(fill_value, complex):
        msg = "Complex fill values are not supported."
        raise TypeError(msg)

    xp = _compat.array_namespace(x) if xp is None else xp

    # for scalars we want to output an array
    y = xp.asarray(x)

    if (
        _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_numpy_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return xp.nan_to_num(y, nan=fill_value)

    return _agnostic._elementwise.nan_to_num(y, fill_value=fill_value, xp=xp)


def sinc(x: Array, /, *, xp: ArrayNamespace | None = None) -> Array:
    r"""
    Return the normalized sinc function.

    The sinc function is equal to :math:`\sin(\pi x)/(\pi x)` for any argument
    :math:`x\ne 0`. ``sinc(0)`` takes the limit value 1, making ``sinc`` not
    only everywhere continuous but also infinitely differentiable.

    .. note::

        Note the normalization factor of ``pi`` used in the definition.
        This is the most commonly used definition in signal processing.
        Use ``sinc(x / xp.pi)`` to obtain the unnormalized sinc function
        :math:`\sin(x)/x` that is more common in mathematics.

    Parameters
    ----------
    x : array
        Array (possibly multi-dimensional) of values for which to calculate
        ``sinc(x)``. Must have a real floating point dtype.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        ``sinc(x)`` calculated elementwise, which has the same shape as the input.

    Notes
    -----
    The name sinc is short for "sine cardinal" or "sinus cardinalis".

    The sinc function is used in various signal processing applications,
    including in anti-aliasing, in the construction of a Lanczos resampling
    filter, and in interpolation.

    For bandlimited interpolation of discrete-time signals, the ideal
    interpolation kernel is proportional to the sinc function.

    References
    ----------
    #. Weisstein, Eric W. "Sinc Function." From MathWorld--A Wolfram Web
       Resource. https://mathworld.wolfram.com/SincFunction.html
    #. Wikipedia, "Sinc function",
       https://en.wikipedia.org/wiki/Sinc_function

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> x = xp.linspace(-4, 4, 41)
    >>> xpx.sinc(x, xp=xp)
    Array([-3.89817183e-17, -4.92362781e-02,
           -8.40918587e-02, -8.90384387e-02,
           -5.84680802e-02,  3.89817183e-17,
            6.68206631e-02,  1.16434881e-01,
            1.26137788e-01,  8.50444803e-02,
           -3.89817183e-17, -1.03943254e-01,
           -1.89206682e-01, -2.16236208e-01,
           -1.55914881e-01,  3.89817183e-17,
            2.33872321e-01,  5.04551152e-01,
            7.56826729e-01,  9.35489284e-01,
            1.00000000e+00,  9.35489284e-01,
            7.56826729e-01,  5.04551152e-01,
            2.33872321e-01,  3.89817183e-17,
           -1.55914881e-01, -2.16236208e-01,
           -1.89206682e-01, -1.03943254e-01,
           -3.89817183e-17,  8.50444803e-02,
            1.26137788e-01,  1.16434881e-01,
            6.68206631e-02,  3.89817183e-17,
           -5.84680802e-02, -8.90384387e-02,
           -8.40918587e-02, -4.92362781e-02,
           -3.89817183e-17], dtype=array_api_strict.float64)
    """

    if xp is None:
        xp = _compat.array_namespace(x)

    if not xp.isdtype(x.dtype, "real floating"):
        err_msg = "`x` must have a real floating data type."
        raise ValueError(err_msg)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_torch_namespace(xp)
        or _compat.is_dask_namespace(xp)
    ):
        return xp.sinc(x)

    return _agnostic._elementwise.sinc(x, xp=xp)
