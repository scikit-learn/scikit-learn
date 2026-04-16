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
from ._lib._utils._helpers import asarrays, eager_shape
from ._lib._utils._typing import Array, DType

__all__ = [
    "atleast_nd",
    "cov",
    "create_diagonal",
    "expand_dims",
    "isclose",
    "nan_to_num",
    "one_hot",
    "pad",
    "sinc",
]


def atleast_nd(x: Array, /, *, ndim: int, xp: ModuleType | None = None) -> Array:
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
        xp = array_namespace(x)

    if 1 <= ndim <= 2 and (
        is_numpy_namespace(xp)
        or is_jax_namespace(xp)
        or is_dask_namespace(xp)
        or is_cupy_namespace(xp)
        or is_torch_namespace(xp)
    ):
        return getattr(xp, f"atleast_{ndim}d")(x)

    return _funcs.atleast_nd(x, ndim=ndim, xp=xp)


def cov(m: Array, /, *, xp: ModuleType | None = None) -> Array:
    """
    Estimate a covariance matrix (or a stack of covariance matrices).

    Covariance indicates the level to which two variables vary together.
    If we examine *N*-dimensional samples, :math:`X = [x_1, x_2, ... x_N]^T`,
    each with *M* observations, then element :math:`C_{ij}` of the
    :math:`N \times N` covariance matrix is the covariance of
    :math:`x_i` and :math:`x_j`. The element :math:`C_{ii}` is the variance
    of :math:`x_i`.

    With the exception of supporting batch input, this provides a subset of
    the functionality of ``numpy.cov``.

    Parameters
    ----------
    m : array
        An array of shape ``(..., N, M)`` whose innermost two dimensions
        contain *M* observations of *N* variables. That is,
        each row of `m` represents a variable, and each column a single
        observation of all those variables.
    xp : array_namespace, optional
        The standard-compatible namespace for `m`. Default: infer.

    Returns
    -------
    array
        An array having shape (..., N, N) whose innermost two dimensions represent
        the covariance matrix of the variables.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx

    Consider two variables, :math:`x_0` and :math:`x_1`, which
    correlate perfectly, but in opposite directions:

    >>> x = xp.asarray([[0, 2], [1, 1], [2, 0]]).T
    >>> x
    Array([[0, 1, 2],
           [2, 1, 0]], dtype=array_api_strict.int64)

    Note how :math:`x_0` increases while :math:`x_1` decreases. The covariance
    matrix shows this clearly:

    >>> xpx.cov(x, xp=xp)
    Array([[ 1., -1.],
           [-1.,  1.]], dtype=array_api_strict.float64)

    Note that element :math:`C_{0,1}`, which shows the correlation between
    :math:`x_0` and :math:`x_1`, is negative.

    Further, note how `x` and `y` are combined:

    >>> x = xp.asarray([-2.1, -1,  4.3])
    >>> y = xp.asarray([3,  1.1,  0.12])
    >>> X = xp.stack((x, y), axis=0)
    >>> xpx.cov(X, xp=xp)
    Array([[11.71      , -4.286     ],
           [-4.286     ,  2.14413333]], dtype=array_api_strict.float64)

    >>> xpx.cov(x, xp=xp)
    Array(11.71, dtype=array_api_strict.float64)

    >>> xpx.cov(y, xp=xp)
    Array(2.14413333, dtype=array_api_strict.float64)

    Input with more than two dimensions is treated as a stack of
    two-dimensional input.

    >>> stack = xp.stack((X, 2*X))
    >>> xpx.cov(stack)
    Array([[[ 11.71      ,  -4.286     ],
            [ -4.286     ,   2.14413333]],

           [[ 46.84      , -17.144     ],
            [-17.144     ,   8.57653333]]], dtype=array_api_strict.float64)
    """

    if xp is None:
        xp = array_namespace(m)

    if (
        is_numpy_namespace(xp)
        or is_cupy_namespace(xp)
        or is_torch_namespace(xp)
        or is_dask_namespace(xp)
        or is_jax_namespace(xp)
    ) and m.ndim <= 2:
        return xp.cov(m)

    return _funcs.cov(m, xp=xp)


def create_diagonal(
    x: Array, /, *, offset: int = 0, xp: ModuleType | None = None
) -> Array:
    """
    Construct a diagonal array.

    Parameters
    ----------
    x : array
        An array having shape ``(*batch_dims, k)``.
    offset : int, optional
        Offset from the leading diagonal (default is ``0``).
        Use positive ints for diagonals above the leading diagonal,
        and negative ints for diagonals below the leading diagonal.
    xp : array_namespace, optional
        The standard-compatible namespace for `x`. Default: infer.

    Returns
    -------
    array
        An array having shape ``(*batch_dims, k+abs(offset), k+abs(offset))`` with `x`
        on the diagonal (offset by `offset`).

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx
    >>> x = xp.asarray([2, 4, 8])

    >>> xpx.create_diagonal(x, xp=xp)
    Array([[2, 0, 0],
           [0, 4, 0],
           [0, 0, 8]], dtype=array_api_strict.int64)

    >>> xpx.create_diagonal(x, offset=-2, xp=xp)
    Array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [2, 0, 0, 0, 0],
           [0, 4, 0, 0, 0],
           [0, 0, 8, 0, 0]], dtype=array_api_strict.int64)
    """
    if xp is None:
        xp = array_namespace(x)

    if x.ndim == 0:
        err_msg = "`x` must be at least 1-dimensional."
        raise ValueError(err_msg)

    if is_torch_namespace(xp):
        return xp.diag_embed(x, offset=offset, dim1=-2, dim2=-1)

    if (is_dask_namespace(xp) or is_cupy_namespace(xp)) and x.ndim < 2:
        return xp.diag(x, k=offset)

    if (is_jax_namespace(xp) or is_numpy_namespace(xp)) and x.ndim < 3:
        batch_dim, n = eager_shape(x)[:-1], eager_shape(x, -1)[0] + abs(offset)
        return xp.reshape(xp.diag(x, k=offset), (*batch_dim, n, n))

    return _funcs.create_diagonal(x, offset=offset, xp=xp)


def expand_dims(
    a: Array, /, *, axis: int | tuple[int, ...] = (0,), xp: ModuleType | None = None
) -> Array:
    """
    Expand the shape of an array.

    Insert (a) new axis/axes that will appear at the position(s) specified by
    `axis` in the expanded array shape.

    This is ``xp.expand_dims`` for `axis` an int *or a tuple of ints*.
    Roughly equivalent to ``numpy.expand_dims`` for NumPy arrays.

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
        xp = array_namespace(a)

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

    if is_numpy_namespace(xp) or is_dask_namespace(xp) or is_jax_namespace(xp):
        return xp.expand_dims(a, axis=axis)

    return _funcs.expand_dims(a, axis=axis, xp=xp)


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


def setdiff1d(
    x1: Array | complex,
    x2: Array | complex,
    /,
    *,
    assume_unique: bool = False,
    xp: ModuleType | None = None,
) -> Array:
    """
    Find the set difference of two arrays.

    Return the unique values in `x1` that are not in `x2`.

    Parameters
    ----------
    x1 : array | int | float | complex | bool
        Input array.
    x2 : array
        Input comparison array.
    assume_unique : bool
        If ``True``, the input arrays are both assumed to be unique, which
        can speed up the calculation. Default is ``False``.
    xp : array_namespace, optional
        The standard-compatible namespace for `x1` and `x2`. Default: infer.

    Returns
    -------
    array
        1D array of values in `x1` that are not in `x2`. The result
        is sorted when `assume_unique` is ``False``, but otherwise only sorted
        if the input is sorted.

    Examples
    --------
    >>> import array_api_strict as xp
    >>> import array_api_extra as xpx

    >>> x1 = xp.asarray([1, 2, 3, 2, 4, 1])
    >>> x2 = xp.asarray([3, 4, 5, 6])
    >>> xpx.setdiff1d(x1, x2, xp=xp)
    Array([1, 2], dtype=array_api_strict.int64)
    """

    if xp is None:
        xp = array_namespace(x1, x2)

    if is_numpy_namespace(xp) or is_cupy_namespace(xp) or is_jax_namespace(xp):
        x1, x2 = asarrays(x1, x2, xp=xp)
        return xp.setdiff1d(x1, x2, assume_unique=assume_unique)

    return _funcs.setdiff1d(x1, x2, assume_unique=assume_unique, xp=xp)


def sinc(x: Array, /, *, xp: ModuleType | None = None) -> Array:
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
        xp = array_namespace(x)

    if not xp.isdtype(x.dtype, "real floating"):
        err_msg = "`x` must have a real floating data type."
        raise ValueError(err_msg)

    if (
        is_numpy_namespace(xp)
        or is_cupy_namespace(xp)
        or is_jax_namespace(xp)
        or is_torch_namespace(xp)
        or is_dask_namespace(xp)
    ):
        return xp.sinc(x)

    return _funcs.sinc(x, xp=xp)


def partition(
    a: Array,
    kth: int,
    /,
    axis: int | None = -1,
    *,
    xp: ModuleType | None = None,
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
        xp = array_namespace(a)
    if a.ndim < 1:
        msg = "`a` must be at least 1-dimensional"
        raise TypeError(msg)
    if axis is None:
        return partition(xp.reshape(a, (-1,)), kth, axis=0, xp=xp)
    (size,) = eager_shape(a, axis)
    if not (0 <= kth < size):
        msg = f"kth(={kth}) out of bounds [0 {size})"
        raise ValueError(msg)

    # Delegate where possible.
    if is_numpy_namespace(xp) or is_cupy_namespace(xp) or is_jax_namespace(xp):
        return xp.partition(a, kth, axis=axis)

    # Use top-k when possible:
    if is_torch_namespace(xp):
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

    return _funcs.partition(a, kth, axis=axis, xp=xp)


def argpartition(
    a: Array,
    kth: int,
    /,
    axis: int | None = -1,
    *,
    xp: ModuleType | None = None,
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
        xp = array_namespace(a)
    if is_pydata_sparse_namespace(xp):
        msg = "Not implemented for sparse backend: no argsort"
        raise NotImplementedError(msg)
    if a.ndim < 1:
        msg = "`a` must be at least 1-dimensional"
        raise TypeError(msg)
    if axis is None:
        return argpartition(xp.reshape(a, (-1,)), kth, axis=0, xp=xp)
    (size,) = eager_shape(a, axis)
    if not (0 <= kth < size):
        msg = f"kth(={kth}) out of bounds [0 {size})"
        raise ValueError(msg)

    # Delegate where possible.
    if is_numpy_namespace(xp) or is_cupy_namespace(xp) or is_jax_namespace(xp):
        return xp.argpartition(a, kth, axis=axis)

    # Use top-k when possible:
    if is_torch_namespace(xp):
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

    return _funcs.argpartition(a, kth, axis=axis, xp=xp)


def isin(
    a: Array,
    b: Array,
    /,
    *,
    assume_unique: bool = False,
    invert: bool = False,
    kind: str | None = None,
    xp: ModuleType | None = None,
) -> Array:
    """
    Determine whether each element in `a` is present in `b`.

    Return a boolean array of the same shape as `a` that is True for elements
    that are in `b` and False otherwise.

    Parameters
    ----------
    a : array
        Input elements.
    b : array
        The elements against which to test each element of `a`.
    assume_unique : bool, optional
        If True, the input arrays are both assumed to be unique which can speed
        up the calculation. Default: False.
    invert : bool, optional
        If True, the values in the returned array are inverted. Default: False.
    kind : str | None, optional
        The algorithm or method to use. This will not affect the final result,
        but will affect the speed and memory use.
        For NumPy the options are {None, "sort", "table"}.
        For Jax the mapped parameter is instead `method` and the options are
        {"compare_all", "binary_search", "sort", and "auto" (default)}
        For CuPy, Dask, Torch and the default case this parameter is not present and
        thus ignored. Default: None.
    xp : array_namespace, optional
        The standard-compatible namespace for `a` and `b`. Default: infer.

    Returns
    -------
    array
        An array having the same shape as that of `a` that is True for elements
        that are in `b` and False otherwise.
    """
    if xp is None:
        xp = array_namespace(a, b)

    if is_numpy_namespace(xp):
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert, kind=kind)
    if is_jax_namespace(xp):
        if kind is None:
            kind = "auto"
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert, method=kind)
    if is_cupy_namespace(xp) or is_torch_namespace(xp) or is_dask_namespace(xp):
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert)

    return _funcs.isin(a, b, assume_unique=assume_unique, invert=invert, xp=xp)
