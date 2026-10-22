"""Delegation layer for statistical functions."""

from . import _agnostic
from ._lib import _compat
from ._lib._typing import Array, ArrayNamespace

__all__ = ["cov", "nanmax", "nanmean", "nanmin", "nansum"]


def cov(
    m: Array,
    /,
    *,
    axis: int = -1,
    correction: float = 1,
    fweights: Array | None = None,
    aweights: Array | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Estimate a covariance matrix (or a stack of covariance matrices).

    Covariance indicates the level to which two variables vary together.
    If we examine *N*-dimensional samples, :math:`X = [x_1, x_2, ... x_N]^T`,
    each with *M* observations, then element :math:`C_{ij}` of the
    :math:`N \\times N` covariance matrix is the covariance of
    :math:`x_i` and :math:`x_j`. The element :math:`C_{ii}` is the variance
    of :math:`x_i`.

    Extends :func:`numpy.cov` with support for batch input.
    Naming follows the array API conventions used elsewhere in
    this library (``axis``, ``correction``) rather than the NumPy spellings
    (``rowvar``, ``bias``, ``ddof``); see Notes for the mapping.

    Parameters
    ----------
    m : array
        An array of shape ``(..., N, M)`` whose innermost two dimensions
        contain *M* observations of *N* variables by default. The axis of
        observations is controlled by `axis`.
    axis : int, optional
        Axis of `m` containing the observations. Default: ``-1`` (the last
        axis), matching the array API convention. Use ``axis=-2`` (or ``0``
        for 2-D input) to treat each column as a variable, which
        corresponds to ``rowvar=False`` in :func:`numpy.cov`.
    correction : int or float, optional
        Degrees of freedom correction: normalization divides by
        ``N - correction`` (for unweighted input). Default: ``1``, which
        gives the unbiased estimate (matches :func:`numpy.cov` default of
        ``bias=False``). Set to ``0`` for the biased estimate (``N``
        normalization). Corresponds to ``ddof`` in :func:`numpy.cov` and to
        ``correction`` in :func:`numpy.var`/:func:`numpy.std` and
        :func:`torch.cov`.
        Non-integer values are allowed for advanced use cases: the
        unbiased correction for weighted observations depends on the
        sum and dispersion of the weights and is generally not an
        integer, and autocorrelated data may also require a fractional
        correction. Non-integer ``correction`` routes through the
        generic implementation because :func:`numpy.cov`'s ``ddof`` and
        :func:`torch.cov`'s ``correction`` both require integers.
    fweights : array, optional
        1-D array of integer frequency weights: the number of times each
        observation is repeated. Same as ``fweights`` in
        :func:`numpy.cov`/:func:`torch.cov`.
    aweights : array, optional
        1-D array of observation-vector weights (analytic weights). Larger
        values mark more important observations. Same as ``aweights`` in
        :func:`numpy.cov`/:func:`torch.cov`.
    xp : array_namespace, optional
        The standard-compatible namespace for `m`. Default: infer.

    Returns
    -------
    array
        An array having shape ``(..., N, N)`` whose innermost two dimensions represent
        the covariance matrix of the variables.

    Notes
    -----
    Mapping from :func:`numpy.cov` to this function::

        numpy.cov(m, rowvar=True)           -> cov(m, axis=-1)   # default
        numpy.cov(m, rowvar=False)          -> cov(m, axis=-2)
        numpy.cov(m, bias=True)             -> cov(m, correction=0)
        numpy.cov(m, ddof=k)                -> cov(m, correction=k)
        numpy.cov(m, fweights=f)            -> cov(m, fweights=f)
        numpy.cov(m, aweights=a)            -> cov(m, aweights=a)

    A ``RuntimeWarning`` is emitted for non-positive effective degrees of
    freedom when the effective normalizer can be checked without materializing
    a lazy array. When the normalizer itself is lazy (e.g. for weighted Dask
    inputs), this check is skipped; choose ``correction`` and weights such that
    it is positive.

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

    The normalization can be adjusted with `correction`, and observations
    can be weighted with integer frequencies `fweights` or importance
    weights `aweights`:

    >>> x = xp.asarray([0., 1., 2., 3., 4.])
    >>> xpx.cov(x, xp=xp)  # unbiased variance: divide by N - 1
    Array(2.5, dtype=array_api_strict.float64)
    >>> xpx.cov(x, correction=0, xp=xp)  # biased variance: divide by N
    Array(2., dtype=array_api_strict.float64)

    Giving the two extreme observations frequency 2 via `fweights` is
    equivalent to repeating them in `x`:

    >>> xpx.cov(x, fweights=xp.asarray([2, 1, 1, 1, 2]), xp=xp)
    Array(3., dtype=array_api_strict.float64)
    >>> xpx.cov(xp.asarray([0., 0., 1., 2., 3., 4., 4.]), xp=xp)
    Array(3., dtype=array_api_strict.float64)

    `aweights` instead adjusts the relative importance of observations,
    here down-weighting the two extremes:

    >>> xpx.cov(x, aweights=xp.asarray([0.5, 1., 1., 1., 0.5]), xp=xp)
    Array(1.92, dtype=array_api_strict.float64)
    """

    if xp is None:
        xp = _compat.array_namespace(m, fweights, aweights)

    # Validate axis against m.ndim.
    ndim = max(m.ndim, 1)
    if not -ndim <= axis < ndim:
        msg = f"axis {axis} is out of bounds for array of dimension {m.ndim}"
        raise IndexError(msg)

    # Normalize: observations on the last axis. After this, every backend
    # sees the same convention and we never need to deal with `rowvar`.
    if m.ndim >= 2 and axis not in (-1, m.ndim - 1):
        m = xp.moveaxis(m, axis, -1)

    # `numpy.cov` (and cupy/dask/jax) require integer `ddof`; `torch.cov`
    # requires integer `correction`. For non-integer-valued `correction`,
    # fall through to the generic implementation.
    integer_correction = float(correction).is_integer()
    has_weights = fweights is not None or aweights is not None

    if m.ndim <= 2 and integer_correction:
        # Not just for static typing: `correction` may be an integer-valued
        # float such as 1.0, which `torch.cov` rejects at runtime.
        int_correction = int(correction)
        if _compat.is_torch_namespace(xp):
            fw = None if fweights is None else xp.asarray(fweights)
            aw = None if aweights is None else xp.asarray(aweights)
            return xp.cov(m, correction=int_correction, fweights=fw, aweights=aw)
        # `dask.array.cov` forces `.compute()` whenever weights are given:
        # its internal `if fact <= 0` check on a lazy 0-D scalar triggers
        # materialization. Route to the generic impl, which is fully lazy
        # because it only does sum/matmul and skips that scalar check.
        if (
            _compat.is_numpy_namespace(xp)
            or _compat.is_cupy_namespace(xp)
            or _compat.is_jax_namespace(xp)
            or (_compat.is_dask_namespace(xp) and not has_weights)
        ):
            return xp.cov(
                m,
                ddof=int_correction,
                fweights=fweights,
                aweights=aweights,
            )

    return _agnostic._statistical.cov(
        m,
        correction=correction,
        fweights=fweights,
        aweights=aweights,
        xp=xp,
    )


def nanmin(
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Return the minimum of the array elements along a given axis, ignoring NaNs.

    Parameters
    ----------
    a : Array
        Input array.
    axis : int or tuple of ints or None, optional
        Axis or axes along which the minimum is computed. The default is to compute
        the minimum of the flattened array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a`. Default: infer.

    Returns
    -------
    array
        An array of minimum values along the given axis, ignoring NaNs.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> a = xp.asarray([[5, 3, xp.nan, 1], [4, xp.nan, 2, xp.nan]])
    >>> xpx.nanmin(a)
    Array(1., dtype=array_api_strict.float64)
    >>> xpx.nanmin(a, axis=0)
    Array([4., 3., 2., 1.], dtype=array_api_strict.float64)
    >>> xpx.nanmin(a, axis=1)
    Array([1., 2.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(a)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.nanmin(a, axis=axis)

    return _agnostic._statistical.nanmin(a, axis=axis, xp=xp)


def nanmax(
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Return the maximum of the array elements along a given axis, ignoring NaNs.

    Parameters
    ----------
    a : Array
        Input array.
    axis : int or tuple of ints or None, optional
        Axis or axes along which the maximum is computed. The default is to compute
        the maximum of the flattened array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a`. Default: infer.

    Returns
    -------
    array
        An array of maximum values along the given axis, ignoring NaNs.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> a = xp.asarray([[5, 3, xp.nan, 6], [4, xp.nan, 2, xp.nan]])
    >>> xpx.nanmax(a)
    Array(6., dtype=array_api_strict.float64)
    >>> xpx.nanmax(a, axis=0)
    Array([5., 3., 2., 6.], dtype=array_api_strict.float64)
    >>> xpx.nanmax(a, axis=1)
    Array([6., 4.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(a)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.nanmax(a, axis=axis)

    return _agnostic._statistical.nanmax(a, axis=axis, xp=xp)


def nanmean(
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Return the mean of the array elements along a given axis, ignoring NaNs.

    Parameters
    ----------
    a : Array
        Input array.
    axis : int or tuple of ints or None, optional
        Axis or axes along which the mean is computed. The default is to compute
        the mean of the flattened array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a`. Default: infer.

    Returns
    -------
    array
        An array of mean values along the given axis, ignoring NaNs.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> a = xp.asarray([[5, 3, xp.nan, 1], [4, xp.nan, 2, xp.nan]])
    >>> xpx.nanmean(a)
    Array(3., dtype=array_api_strict.float64)
    >>> xpx.nanmean(a, axis=0)
    Array([4.5, 3. , 2. , 1. ], dtype=array_api_strict.float64)
    >>> xpx.nanmean(a, axis=1)
    Array([3., 3.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(a)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return xp.nanmean(a, axis=axis)

    return _agnostic._statistical.nanmean(a, axis=axis, xp=xp)


def nansum(
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Return the sum of the array elements along a given axis, ignoring NaNs.

    Parameters
    ----------
    a : Array
        Input array.
    axis : int or tuple of ints or None, optional
        Axis or axes along which the sum is computed. The default is to compute
        the sum of the flattened array.
    xp : array_namespace, optional
        The standard-compatible namespace for `a`. Default: infer.

    Returns
    -------
    array
        An array of sum values along the given axis, ignoring NaNs.

    Examples
    --------
    >>> import array_api_extra as xpx
    >>> import array_api_strict as xp
    >>> a = xp.asarray([[5, 3, xp.nan, 1], [4, xp.nan, 2, xp.nan]])
    >>> xpx.nansum(a)
    Array(15., dtype=array_api_strict.float64)
    >>> xpx.nansum(a, axis=0)
    Array([9., 3., 2., 1.], dtype=array_api_strict.float64)
    >>> xpx.nansum(a, axis=1)
    Array([9., 6.], dtype=array_api_strict.float64)
    """
    if xp is None:
        xp = _compat.array_namespace(a)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
        or _compat.is_torch_namespace(xp)
    ):
        return xp.nansum(a, axis=axis)

    return _agnostic._statistical.nansum(a, axis=axis, xp=xp)
