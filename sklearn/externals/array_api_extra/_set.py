"""Delegation layer for set functions."""

from . import _agnostic
from ._lib import _compat, _helpers
from ._lib._typing import Array, ArrayNamespace

__all__ = ["isin", "nunique", "setdiff1d", "union1d"]


def isin(
    a: Array,
    b: Array,
    /,
    *,
    assume_unique: bool = False,
    invert: bool = False,
    kind: str | None = None,
    xp: ArrayNamespace | None = None,
) -> Array:
    """
    Determine whether each element in `a` is present in `b`.

    This is :func:`array_api.isin`, with additional `assume_unique`
    and `kind` parameters.

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
        xp = _compat.array_namespace(a, b)

    if _compat.is_numpy_namespace(xp):
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert, kind=kind)
    if _compat.is_jax_namespace(xp):
        if kind is None:
            kind = "auto"
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert, method=kind)
    if (
        _compat.is_cupy_namespace(xp)
        or _compat.is_torch_namespace(xp)
        or _compat.is_dask_namespace(xp)
    ):
        return xp.isin(a, b, assume_unique=assume_unique, invert=invert)

    return _agnostic._set.isin(a, b, assume_unique=assume_unique, invert=invert, xp=xp)


def nunique(x: Array, /, *, xp: ArrayNamespace | None = None) -> Array:
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
        xp = _compat.array_namespace(x)

    if _compat.is_jax_array(x):
        # size= is JAX-specific
        # https://github.com/data-apis/array-api/issues/883
        _, counts = xp.unique_counts(x, size=_compat.size(x))
        return (counts > 0).sum()

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or (
            _compat.is_torch_namespace(xp)
            and _helpers.capabilities(xp, x)["data-dependent shapes"]
        )
    ):
        _, counts = xp.unique_counts(x)
        return xp.asarray(_compat.size(counts), device=_compat.device(x))

    return _agnostic._set.nunique(x, xp=xp)


def setdiff1d(
    x1: Array | complex,
    x2: Array | complex,
    /,
    *,
    assume_unique: bool = False,
    xp: ArrayNamespace | None = None,
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
        xp = _compat.array_namespace(x1, x2)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        x1, x2 = _helpers.asarrays(x1, x2, xp=xp)
        return xp.setdiff1d(x1, x2, assume_unique=assume_unique)

    return _agnostic._set.setdiff1d(x1, x2, assume_unique=assume_unique, xp=xp)


def union1d(a: Array, b: Array, /, *, xp: ArrayNamespace | None = None) -> Array:
    """
    Find the union of two arrays.

    Return the unique, sorted array of values that are in either of the two
    input arrays.

    Parameters
    ----------
    a, b : Array
        Input arrays. They are flattened internally if they are not already 1D.

    xp : array_namespace, optional
        The standard-compatible namespace for `a` and `b`. Default: infer.

    Returns
    -------
    Array
        Unique, sorted union of the input arrays.

    See Also
    --------
    jax.numpy.union1d : Corresponding function in JAX.

    Notes
    -----
    This function is not compatible with `jax.jit`.
    See the docstring of the corresponding JAX function for more information.
    """
    if xp is None:
        xp = _compat.array_namespace(a, b)

    if (
        _compat.is_numpy_namespace(xp)
        or _compat.is_cupy_namespace(xp)
        or _compat.is_dask_namespace(xp)
        or _compat.is_jax_namespace(xp)
    ):
        return xp.union1d(a, b)

    return _agnostic._set.union1d(a, b, xp=xp)
