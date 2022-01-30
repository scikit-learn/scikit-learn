from __future__ import annotations

import numpy as np

from pandas._libs import lib
from pandas._typing import (
    ArrayLike,
    Scalar,
    npt,
)
from pandas.compat.numpy import np_percentile_argname

from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
)


def quantile_compat(
    values: ArrayLike, qs: npt.NDArray[np.float64], interpolation: str
) -> ArrayLike:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray or ExtensionArray
    qs : np.ndarray[float64]
    interpolation : str

    Returns
    -------
    np.ndarray or ExtensionArray
    """
    if isinstance(values, np.ndarray):
        fill_value = na_value_for_dtype(values.dtype, compat=False)
        mask = isna(values)
        return quantile_with_mask(values, mask, fill_value, qs, interpolation)
    else:
        return values._quantile(qs, interpolation)


def quantile_with_mask(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    fill_value,
    qs: npt.NDArray[np.float64],
    interpolation: str,
) -> np.ndarray:
    """
    Compute the quantiles of the given values for each quantile in `qs`.

    Parameters
    ----------
    values : np.ndarray
        For ExtensionArray, this is _values_for_factorize()[0]
    mask : np.ndarray[bool]
        mask = isna(values)
        For ExtensionArray, this is computed before calling _value_for_factorize
    fill_value : Scalar
        The value to interpret fill NA entries with
        For ExtensionArray, this is _values_for_factorize()[1]
    qs : np.ndarray[float64]
    interpolation : str
        Type of interpolation

    Returns
    -------
    np.ndarray

    Notes
    -----
    Assumes values is already 2D.  For ExtensionArray this means np.atleast_2d
    has been called on _values_for_factorize()[0]

    Quantile is computed along axis=1.
    """
    assert values.ndim == 2

    is_empty = values.shape[1] == 0

    if is_empty:
        # create the array of na_values
        # 2d len(values) * len(qs)
        flat = np.array([fill_value] * len(qs))
        result = np.repeat(flat, len(values)).reshape(len(values), len(qs))
    else:
        result = _nanpercentile(
            values,
            qs * 100.0,
            na_value=fill_value,
            mask=mask,
            interpolation=interpolation,
        )

        result = np.array(result, copy=False)
        result = result.T

    return result


def _nanpercentile_1d(
    values: np.ndarray,
    mask: npt.NDArray[np.bool_],
    qs: npt.NDArray[np.float64],
    na_value: Scalar,
    interpolation,
) -> Scalar | np.ndarray:
    """
    Wrapper for np.percentile that skips missing values, specialized to
    1-dimensional case.

    Parameters
    ----------
    values : array over which to find quantiles
    mask : ndarray[bool]
        locations in values that should be considered missing
    qs : np.ndarray[float64] of quantile indices to find
    na_value : scalar
        value to return for empty or all-null values
    interpolation : str

    Returns
    -------
    quantiles : scalar or array
    """
    # mask is Union[ExtensionArray, ndarray]
    values = values[~mask]

    if len(values) == 0:
        return np.array([na_value] * len(qs), dtype=values.dtype)

    return np.percentile(values, qs, **{np_percentile_argname: interpolation})


def _nanpercentile(
    values: np.ndarray,
    qs: npt.NDArray[np.float64],
    *,
    na_value,
    mask: npt.NDArray[np.bool_],
    interpolation,
):
    """
    Wrapper for np.percentile that skips missing values.

    Parameters
    ----------
    values : np.ndarray[ndim=2]  over which to find quantiles
    qs : np.ndarray[float64] of quantile indices to find
    na_value : scalar
        value to return for empty or all-null values
    mask : np.ndarray[bool]
        locations in values that should be considered missing
    interpolation : str

    Returns
    -------
    quantiles : scalar or array
    """

    if values.dtype.kind in ["m", "M"]:
        # need to cast to integer to avoid rounding errors in numpy
        result = _nanpercentile(
            values.view("i8"),
            qs=qs,
            na_value=na_value.view("i8"),
            mask=mask,
            interpolation=interpolation,
        )

        # Note: we have to do `astype` and not view because in general we
        #  have float result at this point, not i8
        return result.astype(values.dtype)

    if not lib.is_scalar(mask) and mask.any():
        # Caller is responsible for ensuring mask shape match
        assert mask.shape == values.shape
        result = [
            _nanpercentile_1d(val, m, qs, na_value, interpolation=interpolation)
            for (val, m) in zip(list(values), list(mask))
        ]
        result = np.array(result, dtype=values.dtype, copy=False).T
        return result
    else:
        return np.percentile(
            values, qs, axis=1, **{np_percentile_argname: interpolation}
        )
