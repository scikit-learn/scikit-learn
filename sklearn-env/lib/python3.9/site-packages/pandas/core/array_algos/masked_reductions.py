"""
masked_reductions.py is for reduction algorithms using a mask-based approach
for missing values.
"""

from typing import (
    Callable,
    Optional,
)

import numpy as np

from pandas._libs import missing as libmissing

from pandas.core.nanops import check_below_min_count


def _sumprod(
    func: Callable,
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: Optional[int] = None,
):
    """
    Sum or product for 1D masked array.

    Parameters
    ----------
    func : np.sum or np.prod
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation).
    mask : np.ndarray
        Boolean numpy array (True values indicate missing values).
    skipna : bool, default True
        Whether to skip NA.
    min_count : int, default 0
        The required number of valid values to perform the operation. If fewer than
        ``min_count`` non-NA values are present the result will be NA.
    axis : int, optional, default None
    """
    if not skipna:
        if mask.any(axis=axis) or check_below_min_count(values.shape, None, min_count):
            return libmissing.NA
        else:
            return func(values, axis=axis)
    else:
        if check_below_min_count(values.shape, mask, min_count) and (
            axis is None or values.ndim == 1
        ):
            return libmissing.NA

        return func(values, where=~mask, axis=axis)


def sum(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: Optional[int] = None,
):
    return _sumprod(
        np.sum, values=values, mask=mask, skipna=skipna, min_count=min_count, axis=axis
    )


def prod(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    min_count: int = 0,
    axis: Optional[int] = None,
):
    return _sumprod(
        np.prod, values=values, mask=mask, skipna=skipna, min_count=min_count, axis=axis
    )


def _minmax(
    func: Callable,
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    axis: Optional[int] = None,
):
    """
    Reduction for 1D masked array.

    Parameters
    ----------
    func : np.min or np.max
    values : np.ndarray
        Numpy array with the values (can be of any dtype that support the
        operation).
    mask : np.ndarray
        Boolean numpy array (True values indicate missing values).
    skipna : bool, default True
        Whether to skip NA.
    axis : int, optional, default None
    """
    if not skipna:
        if mask.any() or not values.size:
            # min/max with empty array raise in numpy, pandas returns NA
            return libmissing.NA
        else:
            return func(values)
    else:
        subset = values[~mask]
        if subset.size:
            return func(subset)
        else:
            # min/max with empty array raise in numpy, pandas returns NA
            return libmissing.NA


def min(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    axis: Optional[int] = None,
):
    return _minmax(np.min, values=values, mask=mask, skipna=skipna, axis=axis)


def max(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    skipna: bool = True,
    axis: Optional[int] = None,
):
    return _minmax(np.max, values=values, mask=mask, skipna=skipna, axis=axis)


# TODO: axis kwarg
def mean(values: np.ndarray, mask: np.ndarray, skipna: bool = True):
    if not values.size or mask.all():
        return libmissing.NA
    _sum = _sumprod(np.sum, values=values, mask=mask, skipna=skipna)
    count = np.count_nonzero(~mask)
    mean_value = _sum / count
    return mean_value
