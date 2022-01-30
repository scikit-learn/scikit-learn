"""
Numba 1D sum kernels that can be shared by
* Dataframe / Series
* groupby
* rolling / expanding

Mirrors pandas/_libs/window/aggregation.pyx
"""
from __future__ import annotations

import numba
import numpy as np

from pandas.core._numba.kernels.shared import is_monotonic_increasing


@numba.jit(nopython=True, nogil=True, parallel=False)
def add_sum(
    val: float, nobs: int, sum_x: float, compensation: float
) -> tuple[int, float, float]:
    if not np.isnan(val):
        nobs += 1
        y = val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t
    return nobs, sum_x, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def remove_sum(
    val: float, nobs: int, sum_x: float, compensation: float
) -> tuple[int, float, float]:
    if not np.isnan(val):
        nobs -= 1
        y = -val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t
    return nobs, sum_x, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def sliding_sum(
    values: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    min_periods: int,
) -> np.ndarray:
    N = len(start)
    nobs = 0
    sum_x = 0.0
    compensation_add = 0.0
    compensation_remove = 0.0

    is_monotonic_increasing_bounds = is_monotonic_increasing(
        start
    ) and is_monotonic_increasing(end)

    output = np.empty(N, dtype=np.float64)

    for i in range(N):
        s = start[i]
        e = end[i]
        if i == 0 or not is_monotonic_increasing_bounds:
            for j in range(s, e):
                val = values[j]
                nobs, sum_x, compensation_add = add_sum(
                    val, nobs, sum_x, compensation_add
                )
        else:
            for j in range(start[i - 1], s):
                val = values[j]
                nobs, sum_x, compensation_remove = remove_sum(
                    val, nobs, sum_x, compensation_remove
                )

            for j in range(end[i - 1], e):
                val = values[j]
                nobs, sum_x, compensation_add = add_sum(
                    val, nobs, sum_x, compensation_add
                )

        if nobs == 0 == nobs:
            result = 0.0
        elif nobs >= min_periods:
            result = sum_x
        else:
            result = np.nan

        output[i] = result

        if not is_monotonic_increasing_bounds:
            nobs = 0
            sum_x = 0.0
            compensation_remove = 0.0

    return output
