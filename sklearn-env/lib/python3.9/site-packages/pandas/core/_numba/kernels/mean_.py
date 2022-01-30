"""
Numba 1D mean kernels that can be shared by
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
def add_mean(
    val: float, nobs: int, sum_x: float, neg_ct: int, compensation: float
) -> tuple[int, float, int, float]:
    if not np.isnan(val):
        nobs += 1
        y = val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t
        if val < 0:
            neg_ct += 1
    return nobs, sum_x, neg_ct, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def remove_mean(
    val: float, nobs: int, sum_x: float, neg_ct: int, compensation: float
) -> tuple[int, float, int, float]:
    if not np.isnan(val):
        nobs -= 1
        y = -val - compensation
        t = sum_x + y
        compensation = t - sum_x - y
        sum_x = t
        if val < 0:
            neg_ct -= 1
    return nobs, sum_x, neg_ct, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def sliding_mean(
    values: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    min_periods: int,
) -> np.ndarray:
    N = len(start)
    nobs = 0
    sum_x = 0.0
    neg_ct = 0
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
                nobs, sum_x, neg_ct, compensation_add = add_mean(
                    val, nobs, sum_x, neg_ct, compensation_add
                )
        else:
            for j in range(start[i - 1], s):
                val = values[j]
                nobs, sum_x, neg_ct, compensation_remove = remove_mean(
                    val, nobs, sum_x, neg_ct, compensation_remove
                )

            for j in range(end[i - 1], e):
                val = values[j]
                nobs, sum_x, neg_ct, compensation_add = add_mean(
                    val, nobs, sum_x, neg_ct, compensation_add
                )

        if nobs >= min_periods and nobs > 0:
            result = sum_x / nobs
            if neg_ct == 0 and result < 0:
                result = 0
            elif neg_ct == nobs and result > 0:
                result = 0
        else:
            result = np.nan

        output[i] = result

        if not is_monotonic_increasing_bounds:
            nobs = 0
            sum_x = 0.0
            neg_ct = 0
            compensation_remove = 0.0

    return output
