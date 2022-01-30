"""
Numba 1D var kernels that can be shared by
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
def add_var(
    val: float, nobs: int, mean_x: float, ssqdm_x: float, compensation: float
) -> tuple[int, float, float, float]:
    if not np.isnan(val):
        nobs += 1
        prev_mean = mean_x - compensation
        y = val - compensation
        t = y - mean_x
        compensation = t + mean_x - y
        delta = t
        if nobs:
            mean_x += delta / nobs
        else:
            mean_x = 0
        ssqdm_x += (val - prev_mean) * (val - mean_x)
    return nobs, mean_x, ssqdm_x, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def remove_var(
    val: float, nobs: int, mean_x: float, ssqdm_x: float, compensation: float
) -> tuple[int, float, float, float]:
    if not np.isnan(val):
        nobs -= 1
        if nobs:
            prev_mean = mean_x - compensation
            y = val - compensation
            t = y - mean_x
            compensation = t + mean_x - y
            delta = t
            mean_x -= delta / nobs
            ssqdm_x -= (val - prev_mean) * (val - mean_x)
        else:
            mean_x = 0
            ssqdm_x = 0
    return nobs, mean_x, ssqdm_x, compensation


@numba.jit(nopython=True, nogil=True, parallel=False)
def sliding_var(
    values: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    min_periods: int,
    ddof: int = 1,
) -> np.ndarray:
    N = len(start)
    nobs = 0
    mean_x = 0.0
    ssqdm_x = 0.0
    compensation_add = 0.0
    compensation_remove = 0.0

    min_periods = max(min_periods, 1)
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
                nobs, mean_x, ssqdm_x, compensation_add = add_var(
                    val, nobs, mean_x, ssqdm_x, compensation_add
                )
        else:
            for j in range(start[i - 1], s):
                val = values[j]
                nobs, mean_x, ssqdm_x, compensation_remove = remove_var(
                    val, nobs, mean_x, ssqdm_x, compensation_remove
                )

            for j in range(end[i - 1], e):
                val = values[j]
                nobs, mean_x, ssqdm_x, compensation_add = add_var(
                    val, nobs, mean_x, ssqdm_x, compensation_add
                )

        if nobs >= min_periods and nobs > ddof:
            if nobs == 1:
                result = 0.0
            else:
                result = ssqdm_x / (nobs - ddof)
        else:
            result = np.nan

        output[i] = result

        if not is_monotonic_increasing_bounds:
            nobs = 0
            mean_x = 0.0
            ssqdm_x = 0.0
            compensation_remove = 0.0

    return output
