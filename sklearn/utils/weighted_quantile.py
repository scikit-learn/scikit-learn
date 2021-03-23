"""
authors: Jasper Roebroek <roebroek.jasper@gmail.com>

The calculation is roughly 10 times as slow as np.quantile, which
is not terrible as the data needs to be copied and sorted.
"""

import numpy as np


def weighted_quantile(a, q, weights, axis=-1):
    """
    Returns the weighted quantile on a

    Parameters
    ----------
    a: array-like
        Data on which the quantiles are calculated

    q: float
        Quantile (in the range from 0-1)

    weights: array-like, optional
        Weights corresponding to a

    axis : int, optional
        Axis over which the quantile values are calculated. By default the
        last axis is used.

    Returns
    -------
    quantile: array

    References
    ----------
    1. https://en.wikipedia.org/wiki/Percentile#The_Weighted_Percentile_method
    """
    if q > 1 or q < 0:
        raise ValueError("q should be in-between 0 and 1, "
                         "got %d" % q)

    if weights is None:
        return np.quantile(a, q, axis=-1)
    else:
        a = np.asarray(a, dtype=np.float64)
        weights = np.asarray(weights)
        if a.shape[:-1] == weights.shape:
            return np.quantile(a, q, axis=axis)
        elif a.shape != weights.shape:
            raise IndexError("the data and weights need to be of the same shape")

    a = a.copy()
    zeros = weights == 0
    a[zeros] = np.nan
    zeros_count = zeros.sum(axis=axis, keepdims=True)

    idx_sorted = np.argsort(a, axis=axis)
    a_sorted = np.take_along_axis(a, idx_sorted, axis=axis)
    weights_sorted = np.take_along_axis(weights, idx_sorted, axis=axis)

    # Step 1
    weights_cum = np.cumsum(weights_sorted, axis=axis)
    weights_total = np.expand_dims(np.take(weights_cum, -1, axis=axis), axis=axis)

    # Step 2
    weights_norm = (weights_cum - 0.5 * weights_sorted) / weights_total
    start = np.sum(weights_norm < q, axis=axis, keepdims=True) - 1

    idx_low = (start == -1).squeeze()
    high = a.shape[axis] - zeros_count - 1
    idx_high = (start == high).squeeze()

    start = np.clip(start, 0, high - 1)

    # Step 3.
    left_weight = np.take_along_axis(weights_norm, start, axis=axis)
    right_weight = np.take_along_axis(weights_norm, start + 1, axis=axis)
    left_value = np.take_along_axis(a_sorted, start, axis=axis)
    right_value = np.take_along_axis(a_sorted, start + 1, axis=axis)

    fraction = (q - left_weight) / (right_weight - left_weight)
    quantiles = left_value + fraction * (right_value - left_value)

    if idx_low.sum() > 0:
        quantiles[idx_low] = np.take(a_sorted[idx_low], 0, axis=axis)
    if idx_high.sum() > 0:
        quantiles[idx_high] = np.take(a_sorted[idx_high], a.shape[axis] - zeros_count - 1, axis=axis)

    return quantiles.squeeze()
