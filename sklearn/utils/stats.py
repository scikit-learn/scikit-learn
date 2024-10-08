# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from .extmath import stable_cumsum


def _weighted_percentile(array, sample_weight, percentile=50):
    """Compute the lower value at a given weighted percentile.

    If `array` is a 2D array, the `value` is selected along axis 0.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

        .. versionchanged:: 1.6
            Supports handling of `NaN` values. For `NaN` inputs, their corresponding
            weights are set to 0 in `sample_weight`.
            Percentiles that point to `NaN` values are redirected to the next lower
            value if it exists.

    Parameters
    ----------
    array : 1D or 2D array
        Values to take the weighted percentile of.

    sample_weight: 1D or 2D array
        Weights for each value in `array`. Must be same shape as `array` or of shape
        `(array.shape[0],)`.

    percentile: int or float, default=50
        Percentile to compute. Must be value between 0 and 100.

    Returns
    -------
    value : int if `array` 1D, ndarray if `array` 2D
        Lower value at a given weighted percentile.
    """
    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))
    # When sample_weight 1D, repeat for each array.shape[1]
    if array.shape != sample_weight.shape and array.shape[0] == sample_weight.shape[0]:
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T

    # Sort `array` and `sample_weight` along axis=0; set `sample_weight` for nan input
    # to 0:
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = np.take_along_axis(sample_weight, sorted_idx, axis=0)
    sorted_nan_mask = np.take_along_axis(np.isnan(array), sorted_idx, axis=0)
    sorted_weights[sorted_nan_mask] = 0

    # Compute the weighted cumulative distribution function (CDF) based on
    # `sample_weight` and scale `percentile` along it:
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    adjusted_percentile = percentile / 100 * weight_cdf[-1]

    # For percentile=0, ignore leading observations with sample_weight=0; see PR #20528
    mask = adjusted_percentile == 0
    adjusted_percentile[mask] = np.nextafter(
        adjusted_percentile[mask], adjusted_percentile[mask] + 1
    )

    # Find index `adjusted_percentile` would have in `weight_cdf`:
    percentile_idx = np.array(
        [
            np.searchsorted(weight_cdf[:, i], adjusted_percentile[i])
            for i in range(weight_cdf.shape[1])
        ]
    )

    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(
        lambda x: np.clip(x, 0, max_idx), axis=0, arr=percentile_idx
    )

    col_index = np.arange(array.shape[1])
    percentile_in_sorted = sorted_idx[percentile_idx, col_index]
    value = array[percentile_in_sorted, col_index]

    # Percentiles that point to nan values are redirected to the next lower
    # value unless we have reached the lowest index (0) in `sortex_idx`:
    while (percentile_isnan_mask := np.isnan(value)).any() and (
        percentile_idx[percentile_isnan_mask] > 0
    ).any():
        percentile_idx[percentile_isnan_mask] = np.maximum(
            percentile_idx[percentile_isnan_mask] - 1, 0
        )
        percentile_in_sorted = sorted_idx[percentile_idx, col_index]
        value = array[percentile_in_sorted, col_index]

    return value[0] if n_dim == 1 else value
