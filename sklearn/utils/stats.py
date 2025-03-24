# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from .extmath import stable_cumsum


def _weighted_percentile(array, sample_weight, percentile_rank=50):
    """Compute the weighted percentile with method 'inverted_cdf'.

    When the percentile lies between two data points of `array`, the function returns
    the lower value.

    If `array` is a 2D array, the `values` are selected along axis 0.

    `NaN` values are ignored by setting their weights to 0. If `array` is 2D, this
    is done in a column-isolated manner: a `NaN` in the second column, does not impact
    the percentile computed for the first column even if `sample_weight` is 1D.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

        .. versionchanged:: 1.7
            Supports handling of `NaN` values.

    Parameters
    ----------
    array : 1D or 2D array
        Values to take the weighted percentile of.

    sample_weight: 1D or 2D array
        Weights for each value in `array`. Must be same shape as `array` or of shape
        `(array.shape[0],)`.

    percentile_rank: int or float, default=50
        The probability level of the percentile to compute, in percent. Must be between
        0 and 100.

    Returns
    -------
    percentile : int if `array` 1D, ndarray if `array` 2D
        Weighted percentile at the requested probability level.
    """
    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))
    # When sample_weight 1D, repeat for each array.shape[1]
    if array.shape != sample_weight.shape and array.shape[0] == sample_weight.shape[0]:
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T

    # Sort `array` and `sample_weight` along axis=0:
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = np.take_along_axis(sample_weight, sorted_idx, axis=0)

    # Set NaN values in `sample_weight` to 0. We only perform this operation if NaN
    # values are present at all to avoid temporary allocations of size `(n_samples,
    # n_features)`. If NaN values were present, they would sort to the end (which we can
    # observe from `sorted_idx`).
    n_features = array.shape[1]
    largest_value_per_column = array[sorted_idx[-1, ...], np.arange(n_features)]
    if np.isnan(largest_value_per_column).any():
        sorted_nan_mask = np.take_along_axis(np.isnan(array), sorted_idx, axis=0)
        sorted_weights[sorted_nan_mask] = 0

    # Compute the weighted cumulative distribution function (CDF) based on
    # sample_weight and scale percentile_rank along it:
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    adjusted_percentile_rank = percentile_rank / 100 * weight_cdf[-1]

    # For percentile_rank=0, ignore leading observations with sample_weight=0; see
    # PR #20528:
    mask = adjusted_percentile_rank == 0
    adjusted_percentile_rank[mask] = np.nextafter(
        adjusted_percentile_rank[mask], adjusted_percentile_rank[mask] + 1
    )

    # Find index (i) of `adjusted_percentile` in `weight_cdf`,
    # such that weight_cdf[i-1] < percentile <= weight_cdf[i]
    percentile_idx = np.array(
        [
            np.searchsorted(weight_cdf[:, i], adjusted_percentile_rank[i])
            for i in range(weight_cdf.shape[1])
        ]
    )

    # In rare cases, percentile_idx equals to sorted_idx.shape[0]:
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(
        lambda x: np.clip(x, 0, max_idx), axis=0, arr=percentile_idx
    )

    col_indices = np.arange(array.shape[1])
    percentile_in_sorted = sorted_idx[percentile_idx, col_indices]
    result = array[percentile_in_sorted, col_indices]

    return result[0] if n_dim == 1 else result


# TODO: refactor to do the symmetrisation inside _weighted_percentile to avoid
# sorting the input array twice.
def _averaged_weighted_percentile(array, sample_weight, percentile_rank=50):
    return (
        _weighted_percentile(array, sample_weight, percentile_rank)
        - _weighted_percentile(-array, sample_weight, 100 - percentile_rank)
    ) / 2
