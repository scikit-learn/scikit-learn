# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..utils._array_api import (
    _find_matching_floating_dtype,
    get_namespace_and_device,
)


def _weighted_percentile(
    array, sample_weight, percentile_rank=50, xp=None, symmetrize=False
):
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

    xp : array_namespace, default=None
        The standard-compatible namespace for `array`. Default: infer.

    symmetrize : bool, default=False
        If True, compute the symmetrized weighted percentile using the negative of
        `array` and the fact that the sorted order of the negated array is the
        reverse order of the sorted array.

    Returns
    -------
    percentile : scalar or 0D array if `array` 1D (or 0D), array if `array` 2D
        Weighted percentile at the requested probability level.
    """
    xp, _, device = get_namespace_and_device(array)
    floating_dtype = _find_matching_floating_dtype(array, xp=xp)
    array = xp.asarray(array, dtype=floating_dtype, device=device)
    sample_weight = xp.asarray(sample_weight, dtype=floating_dtype, device=device)

    n_dim = array.ndim
    if n_dim == 0:
        return array
    if array.ndim == 1:
        array = xp.reshape(array, (-1, 1))
    if array.shape != sample_weight.shape and array.shape[0] == sample_weight.shape[0]:
        sample_weight = xp.tile(sample_weight, (array.shape[1], 1)).T

    # Sort `array` and `sample_weight` along axis=0:
    sorted_idx = xp.argsort(array, axis=0)
    sorted_weights = xp.take_along_axis(sample_weight, sorted_idx, axis=0)

    n_features = array.shape[1]
    largest_value_per_column = array[
        sorted_idx[-1, ...], xp.arange(n_features, device=device)
    ]
    # NaN values get sorted to end (largest value)
    if xp.any(xp.isnan(largest_value_per_column)):
        sorted_nan_mask = xp.take_along_axis(xp.isnan(array), sorted_idx, axis=0)
        sorted_weights[sorted_nan_mask] = 0

    weight_cdf = xp.cumulative_sum(sorted_weights.T, axis=1)
    adjusted_percentile_rank = percentile_rank / 100 * weight_cdf[..., -1]
    mask = adjusted_percentile_rank == 0
    adjusted_percentile_rank[mask] = xp.nextafter(
        adjusted_percentile_rank[mask], adjusted_percentile_rank[mask] + 1
    )
    percentile_indices = xp.asarray(
        [
            xp.searchsorted(
                weight_cdf[feature_idx, ...], adjusted_percentile_rank[feature_idx]
            )
            for feature_idx in range(weight_cdf.shape[0])
        ],
        device=device,
    )
    max_idx = sorted_idx.shape[0] - 1
    percentile_indices = xp.clip(percentile_indices, 0, max_idx)
    col_indices = xp.arange(array.shape[1], device=device)
    percentile_in_sorted = sorted_idx[percentile_indices, col_indices]
    result = array[percentile_in_sorted, col_indices]

    if symmetrize:
        # Use reverse of the original sorted order for the negated array
        sorted_idx_neg = xp.flip(sorted_idx, axis=0)
        sorted_weights_neg = xp.take_along_axis(sample_weight, sorted_idx_neg, axis=0)
        array_neg = -array
        largest_value_neg = array_neg[
            sorted_idx_neg[-1, ...], xp.arange(n_features, device=device)
        ]
        if xp.any(xp.isnan(largest_value_neg)):
            sorted_nan_mask_neg = xp.take_along_axis(
                xp.isnan(array_neg), sorted_idx_neg, axis=0
            )
            sorted_weights_neg[sorted_nan_mask_neg] = 0

        weight_cdf_neg = xp.cumulative_sum(sorted_weights_neg.T, axis=1)
        adjusted_pr_neg = (100 - percentile_rank) / 100 * weight_cdf_neg[..., -1]
        mask_neg = adjusted_pr_neg == 0
        adjusted_pr_neg[mask_neg] = xp.nextafter(
            adjusted_pr_neg[mask_neg], adjusted_pr_neg[mask_neg] + 1
        )
        percentile_indices_neg = xp.asarray(
            [
                xp.searchsorted(
                    weight_cdf_neg[feature_idx, ...], adjusted_pr_neg[feature_idx]
                )
                for feature_idx in range(weight_cdf_neg.shape[0])
            ],
            device=device,
        )
        percentile_indices_neg = xp.clip(percentile_indices_neg, 0, max_idx)
        percentile_in_sorted_neg = sorted_idx_neg[percentile_indices_neg, col_indices]
        result_neg = array_neg[percentile_in_sorted_neg, col_indices]

        result = (result - result_neg) / 2

    return result[0] if n_dim == 1 else result


def _averaged_weighted_percentile(array, sample_weight, percentile_rank=50, xp=None):
    return _weighted_percentile(
        array, sample_weight, percentile_rank, xp=xp, symmetrize=True
    )
