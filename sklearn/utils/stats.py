# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..utils._array_api import (
    _cumsum,
    _find_matching_floating_dtype,
    _nextafter,
    _take_along_axis,
    get_namespace_and_device,
)


def _weighted_percentile(array, sample_weight, percentile=50):
    """Compute weighted percentile

    Computes lower weighted percentile. If `array` is a 2D array, the
    `percentile` is computed along the axis 0.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

    Parameters
    ----------
    array : 1D or 2D array
        Values to take the weighted percentile of.

    sample_weight: 1D or 2D array
        Weights for each value in `array`. Must be same shape as `array` or
        of shape `(array.shape[0],)`.

    percentile: int or float, default=50
        Percentile to compute. Must be value between 0 and 100.

    Returns
    -------
    percentile : int if `array` 1D, ndarray if `array` 2D
        Weighted percentile.
    """
    xp, is_array_api_compliant, device = get_namespace_and_device(array)
    sample_weight = xp.asarray(
        sample_weight,
        dtype=_find_matching_floating_dtype(sample_weight, xp=xp),
        device=device,
    )
    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = xp.reshape(array, (-1, 1))
    # When sample_weight 1D, repeat for each array.shape[1]
    if array.shape != sample_weight.shape and array.shape[0] == sample_weight.shape[0]:
        sample_weight = xp.tile(sample_weight, (array.shape[1], 1)).T
    sorted_idx = xp.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, xp)

    # Find index of median prediction for each sample
    weight_cdf = _cumsum(sorted_weights, axis=0)
    adjusted_percentile = percentile / 100 * weight_cdf[-1, :]
    weight_cdf = xp.asarray(
        weight_cdf, dtype=_find_matching_floating_dtype(weight_cdf, xp=xp)
    )
    adjusted_percentile = xp.asarray(
        adjusted_percentile,
        dtype=_find_matching_floating_dtype(adjusted_percentile, xp=xp),
    )

    # For percentile=0, ignore leading observations with sample_weight=0. GH20528
    mask = adjusted_percentile == 0
    adjusted_percentile[mask] = _nextafter(
        adjusted_percentile[mask], adjusted_percentile[mask] + 1, xp=xp
    )

    if adjusted_percentile.ndim == 0:
        percentile_idx = xp.asarray(
            [xp.searchsorted(weight_cdf[:, 0], adjusted_percentile)]
        )
    else:
        percentile_idx = xp.asarray(
            [
                xp.searchsorted(weight_cdf[:, i], adjusted_percentile[i])
                for i in range(weight_cdf.shape[1])
            ]
        )

    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = xp.clip(percentile_idx, 0, max_idx)

    col_index = xp.arange(array.shape[1])
    # percentile_in_sorted = sorted_idx[percentile_idx, col_index]
    percentile_in_sorted = []
    for i in range(percentile_idx.shape[0]):
        percentile_in_sorted.append(sorted_idx[percentile_idx[i], col_index[i]])
    percentile_in_sorted = xp.asarray(percentile_in_sorted)
    array = xp.asarray(array)
    # percentile = array[percentile_in_sorted, col_index]
    percentile = []
    for i in range(percentile_in_sorted.shape[0]):
        percentile.append(array[percentile_in_sorted[i], col_index[i]])
    percentile = xp.asarray(percentile)
    return percentile[0] if n_dim == 1 else percentile
