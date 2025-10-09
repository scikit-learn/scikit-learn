# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._array_api import (
    _find_matching_floating_dtype,
    get_namespace_and_device,
)


def _weighted_percentile(
    array, sample_weight, percentile_rank=50, average=False, xp=None
):
    """Compute the weighted percentile.

    Implement an array API compatible (weighted version) of NumPy's 'inverted_cdf'
    method when `average=False` (default) and 'averaged_inverted_cdf' when
    `average=True`.

    For an array ordered by increasing values, when the percentile lies exactly on a
    data point:

    * 'inverted_cdf' takes the exact data point.
    * 'averaged_inverted_cdf' takes the average of the exact data point and the one
      above it (this means it gives the same result as `median` for unit weights).

    E.g., for the array [1, 2, 3, 4] the percentile rank at each data point would
    be [25, 50, 75, 100]. Percentile rank 50 lies on '2'. 'average_inverted_cdf'
    computes the average of '2' and '3', making it 'symmetrical' because if you
    reverse the array, rank 50 would fall on '3'. It also matches 'median'.
    On the other hand, 'inverted_cdf', which does not satisfy the symmetry property,
    would give '2'.

    When the requested percentile lies between two data points, both methods return
    the higher data point.
    E.g., for the array [1, 2, 3, 4, 5] the percentile rank at each data point would
    be [20, 40, 60, 80, 100]. Percentile rank 50, lies between '2' and '3'. Taking the
    higher data point is symmetrical because if you reverse the array, 50 would lie
    between '4' and '3'. Both methods match median in this case.

    If `array` is a 2D array, the `values` are selected along axis 0.

    `NaN` values are ignored by setting their weights to 0. If `array` is 2D, this
    is done in a column-isolated manner: a `NaN` in the second column, does not impact
    the percentile computed for the first column even if `sample_weight` is 1D.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

        .. versionchanged:: 1.7
            Supports handling of `NaN` values.

        .. versionchanged:: 1.8
            Supports `average`, which calculates percentile using the
            "averaged_inverted_cdf" method.

    Parameters
    ----------
    array : 1D or 2D array
        Values to take the weighted percentile of.

    sample_weight: 1D or 2D array
        Weights for each value in `array`. Must be same shape as `array` or of shape
        `(array.shape[0],)`.

    percentile_rank: scalar or 1D array, default=50
        The probability level(s) of the percentile(s) to compute, in percent. Must be
        between 0 and 100. If a 1D array, computes multiple percentiles.

    average : bool, default=False
        If `True`, uses the "averaged_inverted_cdf" quantile method, otherwise
        defaults to "inverted_cdf". "averaged_inverted_cdf" is symmetrical with
        unit `sample_weight`, such that the total of `sample_weight` below or equal to
        `_weighted_percentile(percentile_rank)` is the same as the total of
        `sample_weight` above or equal to `_weighted_percentile(100-percentile_rank)`.
        This symmetry is not guaranteed with non-unit weights.

    xp : array_namespace, default=None
        The standard-compatible namespace for `array`. Default: infer.

    Returns
    -------
    percentile : scalar, 1D array, or 2D array
        Weighted percentile at the requested probability level(s).
        If `array` is 1D and `percentile_rank` is scalar, returns a scalar.
        If `array` is 2D and `percentile_rank` is scalar, returns a 1D array
            of shape `(array.shape[1],)`
        If `array` is 1D and `percentile_rank` is 1D, returns a 1D array
            of shape `(percentile_rank.shape[0],)`
        If `array` is 2D and `percentile_rank` is 1D, returns a 2D array
            of shape `(array.shape[1], percentile_rank.shape[0])`
    """
    xp, _, device = get_namespace_and_device(array)
    # `sample_weight` should follow `array` for dtypes
    floating_dtype = _find_matching_floating_dtype(array, xp=xp)
    array = xp.asarray(array, dtype=floating_dtype, device=device)
    sample_weight = xp.asarray(sample_weight, dtype=floating_dtype, device=device)
    percentile_rank = xp.asarray(percentile_rank, dtype=floating_dtype, device=device)

    n_dim = array.ndim
    if n_dim == 0:
        return array
    if array.ndim == 1:
        array = xp.reshape(array, (-1, 1))

    n_features = array.shape[1]
    n_dim_percentile = percentile_rank.ndim
    if n_dim_percentile == 0:
        percentile_rank = xp.reshape(percentile_rank, (1,))
    if xp.any(percentile_rank[1:] < percentile_rank[:-1]):
        raise ValueError("percentile_rank must be non-decreasing")

    q = percentile_rank / 100
    n_percentiles = percentile_rank.shape[0]
    result = xp.empty((n_features, n_percentiles), dtype=floating_dtype)

    # Compute weighted percentiles for each feature (column)
    for feature_idx in range(n_features):
        x = array[..., feature_idx]
        # Ignore NaN values by masking them out
        mask_not_nan = ~xp.isnan(x)
        x = x[mask_not_nan]
        if x.shape[0] == 0:
            # If all values are NaN, return NaN for this feature
            result[feature_idx, ...] = xp.nan
            continue
        # Select weights for non-NaN values
        w = (
            sample_weight[..., feature_idx][mask_not_nan]
            if sample_weight.ndim == 2
            else sample_weight[mask_not_nan]
        )
        # Ignore zero weights
        mask_nz = w != 0
        has_zero = not xp.all(mask_nz)
        if has_zero:
            w = w[mask_nz]
        weights_sum = xp.sum(w)
        if weights_sum == 0:
            # If all weights are zero, return max value (consistent with NaN handling)
            result[feature_idx, ...] = xp.max(x)
            continue
        if has_zero:
            x = x[mask_nz]
        # Recursively compute weighted percentiles using partitioning
        w_sorted = False
        if not hasattr(xp, "argpartition"):
            x_sorter = xp.argsort(x, stable=False)
            w = w[x_sorter]
            x = x[x_sorter]
            w_sorted = True
        _weighted_percentile_inner(
            x,
            w,
            target_sums=weights_sum * q,
            out=result[feature_idx, ...],
            average=average,
            w_sorted=w_sorted,
            xp=xp,
        )

    if n_dim_percentile == 0:
        result = result[..., 0]

    return result[0] if n_dim == 1 else result


def _weighted_percentile_inner(x, w, target_sums, out, average, w_sorted, xp):
    n = x.shape[0]
    if n == 1:
        out[:] = x
        return
    i = n // 2
    if w_sorted:
        w_left = w[:i]
        x_left = x[:i]
        w_right = w[i:]
        x_right = x[i:]
    else:
        partitioner = xp.argpartition(x, i)
        w_left = w[partitioner[:i]]
        x_left = w_right = x_right = None
    sum_left = xp.sum(w_left)
    j = xp.searchsorted(target_sums, sum_left)
    target_sums[j:] -= sum_left
    if j > 0:
        # some quantiles are to be found on the left side of the partition
        x_left = x[partitioner[:i]] if x_left is None else x_left
        _weighted_percentile_inner(
            x_left, w_left, target_sums[:j], out[:j], average, w_sorted, xp
        )
    if j >= target_sums.shape[0]:
        return
    zero = xp.sum(target_sums[:0])
    # ^ array API needs an array as argument for searchsorted
    idx_0 = xp.searchsorted(target_sums[j:], zero, side="right")
    if idx_0 > 0:
        # some quantiles are precisely at the index of the partition
        x_left = x[partitioner[:i]] if x_left is None else x_left
        x_right = x[partitioner[i:]] if x_right is None else x_right
        out[j : j + idx_0] = (
            (xp.max(x_left) + xp.min(x_right)) / 2 if average else xp.max(x_left)
        )
        j += idx_0
    if j < target_sums.shape[0]:
        # some quantiles are to be found on the right side of the partition
        x_right = x[partitioner[i:]] if x_right is None else x_right
        w_right = w[partitioner[i:]] if w_right is None else w_right
        _weighted_percentile_inner(
            x_right, w_right, target_sums[j:], out[j:], average, w_sorted, xp
        )


def _old_weighted_percentile(
    array, sample_weight, percentile_rank=50, average=False, xp=None
):
    """Compute the weighted percentile.

    Implement an array API compatible (weighted version) of NumPy's 'inverted_cdf'
    method when `average=False` (default) and 'averaged_inverted_cdf' when
    `average=True`.

    For an array ordered by increasing values, when the percentile lies exactly on a
    data point:

    * 'inverted_cdf' takes the exact data point.
    * 'averaged_inverted_cdf' takes the average of the exact data point and the one
      above it (this means it gives the same result as `median` for unit weights).

    E.g., for the array [1, 2, 3, 4] the percentile rank at each data point would
    be [25, 50, 75, 100]. Percentile rank 50 lies on '2'. 'average_inverted_cdf'
    computes the average of '2' and '3', making it 'symmetrical' because if you
    reverse the array, rank 50 would fall on '3'. It also matches 'median'.
    On the other hand, 'inverted_cdf', which does not satisfy the symmetry property,
    would give '2'.

    When the requested percentile lies between two data points, both methods return
    the higher data point.
    E.g., for the array [1, 2, 3, 4, 5] the percentile rank at each data point would
    be [20, 40, 60, 80, 100]. Percentile rank 50, lies between '2' and '3'. Taking the
    higher data point is symmetrical because if you reverse the array, 50 would lie
    between '4' and '3'. Both methods match median in this case.

    If `array` is a 2D array, the `values` are selected along axis 0.

    `NaN` values are ignored by setting their weights to 0. If `array` is 2D, this
    is done in a column-isolated manner: a `NaN` in the second column, does not impact
    the percentile computed for the first column even if `sample_weight` is 1D.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

        .. versionchanged:: 1.7
            Supports handling of `NaN` values.

        .. versionchanged:: 1.8
            Supports `average`, which calculates percentile using the
            "averaged_inverted_cdf" method.

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

    average : bool, default=False
        If `True`, uses the "averaged_inverted_cdf" quantile method, otherwise
        defaults to "inverted_cdf". "averaged_inverted_cdf" is symmetrical with
        unit `sample_weight`, such that the total of `sample_weight` below or equal to
        `_weighted_percentile(percentile_rank)` is the same as the total of
        `sample_weight` above or equal to `_weighted_percentile(100-percentile_rank)`.
        This symmetry is not guaranteed with non-unit weights.

    xp : array_namespace, default=None
        The standard-compatible namespace for `array`. Default: infer.

    Returns
    -------
    percentile : scalar or 0D array if `array` 1D (or 0D), array if `array` 2D
        Weighted percentile at the requested probability level.
    """
    xp, _, device = get_namespace_and_device(array)
    # `sample_weight` should follow `array` for dtypes
    floating_dtype = _find_matching_floating_dtype(array, xp=xp)
    array = xp.asarray(array, dtype=floating_dtype, device=device)
    sample_weight = xp.asarray(sample_weight, dtype=floating_dtype, device=device)

    n_dim = array.ndim
    if n_dim == 0:
        return array
    if array.ndim == 1:
        array = xp.reshape(array, (-1, 1))
    # When sample_weight 1D, repeat for each array.shape[1]
    if array.shape != sample_weight.shape and array.shape[0] == sample_weight.shape[0]:
        sample_weight = xp.tile(sample_weight, (array.shape[1], 1)).T
    # Sort `array` and `sample_weight` along axis=0:
    sorted_idx = xp.argsort(array, axis=0, stable=False)
    sorted_weights = xp.take_along_axis(sample_weight, sorted_idx, axis=0)

    # Set NaN values in `sample_weight` to 0. Only perform this operation if NaN
    # values present to avoid temporary allocations of size `(n_samples, n_features)`.
    n_features = array.shape[1]
    largest_value_per_column = array[
        sorted_idx[-1, ...], xp.arange(n_features, device=device)
    ]
    # NaN values get sorted to end (largest value)
    if xp.any(xp.isnan(largest_value_per_column)):
        sorted_nan_mask = xp.take_along_axis(xp.isnan(array), sorted_idx, axis=0)
        sorted_weights[sorted_nan_mask] = 0

    # Compute the weighted cumulative distribution function (CDF) based on
    # `sample_weight` and scale `percentile_rank` along it.
    #
    # Note: we call `xp.cumulative_sum` on the transposed `sorted_weights` to
    # ensure that the result is of shape `(n_features, n_samples)` so
    # `xp.searchsorted` calls take contiguous inputs as a result (for
    # performance reasons).
    weight_cdf = xp.cumulative_sum(sorted_weights.T, axis=1)
    adjusted_percentile_rank = percentile_rank / 100 * weight_cdf[..., -1]

    # Ignore leading `sample_weight=0` observations when `percentile_rank=0` (#20528)
    mask = adjusted_percentile_rank == 0
    adjusted_percentile_rank[mask] = xp.nextafter(
        adjusted_percentile_rank[mask], adjusted_percentile_rank[mask] + 1
    )
    # For each feature with index j, find sample index i of the scalar value
    # `adjusted_percentile_rank[j]` in 1D array `weight_cdf[j]`, such that:
    # weight_cdf[j, i-1] < adjusted_percentile_rank[j] <= weight_cdf[j, i].
    # Note `searchsorted` defaults to equality on the right, whereas Hyndman and Fan
    # reference equation has equality on the left.
    percentile_indices = xp.stack(
        [
            xp.searchsorted(
                weight_cdf[feature_idx, ...], adjusted_percentile_rank[feature_idx]
            )
            for feature_idx in range(weight_cdf.shape[0])
        ],
    )
    # `percentile_indices` may be equal to `sorted_idx.shape[0]` due to floating
    # point error (see #11813)
    max_idx = sorted_idx.shape[0] - 1
    percentile_indices = xp.clip(percentile_indices, 0, max_idx)

    col_indices = xp.arange(array.shape[1], device=device)
    percentile_in_sorted = sorted_idx[percentile_indices, col_indices]

    if average:
        # From Hyndman and Fan (1996), `fraction_above` is `g`
        fraction_above = (
            weight_cdf[col_indices, percentile_indices] - adjusted_percentile_rank
        )
        is_fraction_above = fraction_above > xp.finfo(floating_dtype).eps
        percentile_plus_one_indices = xp.clip(percentile_indices + 1, 0, max_idx)
        percentile_plus_one_in_sorted = sorted_idx[
            percentile_plus_one_indices, col_indices
        ]
        # Handle case when next index ('plus one') has sample weight of 0
        zero_weight_cols = col_indices[
            sample_weight[percentile_plus_one_in_sorted, col_indices] == 0
        ]
        for col_idx in zero_weight_cols:
            cdf_val = weight_cdf[col_idx, percentile_indices[col_idx]]
            # Search for next index where `weighted_cdf` is greater
            next_index = xp.searchsorted(
                weight_cdf[col_idx, ...], cdf_val, side="right"
            )
            # Handle case where there are trailing 0 sample weight samples
            # and `percentile_indices` is already max index
            if next_index >= max_idx:
                # use original `percentile_indices` again
                next_index = percentile_indices[col_idx]

            percentile_plus_one_in_sorted[col_idx] = sorted_idx[next_index, col_idx]

        result = xp.where(
            is_fraction_above,
            array[percentile_in_sorted, col_indices],
            (
                array[percentile_in_sorted, col_indices]
                + array[percentile_plus_one_in_sorted, col_indices]
            )
            / 2,
        )
    else:
        result = array[percentile_in_sorted, col_indices]

    return result[0] if n_dim == 1 else result
