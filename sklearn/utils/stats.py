from collections.abc import Iterable

import numpy as np

from .extmath import stable_cumsum
from .fixes import _take_along_axis


def _squeeze_arr(arr, n_dim):
    return arr[0] if n_dim == 1 else arr


def _nearest_rank(array, sorted_idx, percentile, cum_weights):
    adjusted_percentile = percentile * cum_weights[-1]
    percentile_idx = np.array([
        np.searchsorted(cum_weights[:, i], adjusted_percentile[i])
        for i in range(cum_weights.shape[1])
    ])
    percentile_idx = np.array(percentile_idx)
    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = sorted_idx.shape[0] - 1
    percentile_idx = np.apply_along_axis(
        lambda x: np.clip(x, 0, max_idx), axis=0, arr=percentile_idx
    )

    col_index = np.arange(array.shape[1])
    percentile_in_sorted = sorted_idx[percentile_idx, col_index]
    percentile_value = array[percentile_in_sorted, col_index]
    return _squeeze_arr(percentile_value, array.ndim)


def _interpolation_closest_ranks(
    array, sorted_idx, sample_weight, sorted_weights, percentile, cum_weights,
    interpolation,
):
    # Percentile can be computed with 3 different alternatives:
    # https://en.wikipedia.org/wiki/Percentile
    # These 3 alternatives depend of the value of a parameter C. NumPy uses
    # the variant where C=0 which allows to obtain a strictly monotonically
    # increasing function which is defined as:
    # P = (x - 1) / (N - 1); x in [1, N]
    # Weighted percentile change this formula by taking into account the
    # weights instead of the data frequency.
    # P_w = (x - w) / (S_w - w), x in [1, N], w being the weight and S_w being
    # the sum of the weights.
    n_rows, n_cols = array.shape
    n_dim = array.ndim

    adjusted_percentile = (cum_weights - sorted_weights)
    with np.errstate(invalid="ignore"):
        adjusted_percentile /= cum_weights[-1] - sorted_weights
        nan_mask = np.isnan(adjusted_percentile)
        adjusted_percentile[nan_mask] = 1

    if interpolation in ("lower", "higher", "nearest"):
        percentile_idx = np.array([
            np.searchsorted(adjusted_percentile[:, col], percentile[col],
                            side="left")
            for col in range(n_cols)
        ])

        if interpolation == "lower" and np.all(percentile < 1):
            # P = 100 is a corner case for "lower"
            percentile_idx -= 1
        elif interpolation == "nearest" and np.all(percentile < 1):
            for col in range(n_cols):
                error_higher = abs(
                    adjusted_percentile[percentile_idx[col], col] -
                    percentile[col]
                )
                error_lower = abs(
                    adjusted_percentile[percentile_idx[col] - 1, col] -
                    percentile[col]
                )
                if error_higher >= error_lower:
                    percentile_idx[col] -= 1

        percentile_idx = np.apply_along_axis(
            lambda x: np.clip(x, 0, n_rows - 1), axis=0,
            arr=percentile_idx
        )

        percentile_value = array[
            sorted_idx[percentile_idx, np.arange(n_cols)],
            np.arange(n_cols)
        ]
        percentile_value = _squeeze_arr(percentile_value, n_dim)

    else:  # interpolation == "linear"
        percentile_value = np.array([
            np.interp(
                x=percentile[col],
                xp=adjusted_percentile[:, col],
                fp=array[sorted_idx[:, col], col],
            )
            for col in range(n_cols)
        ])

        percentile_value = _squeeze_arr(percentile_value, n_dim)

    single_sample_weight = np.count_nonzero(sample_weight, axis=0)
    if np.any(single_sample_weight == 1):
        # edge case where a single weight is non-null in which case the
        # previous methods will fail
        if not isinstance(percentile_value, Iterable):
            percentile_value = _squeeze_arr(
                array[np.nonzero(sample_weight)], n_dim
            )
        else:
            percentile_value = np.array([
                array[np.flatnonzero(sample_weight[:, col])[0], col]
                if n_nonzero == 1 else percentile_value[col]
                for col, n_nonzero in enumerate(single_sample_weight)
            ])

    return percentile_value


def _weighted_percentile(array, sample_weight, percentile=50,
                         interpolation=None):
    """Compute weighted percentile

    Computes lower weighted percentile. If `array` is a 2D array, the
    `percentile` is computed along the axis 0.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

    Parameters
    ----------
    array : ndarray of shape (n,) or (n, m)
        Values to take the weighted percentile of.

    sample_weight: ndarray of shape (n,) or (n, m)
        Weights for each value in `array`. Must be same shape as `array` or
        of shape `(array.shape[0],)`.

    percentile: int or float, default=50
        Percentile to compute. Must be value between 0 and 100.

    interpolation : {"linear", "lower", "higher", "nearest"}, default=None
        The interpolation method to use when the percentile lies between
        data points `i` and `j`:

        * None: no interpolation will be done and the "nearest-rank" method
          will be used.
        * "linear": linearly interpolate between `i` and `j` using `np.interp`.
        * "lower": i`;
        * "higher": `j`;
        * "nearest": `i` or `j`, whichever is nearest (default).

        .. versionadded: 0.24

    Returns
    -------
    percentile_value : float or int if `array` of shape (n,), otherwise \
            ndarray of shape (m,)
        Weighted percentile.
    """
    possible_interpolation = ("linear", "lower", "higher", "nearest", None)
    if interpolation not in possible_interpolation:
        raise ValueError(
            f"'interpolation' should be one of "
            f"{', '.join(possible_interpolation)}. "
            f"Got '{str(interpolation)}' instead."
        )

    if np.any(np.count_nonzero(sample_weight, axis=0) < 1):
        raise ValueError(
            "All weights cannot be null when computing a weighted percentile."
        )

    if array.ndim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))

    if (array.shape != sample_weight.shape and
            array.shape[0] == sample_weight.shape[0]):
        # when `sample_weight` is 1D, we repeat it for each column of `array`
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T

    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)
    percentile = np.array([percentile / 100] * array.shape[1])
    cum_weights = stable_cumsum(sorted_weights, axis=0)

    if interpolation is None:
        return _nearest_rank(array, sorted_idx, percentile, cum_weights)
    return _interpolation_closest_ranks(
        array, sorted_idx, sample_weight, sorted_weights, percentile,
        cum_weights, interpolation,
    )
