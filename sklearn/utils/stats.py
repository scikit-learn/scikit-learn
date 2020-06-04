from collections.abc import Iterable

import numpy as np

from .extmath import stable_cumsum
from .fixes import _take_along_axis


def _weighted_percentile(array, sample_weight, percentile=50,
                         interpolation="linear"):
    """Compute weighted percentile

    Computes lower weighted percentile. If `array` is a 2D array, the
    `percentile` is computed along the axis 0.

        .. versionchanged:: 0.24
            Accepts 2D `array`.

    Parameters
    ----------
    array : ndarray of shape (n,) or (n, m)
        Values to take the weighted percentile of.

    sample_weight: ndarray of (n,) or (n, m)
        Weights for each value in `array`. Must be same shape as `array` or
        of shape `(array.shape[0],)`.

    percentile: inr or float, default=50
        Percentile to compute. Must be value between 0 and 100.

    interpolation : {"linear", "lower", "higher"}, default="linear"
        The interpolation method to use when the percentile lies between
        data points `i` and `j`:

        * `"linear"`: `i + (j - i) * fraction`, where `fraction` is the
          fractional part of the index surrounded by `i` and `j`;
        * `"lower"`: i`;
        * `"higher"`: `j`.

        .. versionadded: 0.24

    Returns
    -------
    percentile_value : float or int if `array` of shape (n,), otherwise\
            ndarray of shape (m,)
        Weighted percentile.
    """
    possible_interpolation = ("linear", "lower", "higher")
    if interpolation not in possible_interpolation:
        raise ValueError(
            f"'interpolation' should be one of "
            f"{', '.join(possible_interpolation)}. Got '{interpolation}' "
            f"instead."
        )

    if np.any(np.count_nonzero(sample_weight, axis=0) < 1):
        raise ValueError(
            "All weights cannot be null when computing a weighted percentile."
        )

    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))

    if (array.shape != sample_weight.shape and
            array.shape[0] == sample_weight.shape[0]):
        # when `sample_weight` is 1D, we repeat it for each column of `array`
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T

    n_rows, n_cols = array.shape

    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)
    percentile = np.array([percentile / 100] * n_cols)
    cum_weigths = stable_cumsum(sorted_weights, axis=0)

    def _squeeze_arr(arr, n_dim):
        return arr[0] if n_dim == 1 else arr

    # Percentile can be computed with 3 different alternative:
    # https://en.wikipedia.org/wiki/Percentile
    # These 3 alternatives depend of the value of a parameter C. NumPy uses
    # the variant where C=0 which allows to obtained a strictly monotically
    # increasing function which is defined as:
    # P = (x - 1) / (N - 1); x in [1, N]
    # Weighted percentile change this formula by taking into account the
    # weights instead of the data frequency.
    # P_w = (x - w) / (S_w - w), x in [1, N], w being the weight and S_n being
    # the sum of the weights.
    adjusted_percentile = (cum_weigths - sorted_weights)
    with np.errstate(invalid="ignore"):
        adjusted_percentile /= cum_weigths[-1] - sorted_weights
        nan_mask = np.isnan(adjusted_percentile)
        adjusted_percentile[nan_mask] = 1

    if interpolation in ("lower", "higher"):
        percentile_idx = np.array([
            np.searchsorted(adjusted_percentile[:, col], percentile[col],
                            side="left")
            for col in range(n_cols)
        ])

        if interpolation == "lower" and np.all(percentile < 1):
            # P = 100 is a corner case for "lower"
            percentile_idx -= 1

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
