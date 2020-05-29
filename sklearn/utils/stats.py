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
    array : 1D or 2D array
        Values to take the weighted percentile of.

    sample_weight: 1D or 2D array
        Weights for each value in `array`. Must be same shape as `array` or
        of shape `(array.shape[0],)`.

    percentile: int, default=50
        Percentile to compute. Must be value between 0 and 100.

    interpolation : {"linear", "lower", "higher"}, default="linear"
        FIXME

    Returns
    -------
    percentile_value : int if `array` 1D, ndarray if `array` 2D
        Weighted percentile.
    """
    possible_interpolation = ("linear", "lower", "higher")
    if interpolation not in possible_interpolation:
        raise ValueError(
            f"'interpolation' should be one of "
            f"{', '.join(possible_interpolation)}. Got '{interpolation}' "
            f"instead."
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
    weight_cdf = stable_cumsum(sorted_weights, axis=0)

    def _squeeze_arr(arr, n_dim):
        return arr[0] if n_dim == 1 else arr

    adjusted_percentile = (weight_cdf - sorted_weights)
    with np.errstate(invalid="ignore"):
        adjusted_percentile /= ((weight_cdf[-1] * (n_rows - 1)) / (n_rows))

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
        return _squeeze_arr(percentile_value, n_dim)

    else:  # interpolation == "linear"
        percentile_value = np.array([
            np.interp(
                x=percentile[col],
                xp=adjusted_percentile[:, col],
                fp=array[sorted_idx[:, col], col],
            )
            for col in range(n_cols)
        ])
        return _squeeze_arr(percentile_value, n_dim)
