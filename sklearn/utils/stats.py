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
        This optional parameter specifies the interpolation method to
        use when the desired percentile lies between two data points
        ``i < j``:
        * 'linear': ``i + (j - i) * fraction``, where ``fraction``
          is the fractional part of the index surrounded by ``i``
          and ``j``.
        * 'lower': ``i``.
        * 'higher': ``j``.
        * 'nearest': ``i`` or ``j``, whichever is nearest.

    Returns
    -------
    percentile_value : int if `array` 1D, ndarray if `array` 2D
        Weighted percentile.
    """
    possible_interpolation = ("linear", "lower", "higher", "nearest")
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
    percentile = [percentile / 100] * n_cols
    weight_cdf = stable_cumsum(sorted_weights, axis=0)

    def _squeeze_arr(arr, n_dim):
        return arr[0] if n_dim == 1 else arr

    if interpolation == "nearest":
        # compute by nearest-rank method
        adjusted_percentile = percentile * weight_cdf[-1]
        percentile_value_idx = np.array([
            np.searchsorted(weight_cdf[:, i], adjusted_percentile[i])
            for i in range(weight_cdf.shape[1])
        ])

        percentile_value_idx = np.apply_along_axis(
            lambda x: np.clip(x, 0, n_rows - 1), axis=0,
            arr=percentile_value_idx
        )
        percentile_value = array[
            sorted_idx[percentile_value_idx, np.arange(n_cols)],
            np.arange(n_cols)
        ]
        return _squeeze_arr(percentile_value, n_dim)

    elif interpolation in ("linear", "lower", "higher"):
        # compute by linear interpolation between closest ranks method
        # adjusted_percentile = (weight_cdf - 0.5 * sorted_weights)
        # with np.errstate(invalid="ignore"):
        #     adjusted_percentile /= weight_cdf[-1]
        adjusted_percentile = (weight_cdf - sorted_weights)
        adjusted_percentile /= ((weight_cdf[-1] * (n_rows - 1)) / (n_rows))

        if interpolation in ("lower", "higher"):
            percentile_idx = np.array([
                np.searchsorted(adjusted_percentile[:, col], percentile[col],
                                side="left")
                for col in range(adjusted_percentile.shape[1])
            ])
            if interpolation == "lower":
                percentile_idx -= 1
            percentile_idx = np.apply_along_axis(
                lambda x: np.clip(x, 0, n_rows - 1), axis=0,
                arr=percentile_idx
            )
            percentile = np.nan_to_num([
                adjusted_percentile[percentile_idx[i], i]
                for i in range(adjusted_percentile.shape[1])
            ], nan=percentile)

        percentile_value = np.array([
            np.interp(
                x=percentile[col],
                xp=adjusted_percentile[:, col],
                fp=array[sorted_idx[:, col], col],
            )
            for col in range(array.shape[1])
        ])
        return _squeeze_arr(percentile_value, n_dim)
