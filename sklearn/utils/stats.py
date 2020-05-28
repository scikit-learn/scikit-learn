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

    Returns
    -------
    percentile_value : int if `array` 1D, ndarray if `array` 2D
        Weighted percentile.
    """
    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    if array.ndim == 1:
        array = array.reshape((-1, 1))

    if (array.shape != sample_weight.shape and
            array.shape[0] == sample_weight.shape[0]):
        # when `sample_weight` is 1D, we repeat it for each column of `array`
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)

    # find the lower percentile value indices
    weight_cdf = stable_cumsum(sorted_weights, axis=0)
    adjusted_percentile = percentile / 100 * weight_cdf[-1]
    percentile_value_lower_idx = np.array([
        np.searchsorted(weight_cdf[:, i], adjusted_percentile[i])
        for i in range(weight_cdf.shape[1])
    ])

    # clip the indices in case the indices is the last element found
    max_idx = sorted_idx.shape[0] - 1
    col_index = np.arange(array.shape[1])

    def _reduce_dim(arr, n_dim):
        return arr[0] if n_dim == 1 else arr

    if interpolation in ("lower", "linear"):
        percentile_value_lower_idx = np.apply_along_axis(
            lambda x: np.clip(x, 0, max_idx), axis=0,
            arr=percentile_value_lower_idx,
        )
        percentile_value_lower = array[
            sorted_idx[percentile_value_lower_idx, col_index]
        ]
        if interpolation == "lower":
            return _reduce_dim(percentile_value_lower, n_dim)

    if interpolation in ("higher", "linear"):
        percentile_value_higher_idx = np.apply_along_axis(
            lambda x: np.clip(x, 0, max_idx), axis=0,
            arr=percentile_value_lower_idx + 1,
        )
        percentile_value_higher = array[
            sorted_idx[percentile_value_higher_idx, col_index]
        ]
        if interpolation == "higher":
            return _reduce_dim(percentile_value_higher, n_dim)

    ratio = percentile / adjusted_percentile
    percentile_lower = weight_cdf[percentile_value_lower_idx] * ratio
    percentile_higher = weight_cdf[percentile_value_higher_idx] * ratio

    # interpolate linearly for the given percentile
    percentile_value = (
        percentile_value_lower + (percentile - percentile_lower) *
        ((percentile_value_higher - percentile_value_lower) /
         (percentile_higher - percentile_lower))
    )
    print(percentile_higher, percentile_lower)
    print(percentile_value_higher, percentile_value_lower)
    print(percentile_value)
    print(weight_cdf)

    return _reduce_dim(percentile_value, n_dim)
