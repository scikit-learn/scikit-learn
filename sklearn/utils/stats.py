from collections.abc import Iterable

import numpy as np

from .extmath import stable_cumsum
from .fixes import _take_along_axis


def _squeeze_arr(arr, n_dim):
    return arr[0] if n_dim == 1 else arr


def _nearest_rank(
    array, sorted_idx, percentile, cum_weights, n_rows, n_cols, n_dim,
):
    """Compute weighted percentile using the nearest-rank approach.

    Percentile can be computed using the nearest-rank:
    https://en.wikipedia.org/wiki/Percentile

    n = ceil(P / 100 * N)

    where n is the index of the percentile value of the percentile P and N is
    the number of samples.

    This definition can be adapted with weights:

    n_w = ceil(P / 100 * S_w) where S_w is the sum of the weights.
    """
    adjusted_percentile = percentile * cum_weights[-1]
    percentile_idx = np.array([
        np.searchsorted(cum_weights[:, i], adjusted_percentile[i])
        for i in range(cum_weights.shape[1])
    ])
    percentile_idx = np.array(percentile_idx)
    # In rare cases, percentile_idx equals to sorted_idx.shape[0]
    max_idx = n_rows - 1
    percentile_idx = np.apply_along_axis(
        lambda x: np.clip(x, 0, max_idx), axis=0, arr=percentile_idx
    )

    col_index = np.arange(n_cols)
    percentile_in_sorted = sorted_idx[percentile_idx, col_index]
    percentile_value = array[percentile_in_sorted, col_index]
    return _squeeze_arr(percentile_value, n_dim)


def _interpolation_closest_ranks(
    array, sorted_idx, sample_weight, sorted_weights, percentile, cum_weights,
    interpolation, n_rows, n_cols, n_dim,
):
    """Compute the weighted percentile by interpolating between the adjacent
    ranks.

    Percentile can be computed with 3 different alternatives when using the
    interpolation between adjacent ranks:
    https://en.wikipedia.org/wiki/Percentile

    These 3 alternatives depend of the value of a parameter C. NumPy uses
    the variant where C=1 which allows to obtain a strictly monotonically
    increasing function which is defined as:

    P = (x - 1) / (N - 1),

    where x in [1, N]

    Weighted percentile change this formula by taking into account the
    weights instead of the data frequency.

    P_w_i = (x - w_i) / (S_w - w_i),

    where x in [1, N], w_i being the weight of sample i and S_w being the sum
    of the weights.
    """
    def _fill_nan(array, value):
        nan_mask = np.isnan(array)
        array[nan_mask] = value
        return array

    adjusted_percentile = (cum_weights - sorted_weights)
    with np.errstate(invalid="ignore"):
        adjusted_percentile /= cum_weights[-1] - sorted_weights
        adjusted_percentile = _fill_nan(adjusted_percentile, 1)

    percentile_higher_idx = np.array([
        np.searchsorted(
            adjusted_percentile[:, col], percentile[col], side="left",
        )
        for col in range(n_cols)
    ])

    if interpolation != "higher":
        # other interpolations which requires the lower percentile indice
        percentile_lower_idx = np.array([
            idx - 1 if val < 1 else idx
            for idx, val in zip(percentile_higher_idx, percentile)
        ])

        if interpolation == "nearest":

            def _compute_abs_error(p_measured, p_desired, idx, col):
                return abs(p_measured[idx[col], col] - p_desired[col])

            error_higher = np.array([
                _compute_abs_error(
                    p_measured=adjusted_percentile,
                    p_desired=percentile,
                    idx=percentile_higher_idx,
                    col=col,
                )
                for col in range(n_cols)
            ])
            error_lower = np.array([
                _compute_abs_error(
                    p_measured=adjusted_percentile,
                    p_desired=percentile,
                    idx=percentile_lower_idx,
                    col=col,
                )
                for col in range(n_cols)
            ])

            percentile_idx = np.where(
                error_higher >= error_lower,
                percentile_lower_idx,
                percentile_higher_idx
            )
        elif interpolation == "lower":
            percentile_idx = percentile_lower_idx
    else:
        percentile_idx = percentile_higher_idx

    if interpolation in ("weighted", "linear"):
        percentile_value_lower = array[
            sorted_idx[percentile_lower_idx, np.arange(n_cols)],
            np.arange(n_cols)
        ]
        percentile_value_higher = array[
            sorted_idx[percentile_higher_idx, np.arange(n_cols)],
            np.arange(n_cols)
        ]

        if interpolation == "weighted":
            weight_higher = sorted_weights[
                percentile_higher_idx, np.arange(n_cols)
            ]
            weight_lower = sorted_weights[
                percentile_lower_idx, np.arange(n_cols)
            ]
            weight_normalizer = weight_higher + weight_lower
            with np.errstate(invalid="ignore"):
                weight_higher = weight_higher / weight_normalizer
                weight_higher = _fill_nan(weight_higher, 0)
                weight_lower = weight_lower / weight_normalizer
                weight_lower = _fill_nan(weight_lower, 0)

        else:
            # linear interpolation does not use weights
            weight_lower = np.zeros(shape=(n_cols,))
            weight_higher = np.zeros(shape=(n_cols,))

        def _weighted_interpolation(
            x, x_low, x_high, y_low, y_high, w_low, w_high,
        ):
            """Weighted interpolation between 2 values.

            If the weights are null, then it is equivalent to a linear
            interpolation.
            """
            with np.errstate(invalid="ignore"):
                percentile_value = (
                    y_low + ((1 - w_low) * (x - x_low)) * (
                        (y_high - y_low) / (
                            (1 - w_high) * (x_high - x) +
                            (1 - w_low) * (x - x_low)
                        )
                    )
                )
            if np.isnan(percentile_value):
                # case that percentile is 100
                return y_high
            return percentile_value

        percentile_value = np.array([
            _weighted_interpolation(
                x=percentile[col],
                x_low=adjusted_percentile[percentile_lower_idx[col], col],
                x_high=adjusted_percentile[percentile_higher_idx[col], col],
                y_low=percentile_value_lower[col],
                y_high=percentile_value_higher[col],
                w_low=weight_lower[col],
                w_high=weight_higher[col],
            )
            for col in range(n_cols)
        ])
    else:
        # no need to interpolate
        percentile_idx = np.apply_along_axis(
            lambda x: np.clip(x, 0, n_rows - 1), axis=0,
            arr=percentile_idx
        )

        percentile_value = array[
            sorted_idx[percentile_idx, np.arange(n_cols)],
            np.arange(n_cols)
        ]

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


def _weighted_percentile(
    array, sample_weight, percentile=50, interpolation=None,
):
    """Compute weighted percentile.

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

    percentile: float, default=50
        Percentile to compute. Must be value between 0 and 100.

    interpolation : {"linear", "lower", "higher", "nearest"}, default=None
        The interpolation method to use when the percentile lies between
        data points `i` and `j`:

        * None: no interpolation will be done and the "nearest-rank" method
          will be used (default).
        * "linear": linearly interpolate between `i` and `j`;
        * "weighted": interpolate between `i` and `j` taking into account the
          weights assigned to `i` and `j`.
        * "lower": i`;
        * "higher": `j`;
        * "nearest": `i` or `j`, whichever is nearest.

        .. versionadded: 0.24

    Returns
    -------
    percentile_value : float or int if `array` of shape (n,), otherwise \
            ndarray of shape (m,)
        Weighted percentile.
    """
    possible_interpolation = (
        "linear", "lower", "higher", "nearest", "weighted", None
    )
    if interpolation not in possible_interpolation:
        raise ValueError(
            f"'interpolation' should be one of "
            f"{', '.join([str(val) for val in possible_interpolation])}. "
            f"Got '{interpolation}' instead."
        )

    if np.any(np.count_nonzero(sample_weight, axis=0) < 1):
        raise ValueError(
            "All weights cannot be null when computing a weighted percentile."
        )

    n_dim = array.ndim
    if n_dim == 0:
        return array[()]
    elif n_dim == 1:
        array = array.reshape((-1, 1))

    if (array.shape != sample_weight.shape and
            array.shape[0] == sample_weight.shape[0]):
        # when `sample_weight` is 1D, we repeat it for each column of `array`
        sample_weight = np.tile(sample_weight, (array.shape[1], 1)).T

    n_rows, n_cols = array.shape
    sorted_idx = np.argsort(array, axis=0)
    sorted_weights = _take_along_axis(sample_weight, sorted_idx, axis=0)
    percentile = np.array([percentile / 100] * n_cols)
    cum_weights = stable_cumsum(sorted_weights, axis=0)

    if interpolation is None:
        return _nearest_rank(
            array, sorted_idx, percentile, cum_weights, n_rows, n_cols, n_dim,
        )
    return _interpolation_closest_ranks(
        array, sorted_idx, sample_weight, sorted_weights, percentile,
        cum_weights, interpolation, n_rows, n_cols, n_dim,
    )
