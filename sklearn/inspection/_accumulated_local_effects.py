"""Accumulated local effects for regression and classification models."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier, is_regressor
from sklearn.inspection._pd_utils import _check_feature_names
from sklearn.utils import Bunch, _safe_indexing, check_array
from sklearn.utils._indexing import (
    _determine_key_type,
    _get_column_indices,
    _safe_assign,
)
from sklearn.utils._mask import _get_mask
from sklearn.utils._param_validation import (
    HasMethods,
    Integral,
    Interval,
    StrOptions,
    validate_params,
)
from sklearn.utils._response import _get_response_values
from sklearn.utils.validation import _check_categorical_features, check_is_fitted


def _replace_feature_values_and_predict(estimator, X, feature, values, response_method):
    """Replace `feature` in a copy of `X` by `values` and return the predictions.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data where `feature` is replaced by `values`.

    feature : int
        Index of the feature to replace.

    values : ndarray of shape (n_samples,)
        New value of the feature for each sample.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions, passed
        to :func:`~sklearn.utils._response._get_response_values`.

    Returns
    -------
    predictions : ndarray of shape (n_samples, n_outputs)
        The predictions on the modified data. `n_outputs` is 1 for regression
        and binary classification, or greater for multi-class and multi-output
        regression.
    """
    # TODO: decide how to handle libraries other than numpy and pandas
    X_eval = X.copy()
    _safe_assign(X_eval, values, column_indexer=feature)
    pred, _ = _get_response_values(estimator, X_eval, response_method=response_method)
    return pred.reshape(pred.shape[0], -1)


def _mean_local_effect_per_bin(local_effects, bins, n_bins):
    """Average the local effects of the samples that fall in each bin.

    Each row has one local effect: the difference of its predictions when its
    feature value is moved from one end of its bin to the other. A bin is an
    interval of a continuous feature, or the step between two adjacent
    integers/values/categories. A bin contains several samples and therefore
    several local effects. This function returns, for each bin, the mean of those
    local effects: the sum of the local effects in bin `k` divided by the number
    of samples in bin `k`, as in the ALE estimator of Apley and Zhu (2020), Eq. 15.
    The caller then accumulates these means over the bins.

    For bins with no samples, it returns 0.

    Parameters
    ----------
    local_effects : ndarray of shape (n_rows, n_outputs)
        Local effect of each row. One column for regression or binary classif,
        one per class/task for multi-class or multioutput regression.

    bins : ndarray of shape (n_rows,), dtype=int
        Index of the bin of each row, in `[0, n_bins)`.

    n_bins : int
        Number of bins. Bins without rows get a mean of 0.

    Returns
    -------
    means : ndarray of shape (n_outputs, n_bins)
        Mean local effect in each bin, for each output.

    bin_n_samples : ndarray of shape (n_bins,)
        Number of samples in each bin. It is used by the caller to center the
        effects.
    """
    # The number of samples per bin
    bin_n_samples = np.bincount(bins, minlength=n_bins)
    non_empty = bin_n_samples > 0
    means = np.zeros((local_effects.shape[1], n_bins))
    for output in range(local_effects.shape[1]):
        # The sum of the local effects per bin
        sums = np.bincount(bins, weights=local_effects[:, output], minlength=n_bins)
        # The mean local effect per bin
        means[output, non_empty] = sums[non_empty] / bin_n_samples[non_empty]
    return means, bin_n_samples


def _ale_continuous(estimator, X, feature, x, edges, response_method):
    """ALE of a continuous feature, evaluated at the bin `edges`.

    Each sample is assigned to the limits of the interval it falls into:
    `(edges[k], edges[k + 1]]`. The local effect of a sample is the difference
    of the predictions when its feature is set to the upper and to the lower
    bound of its interval. The local effects are averaged within each interval,
    accumulated from the first interval, and centered by subtracting the mean of
    the interval midpoints, weighted by the number of samples, as in ALEPlot
    (Eq. 16 of Apley and Zhu (2020) uses the upper edges instead).

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data.

    feature : int
        Index of the feature of interest in `X`.

    x : ndarray of shape (n_samples,)
        Values of the feature of interest.

    edges : ndarray of shape (n_bins + 1,)
        Sorted bounds of the intervals.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions.

    Returns
    -------
    ale : ndarray of shape (n_outputs, n_bins + 1)
        Centered accumulated local effects at each value of `edges`.

    bin_sizes : ndarray of shape (n_bins,)
        Number of samples in each interval.
    """
    n_bins = edges.shape[0] - 1
    bins = np.clip(np.searchsorted(edges, x, side="left"), 1, n_bins) - 1

    pred_upper = _replace_feature_values_and_predict(
        estimator, X, feature, edges[bins + 1], response_method
    )
    pred_lower = _replace_feature_values_and_predict(
        estimator, X, feature, edges[bins], response_method
    )

    mean_local_effects, bin_sizes = _mean_local_effect_per_bin(
        pred_upper - pred_lower, bins, n_bins
    )
    ale = np.hstack(
        [
            np.zeros((mean_local_effects.shape[0], 1)),
            np.cumsum(mean_local_effects, axis=1),
        ]
    )
    midpoints = (ale[:, :-1] + ale[:, 1:]) / 2
    ale -= np.average(midpoints, axis=1, weights=bin_sizes)[:, np.newaxis]
    return ale, bin_sizes


def _ale_categorical(estimator, X, feature, x, levels, response_method):
    """ALE of a discrete or categorical feature with ordered `levels`.

    The changes between neighbouring values are called bridges. For example, a
    feature with values 1, 3, 5 and 7 has three bridges: 1 -> 3, 3 -> 5 and
    5 -> 7. Each bridge is measured with the samples at its two ends::

        bridge 1 -> 3:  samples at 1 moved up to 3,  samples at 3 moved down to 1
        bridge 3 -> 5:  samples at 3 moved up to 5,  samples at 5 moved down to 3
        bridge 5 -> 7:  samples at 5 moved up to 7,  samples at 7 moved down to 5

    The local effect of each move is the prediction at the upper end of the bridge
    minus the prediction at the lower end, keeping the other features unchanged.
    The effect of a bridge is the mean of the local effects of its moves. Samples
    at 3 and 5 are used in two bridges, once in each.

    The bridges are then accumulated from the first value::

        ALE at 1 = 0
        ALE at 3 = bridge 1 -> 3
        ALE at 5 = bridge 1 -> 3 + bridge 3 -> 5
        ALE at 7 = bridge 1 -> 3 + bridge 3 -> 5 + bridge 5 -> 7

    and centered, so that the mean of the ALE over the samples is zero.

    In the continuous case, each sample lies inside one interval and is moved to
    its two edges, so each local effect covers exactly one interval. In the
    categorical case, the bridges play the role of the intervals, and a sample
    sits on the edge between two bridges, so it gives one local effect for each.
    Every local effect still covers exactly one bridge.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data.

    feature : int
        Index of the feature of interest in `X`.

    x : ndarray of shape (n_samples,)
        Values of the feature of interest.

    levels : ndarray of shape (n_levels,)
        Unique values of the feature.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions.

    Returns
    -------
    ale : ndarray of shape (n_outputs, n_levels)
        Centered accumulated local effects at each value of `levels`.

    level_sizes : ndarray of shape (n_levels,)
        Number of samples at each level.

    References
    ----------
    Local effects of categorical features in the ALEPlot R package by Apley:
    https://github.com/cran/ALEPlot/blob/d283def19a5d7f1840f18dfca82c073210df0d81/R/ALEPlot.R#L53-L67
    """
    n_levels = levels.shape[0]
    sorter = np.argsort(levels)
    level_idx = sorter[np.searchsorted(levels, x, sorter=sorter)]

    pred = _replace_feature_values_and_predict(
        estimator, X, feature, levels[level_idx], response_method
    )
    pred_up = _replace_feature_values_and_predict(
        estimator,
        X,
        feature,
        levels[np.minimum(level_idx + 1, n_levels - 1)],
        response_method,
    )
    pred_down = _replace_feature_values_and_predict(
        estimator, X, feature, levels[np.maximum(level_idx - 1, 0)], response_method
    )

    # Samples that can move up or down (i.e., not the edges of the
    # unique value range).
    up = level_idx < n_levels - 1
    down = level_idx > 0
    # Step k goes from levels[k] to levels[k + 1]: moving up from k, down from k + 1.
    steps = np.concatenate([level_idx[up], level_idx[down] - 1])
    # Average the prediction differences of both groups within each step.
    mean_local_effects, _ = _mean_local_effect_per_bin(
        np.concatenate([(pred_up - pred)[up], (pred - pred_down)[down]]),
        steps,
        n_levels - 1,
    )
    # Accumulate the step effects from the first level, which starts at 0.
    ale = np.hstack(
        [
            np.zeros((mean_local_effects.shape[0], 1)),
            np.cumsum(mean_local_effects, axis=1),
        ]
    )
    # Number of samples at each level.
    level_sizes = np.bincount(level_idx, minlength=n_levels)
    # Center so that the average effect over the samples is zero.
    ale -= np.average(ale, axis=1, weights=level_sizes)[:, np.newaxis]
    return ale, level_sizes


def _ale_continuous_2d(estimator, X, features, x, edges, response_method):
    """Second-order ALE of two continuous features, evaluated at the grid `edges`.

    The two features are split into intervals, which form a grid of cells. The
    local effect of a sample is its four-corner difference in the cell it falls
    into: f(upper, upper) - f(upper, lower) - f(lower, upper) + f(lower, lower),
    keeping the other features unchanged. The local effects are averaged within
    each cell and accumulated over both features. Then, the first-order effects
    of each feature are removed and the result is centered, as in ALEPlot.

    Cells without samples get a mean local effect of 0. Their ALE is not
    imputed from neighbouring cells, and they can be identified with the
    returned `cell_sizes`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data.

    features : tuple of int
        Indices of the two features of interest in `X`.

    x : tuple of two ndarrays of shape (n_samples,)
        Values of the two features of interest.

    edges : tuple of two ndarrays of shape (n_bins_0 + 1,) and (n_bins_1 + 1,)
        Sorted bounds of the intervals of each feature.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions.

    Returns
    -------
    ale : ndarray of shape (n_outputs, n_bins_0 + 1, n_bins_1 + 1)
        Centered second-order accumulated local effects at each grid point.

    cell_sizes : ndarray of shape (n_bins_0, n_bins_1)
        Number of samples in each cell.

    References
    ----------
    Second-order effects of two numeric features in the ALEPlot R package by Apley:
    https://github.com/cran/ALEPlot/blob/d283def19a5d7f1840f18dfca82c073210df0d81/R/ALEPlot.R#L240-L309
    """
    n_bins = [e.shape[0] - 1 for e in edges]
    bins = [
        np.clip(np.searchsorted(e, xi, side="left"), 1, n) - 1
        for e, xi, n in zip(edges, x, n_bins)
    ]

    # Predictions at the four corners of the cell of each sample.
    pred = {}
    for i in (0, 1):
        X_first = X.copy()
        _safe_assign(X_first, edges[0][bins[0] + i], column_indexer=features[0])
        for j in (0, 1):
            pred[i, j] = _replace_feature_values_and_predict(
                estimator, X_first, features[1], edges[1][bins[1] + j], response_method
            )
    local_effects = (pred[1, 1] - pred[1, 0]) - (pred[0, 1] - pred[0, 0])

    # Mean local effect per cell, with the cells numbered row by row.
    mean_local_effects, cell_sizes = _mean_local_effect_per_bin(
        local_effects, bins[0] * n_bins[1] + bins[1], n_bins[0] * n_bins[1]
    )
    mean_local_effects = mean_local_effects.reshape(-1, *n_bins)
    cell_sizes = cell_sizes.reshape(n_bins)

    # Accumulate over both features, starting at 0 at the first grid point.
    ale = np.zeros((mean_local_effects.shape[0], n_bins[0] + 1, n_bins[1] + 1))
    ale[:, 1:, 1:] = np.cumsum(np.cumsum(mean_local_effects, axis=1), axis=2)

    # Remove the first-order effect of the first feature.
    diff_0 = ale[:, 1:, :] - ale[:, :-1, :]
    weighted = cell_sizes * (diff_0[:, :, :-1] + diff_0[:, :, 1:]) / 2
    effect_0 = weighted.sum(axis=2) / cell_sizes.sum(axis=1)
    ale_0 = np.hstack([np.zeros((ale.shape[0], 1)), np.cumsum(effect_0, axis=1)])

    # Remove the first-order effect of the second feature.
    diff_1 = ale[:, :, 1:] - ale[:, :, :-1]
    weighted = cell_sizes * (diff_1[:, :-1, :] + diff_1[:, 1:, :]) / 2
    effect_1 = weighted.sum(axis=1) / cell_sizes.sum(axis=0)
    ale_1 = np.hstack([np.zeros((ale.shape[0], 1)), np.cumsum(effect_1, axis=1)])

    ale -= ale_0[:, :, np.newaxis] + ale_1[:, np.newaxis, :]

    # Center with the mean of the four corners of each cell, weighted by the
    # number of samples in the cell.
    corners = (
        ale[:, :-1, :-1] + ale[:, :-1, 1:] + ale[:, 1:, :-1] + ale[:, 1:, 1:]
    ) / 4
    ale -= (np.sum(cell_sizes * corners, axis=(1, 2)) / cell_sizes.sum())[
        :, np.newaxis, np.newaxis
    ]
    return ale, cell_sizes


@validate_params(
    {
        "estimator": [
            HasMethods(["fit", "predict"]),
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
        ],
        "X": ["array-like"],
        "features": ["array-like", Integral, str],
        "categorical_features": ["array-like", StrOptions({"from_dtype"}), None],
        "feature_names": ["array-like", None],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function"})],
        "grid_resolution": [Interval(Integral, 2, None, closed="left")],
    },
    prefer_skip_nested_validation=True,
)
def accumulated_local_effects(
    estimator,
    X,
    features,
    *,
    categorical_features="from_dtype",
    feature_names=None,
    response_method="auto",
    grid_resolution=100,
):
    """Accumulated local effects (ALE) of ``features``.

    The range of the feature is split into intervals. Within each interval, the
    feature of the samples that fall into it is set to the lower and upper
    interval bounds, and the difference in the predictions is averaged. These
    local effects are accumulated over the intervals and centered so that their
    mean over the samples is zero [1]_.

    Unlike partial dependence, ALE only evaluates the estimator in the vicinity
    of observed samples and therefore prevents or minimizes the occurrence of
    feature-value combinations that would not occur in the real-world.

    Read more in the :ref:`User Guide <accumulated_local_effects>`.

    .. versionadded:: 1.10

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.
        Multioutput-multiclass classifiers are not supported.

    X : {array-like, dataframe} of shape (n_samples, n_features)
        Data used to define the intervals and to compute the local effects.
        Samples with a missing value in the target feature are ignored.

    features : int, str or array-like of {int, str, bool}
        The feature (e.g. `0`) or pair of interacting features (e.g. `[0, 1]`)
        for which the accumulated local effects are computed. Pairs are only
        supported for continuous features.

    categorical_features : array-like of {bool, int, str} of shape (n_features,) or \
            shape (n_categorical_features,), "from_dtype" or None, \
            default="from_dtype"
        Indicates the categorical or discrete features.

        - `None`: no feature is considered categorical.
        - boolean array-like: boolean mask of shape `(n_features,)` indicating
          which features are categorical.
        - integer array-like: integer indices indicating categorical features.
        - str array-like: names of categorical features (`X` must be a
          dataframe).
        - `"from_dtype"`: dataframe columns with dtype "Categorical" and "Enum"
          are considered categorical. For other inputs, no feature is
          categorical.

        For these features, the local effect between two consecutive observed
        values is computed by moving the samples at the lower value to the
        upper value and the samples at the upper value to the lower value.

        .. warning::

            Categories are ordered as returned by :func:`numpy.unique`. For
            nominal features, this order is arbitrary and the accumulated local
            effects depend on it.

    feature_names : array-like of shape (n_features,), dtype=str, default=None
        Name of each feature; `feature_names[i]` holds the name of the feature
        with index `i`.
        By default, the name of the feature corresponds to their numerical
        index for NumPy array and their column name for dataframes.

    response_method : {'auto', 'predict_proba', 'decision_function'}, \
            default='auto'
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist.

    grid_resolution : int, default=100
        Maximum number of points at which the accumulated local effects of a
        continuous feature are computed. The range of the feature is split into
        `grid_resolution - 1` intervals with roughly the same number of samples.

    Returns
    -------
    result : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        effects : ndarray of shape (n_outputs, n_grid_values) or \
                (n_outputs, n_grid_values_0, n_grid_values_1)
            The accumulated local effects at each value in `grid_values`. For a
            pair of features, the second-order effects at each grid point.

        grid_values : list of 1d ndarrays
            One array per feature in `features`: the interval bounds for
            continuous features, or the categories for categorical features.

        bin_sizes : ndarray
            The number of samples in each interval, of shape
            `(n_grid_values - 1,)`, or in each category, of shape
            `(n_grid_values,)`. For a pair of features, the number of samples in
            each cell of the grid, of shape
            `(n_grid_values_0 - 1, n_grid_values_1 - 1)`. Cells without samples
            have no data to estimate their effects.

        `n_outputs` corresponds to the number of classes in a multi-class
        setting, or to the number of tasks for multi-output regression.
        For classical regression and binary classification `n_outputs==1`.

    See Also
    --------
    AccumulatedLocalEffectsDisplay.from_estimator : Plot accumulated local effects.
    partial_dependence : Compute partial dependence.

    References
    ----------
    .. [1] D. W. Apley and J. Zhu, "Visualizing the effects of predictor
       variables in black box supervised learning models", Journal of the
       Royal Statistical Society Series B, 82(4), 1059-1086, 2020.
       :doi:`10.1111/rssb.12377`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.inspection import accumulated_local_effects
    >>> from sklearn.linear_model import LinearRegression
    >>> X = np.array([[0.0, 1.0], [1.0, 0.0], [2.0, 1.0], [3.0, 0.0]])
    >>> y = 2 * X[:, 0]
    >>> reg = LinearRegression().fit(X, y)
    >>> result = accumulated_local_effects(reg, X, features=0, grid_resolution=4)
    >>> result.grid_values
    [array([0., 1., 2., 3.])]
    >>> result.effects
    array([[-2.5, -0.5,  1.5,  3.5]])
    """
    check_is_fitted(estimator)

    if not (is_classifier(estimator) or is_regressor(estimator)):
        raise ValueError("'estimator' must be a fitted regressor or classifier.")

    if is_classifier(estimator) and isinstance(estimator.classes_[0], np.ndarray):
        raise ValueError("Multiclass-multioutput estimators are not supported")

    # Use check_array only on lists and other non-array-likes / sparse. Do not
    # convert DataFrame into a NumPy array.
    if not hasattr(X, "__array__"):
        X = check_array(X, ensure_all_finite="allow-nan", dtype=object)

    if is_regressor(estimator) and response_method != "auto":
        raise ValueError(
            "The response_method parameter is ignored for regressors and "
            "must be 'auto'."
        )
    if response_method == "auto":
        response_method = (
            "predict"
            if is_regressor(estimator)
            else ["predict_proba", "decision_function"]
        )

    if _determine_key_type(features, accept_slice=False) == "int":
        # _get_column_indices() supports negative indexing. Here, we limit
        # the indexing to be positive. The upper bound will be checked
        # by _get_column_indices()
        if np.any(np.less(features, 0)):
            raise ValueError(f"all features must be in [0, {X.shape[1] - 1}]")

    features_indices = np.asarray(
        _get_column_indices(X, features), dtype=np.intp, order="C"
    ).ravel()

    if features_indices.size > 2:
        raise ValueError(
            "`features` must be a single feature or a pair of features. Got "
            f"{features_indices.size} features."
        )

    feature_names = _check_feature_names(X, feature_names)

    categorical_mask = _check_categorical_features(X, categorical_features)
    if categorical_mask is None:
        is_categorical = [False] * features_indices.size
    else:
        is_categorical = [bool(categorical_mask[idx]) for idx in features_indices]

    X_subset = _safe_indexing(X, features_indices, axis=1)

    # A list with one 1d array per feature of interest, holding its values for all
    # samples.
    observed_values = [
        np.asarray(_safe_indexing(X_subset, feature, axis=1))
        for feature in range(len(is_categorical))
    ]

    # Ignore the samples with a missing value in any of the features of interest.
    valid = ~np.any([_get_mask(x, np.nan) for x in observed_values], axis=0)
    if not np.any(valid):
        raise ValueError(
            "No sample is left to compute the accumulated local effects of "
            f"{[feature_names[idx] for idx in features_indices]!r}: all samples "
            "have a missing value."
        )

    # Build the grid of each feature from the samples that are left.
    grid_values = []
    for x, is_cat in zip(observed_values, is_categorical):
        if is_cat:
            # Categorical: the categories observed in the samples that are left.
            axis = np.unique(x[valid])
        else:
            # Continuous: quantiles.
            axis = np.percentile(
                x[valid].astype(np.float64),
                np.linspace(0, 100, grid_resolution),
                method="inverted_cdf",
            )
            # Merge repeated quantiles and keep the dtype of the feature.
            axis = np.unique(axis).astype(x.dtype)
        grid_values.append(axis)

    for feature_idx, axis in zip(features_indices, grid_values):
        if axis.shape[0] < 2:
            raise ValueError(
                f"Feature {feature_names[feature_idx]!r} must have at least two "
                "distinct values to compute accumulated local effects. Got "
                f"{axis!r}."
            )

    if features_indices.size == 1:
        column = features_indices[0]
        # Keep only the samples that are not ignored.
        X_valid = _safe_indexing(X, np.flatnonzero(valid), axis=0)
        x = observed_values[0][valid]
        if is_categorical[0]:
            effects, bin_sizes = _ale_categorical(
                estimator, X_valid, column, x, grid_values[0], response_method
            )
        else:
            effects, bin_sizes = _ale_continuous(
                estimator, X_valid, column, x, grid_values[0], response_method
            )
    else:
        if any(is_categorical):
            # TODO: add support for pairs with categorical features (ALEPlot
            # supports a categorical first feature with a continuous second one).
            raise NotImplementedError(
                "Second-order accumulated local effects are only supported for "
                "two continuous features."
            )
        X_valid = _safe_indexing(X, np.flatnonzero(valid), axis=0)
        effects, bin_sizes = _ale_continuous_2d(
            estimator,
            X_valid,
            tuple(features_indices),
            tuple(x[valid] for x in observed_values),
            tuple(grid_values),
            response_method,
        )

    return Bunch(effects=effects, grid_values=grid_values, bin_sizes=bin_sizes)
