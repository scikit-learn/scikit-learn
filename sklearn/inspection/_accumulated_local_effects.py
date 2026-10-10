"""Accumulated local effects for regression and classification models."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier, is_regressor
from sklearn.inspection._pd_utils import _check_feature_names, _get_is_categorical
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
from sklearn.utils.stats import _weighted_percentile
from sklearn.utils.validation import _check_sample_weight, check_is_fitted


def _replace_values_and_predict(estimator, X, column, values, response_method):
    """Predict on a copy of `X` where the feature `column` is replaced by `values`.

    ALE compares the predictions of each sample when only the feature of interest
    is moved, e.g. to the lower and upper bounds of its interval. The other
    features keep their observed values.

    This function replaces the value of the feature `column` in `X` by `values` and
    returns the predictions of `estimator` on this modified data.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data where feature `column` is replaced by `values`.

    column : int
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
        or greater for classification and multi-output regression.
    """
    # TODO: decide how to handle libraries other than numpy and pandas
    X_eval = X.copy()
    _safe_assign(X_eval, values, column_indexer=column)
    pred, _ = _get_response_values(estimator, X_eval, response_method=response_method)
    return pred.reshape(pred.shape[0], -1)


def _mean_local_effect_per_bin(local_effects, bins, sample_weight, n_bins):
    """Average the local effects of the samples that fall in each bin.

    Each sample has one local effect: the difference of its predictions when its
    feature value is moved from one end of its bin to the other. A bin is an interval
    of a continuous feature, or the step between two adjacent integers/categories.
    A bin contains several samples and therefore several local effects. This function
    returns, for each bin, the (weighted) mean of those local effects: the sum of
    the local effects in bin `k` divided by the number of samples in bin `k`, as
    in the ALE estimator of Apley and Zhu (2020), Eq. 15. The caller then
    accumulates these means over the bins.

    For bins with no samples, it returns 0.

    Parameters
    ----------
    local_effects : ndarray of shape (n_rows, n_outputs)
        Local effect of each row, i.e. its prediction difference. One column for
        regression, one per class/task for classification or multioutput regression.

    bins : ndarray of shape (n_rows,), dtype=int
        Index of the bin of each row, in `[0, n_bins)`.

    sample_weight : ndarray of shape (n_rows,)
        Sample weight of each row. With unit weights, the result is the mean.

    n_bins : int
        Number of bins. Bins without rows get a mean of 0.

    Returns
    -------
    means : ndarray of shape (n_outputs, n_bins)
        Mean local effect in each bin, for each output.

    bin_weights : ndarray of shape (n_bins,)
        Sum of the weights in each bin, i.e. the number of samples when all
        weights are 1. It is used by the caller to center the effects.
    """
    bin_weights = np.bincount(bins, weights=sample_weight, minlength=n_bins)
    non_empty = bin_weights > 0
    means = np.zeros((local_effects.shape[1], n_bins))
    for output in range(local_effects.shape[1]):
        sums = np.bincount(
            bins, weights=sample_weight * local_effects[:, output], minlength=n_bins
        )
        means[output, non_empty] = sums[non_empty] / bin_weights[non_empty]
    return means, bin_weights


def _ale_continuous(estimator, X, column, x, sample_weight, edges, response_method):
    """First-order ALE of a continuous feature, evaluated at the bin `edges`.

    Each sample is assigned to the limits of the interval it falls into:
    `(edges[k], edges[k + 1]]`. The local effect of a sample is the difference
    of the predictions when its feature is set to the upper and to the lower
    bound of its interval. The local effects are averaged within each interval,
    accumulated from the first interval, and centered so that their weighted mean
    over the samples is zero.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data.

    column : int
        Index of the feature of interest in `X`.

    x : ndarray of shape (n_samples,)
        Values of the feature of interest.

    sample_weight : ndarray of shape (n_samples,)
        Sample weights.

    edges : ndarray of shape (n_bins + 1,)
        Sorted bounds of the intervals.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions.

    Returns
    -------
    ale : ndarray of shape (n_outputs, n_bins + 1)
        Centered accumulated local effects at each value of `edges`.

    bin_sizes : ndarray of shape (n_bins,)
        Sum of the sample weights in each interval.
    """
    n_bins = edges.shape[0] - 1
    bins = np.clip(np.searchsorted(edges, x, side="left"), 1, n_bins) - 1

    pred_upper = _replace_values_and_predict(
        estimator, X, column, edges[bins + 1], response_method
    )
    pred_lower = _replace_values_and_predict(
        estimator, X, column, edges[bins], response_method
    )

    local_effects, bin_sizes = _mean_local_effect_per_bin(
        pred_upper - pred_lower, bins, sample_weight, n_bins
    )
    ale = np.hstack(
        [np.zeros((local_effects.shape[0], 1)), np.cumsum(local_effects, axis=1)]
    )
    midpoints = (ale[:, :-1] + ale[:, 1:]) / 2
    ale -= np.average(midpoints, axis=1, weights=bin_sizes)[:, np.newaxis]
    return ale, bin_sizes


# TODO: decide if we want to support categorical/discrete calculation of ale effects
# in short, this is not described in the paper, it's implemented in R by the author
# of the ale paper, but I am not convinced it makes a lot of sense.
def _ale_categorical(estimator, X, column, x, sample_weight, levels, response_method):
    """First-order ALE of a discrete or categorical feature with ordered `levels`.

    For discrete and categorical features, the local effect of a sample is
    calculated by setting its value to the value immediately smaller and to the
    value immediately larger. This gives two local effects per sample: the
    prediction at the larger value minus the prediction at its original value,
    and the prediction at its original value minus the prediction at the smaller
    value. The exception are the samples at the limits (i.e., the smallest and
    the largest values), which can only be set to the value immediately larger
    or smaller, respectively, so they have only one local effect.

    Each local effect describes the change from one value to the next one, and
    the local effects are averaged for each of these changes. So for example, if
    the variable contains values 1, 3 and 7, there are two changes: from 1 to 3,
    and from 3 to 7. The change from 1 to 3 is the average of the local effects
    of the samples with value 1 when set to 3, and of the samples with value 3
    when set to 1. The change from 3 to 7 is the average of the local effects of
    the samples with value 3 when set to 7, and of the samples with value 7 when
    set to 3.

    Finally, the average changes are added up from the smallest value to the
    largest, which gives the accumulated local effect at each value. The result
    is then shifted up or down so that its average over the samples is zero.


    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator.

    X : {ndarray, dataframe} of shape (n_samples, n_features)
        The data.

    column : int
        Index of the feature of interest in `X`.

    x : ndarray of shape (n_samples,)
        Values of the feature of interest.

    sample_weight : ndarray of shape (n_samples,)
        Sample weights.

    levels : ndarray of shape (n_levels,)
        Unique values of the feature.

    response_method : str or list of str
        The method(s) of `estimator` used to compute the predictions.

    Returns
    -------
    ale : ndarray of shape (n_outputs, n_levels)
        Centered accumulated local effects at each value of `levels`.

    level_sizes : ndarray of shape (n_levels,)
        Sum of the sample weights at each level.
    """
    n_levels = levels.shape[0]
    sorter = np.argsort(levels)
    level_idx = sorter[np.searchsorted(levels, x, sorter=sorter)]

    pred = _replace_values_and_predict(
        estimator, X, column, levels[level_idx], response_method
    )
    pred_up = _replace_values_and_predict(
        estimator,
        X,
        column,
        levels[np.minimum(level_idx + 1, n_levels - 1)],
        response_method,
    )
    pred_down = _replace_values_and_predict(
        estimator, X, column, levels[np.maximum(level_idx - 1, 0)], response_method
    )

    # Samples that can move up or down (i.e., not the edges of the
    # unique value range).
    up = level_idx < n_levels - 1
    down = level_idx > 0
    # Step k goes from levels[k] to levels[k + 1]: moving up from k, down from k + 1.
    steps = np.concatenate([level_idx[up], level_idx[down] - 1])
    # Average the prediction differences of both groups within each step.
    local_effects, _ = _mean_local_effect_per_bin(
        np.concatenate([(pred_up - pred)[up], (pred - pred_down)[down]]),
        steps,
        np.concatenate([sample_weight[up], sample_weight[down]]),
        n_levels - 1,
    )
    # Accumulate the step effects from the first level, which starts at 0.
    ale = np.hstack(
        [np.zeros((local_effects.shape[0], 1)), np.cumsum(local_effects, axis=1)]
    )
    # Number of samples (sum of weights) at each level.
    level_sizes = np.bincount(level_idx, weights=sample_weight, minlength=n_levels)
    # Center so that the average effect over the samples is zero.
    ale -= np.average(ale, axis=1, weights=level_sizes)[:, np.newaxis]
    return ale, level_sizes


@validate_params(
    {
        "estimator": [
            HasMethods(["fit", "predict"]),
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
        ],
        "X": ["array-like"],
        "features": ["array-like", Integral, str],
        "sample_weight": ["array-like", None],
        "categorical_features": ["array-like", None],
        "feature_names": ["array-like", None],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function"})],
        "grid_resolution": [Interval(Integral, 2, None, closed="left")],
        "custom_values": [dict, None],
    },
    prefer_skip_nested_validation=True,
)
def accumulated_local_effects(
    estimator,
    X,
    features,
    *,
    sample_weight=None,
    categorical_features=None,
    feature_names=None,
    response_method="auto",
    grid_resolution=100,
    custom_values=None,
):
    """Accumulated local effects (ALE) of ``features``.

    The range of the feature is split into intervals. Within each interval, the
    feature of the samples that fall into it is set to the lower and upper
    interval bounds, and the difference in the predictions is averaged. These
    local effects are accumulated over the intervals and centered so that their
    weighted mean is zero [1]_.

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
        The feature (e.g. `[0]`) or pair of interacting features
        (e.g. `[(0, 1)]`) for which the partial dependency should be computed.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights used to compute the interval bounds, the average local
        effects and the centering. If `None`, samples are equally weighted.

    categorical_features : array-like of shape (n_features,) or shape \
            (n_categorical_features,), dtype={bool, int, str}, default=None
        Indicates the categorical or discrete features.

        - `None`: no feature will be considered categorical;
        - boolean array-like: boolean mask of shape `(n_features,)`
          indicating which features are categorical;
        - integer or string array-like: integer indices or strings
          indicating categorical features.

        For these features, the local effect between two consecutive observed
        values is computed by moving the samples at the lower value to the
        upper value and the samples at the upper value to the lower value.

        .. warning::

            Categories are ordered as returned by :func:`numpy.unique`, unless
            an order is given in `custom_values`. For nominal features, this
            order is arbitrary and the accumulated local effects depend on it.

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
        Maximum number of interval bounds for continuous features. The bounds
        are the `grid_resolution` quantiles of the feature computed with the
        `"inverted_cdf"` method, so they are always observed values and no
        interval is empty. Repeated quantiles are merged, so fewer bounds may
        be returned. Ignored for categorical features and for features in
        `custom_values`.

    custom_values : dict, default=None
        A dictionary mapping a feature to the interval bounds (continuous
        features) or to the categories in the order in which their effects are
        accumulated (categorical features). Samples outside the bounds or not
        in the categories are ignored.

    Returns
    -------
    result : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        effects : ndarray of shape (n_outputs, n_grid_values)
            The accumulated local effects at each value in `grid_values`.

        grid_values : list of 1d ndarrays
            One array per feature in `features`: the interval bounds for
            continuous features, or the categories for categorical features.

        bin_sizes : ndarray
            The (weighted) number of samples in each interval, of shape
            `(n_grid_values - 1,)`, or in each category, of shape
            `(n_grid_values,)`.

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

    # TODO: fix negative sample weights are silently accepted. Consider
    # `_check_sample_weight(..., ensure_non_negative=True)`.
    sample_weight = _check_sample_weight(sample_weight, X)

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

    is_categorical = _get_is_categorical(
        categorical_features, features_indices, feature_names, X.shape[1]
    )

    custom_values = custom_values or {}
    # TODO: fix numpy integer scalars (e.g. `np.int64(0)`) are not wrapped in a list
    # and raise a TypeError below. `partial_dependence` has the same bug.
    if isinstance(features, (str, int)):
        features = [features]

    X_subset = _safe_indexing(X, features_indices, axis=1)

    custom_values_for_X_subset = {
        index: custom_values.get(feature)
        for index, feature in enumerate(features)
        if feature in custom_values
    }

    # One 1d array with the values of each feature of interest.
    columns = [
        np.asarray(_safe_indexing(X_subset, feature, axis=1))
        for feature in range(len(is_categorical))
    ]

    # First, find the samples to ignore, looking at all the features of interest.
    custom_axes = {}
    for feature, (x, is_cat) in enumerate(zip(columns, is_categorical)):
        # Ignore samples with a missing value in this feature.
        sample_weight = np.where(_get_mask(x, np.nan), 0.0, sample_weight)
        # Nothing else to ignore unless the user gave custom values for this feature.
        if feature in custom_values_for_X_subset:
            # TODO: fix empty custom values raise an IndexError, and 2D custom values
            # are silently flattened. Validate them as `partial_dependence` does.
            axis = np.asarray(custom_values_for_X_subset[feature])
            if is_cat:
                # The custom categories define the order, so they can't repeat.
                if np.unique(axis).shape[0] != axis.shape[0]:
                    raise ValueError(
                        "The custom values of a categorical feature must be unique."
                    )
                # Samples whose category is not in the custom categories.
                outside = ~np.isin(x, axis)
            else:
                # Custom interval bounds, sorted and without repeats.
                axis = np.unique(axis)
                # Samples below the first bound or above the last bound.
                outside = (x < axis[0]) | (x > axis[-1])
            # Ignore those samples.
            sample_weight = np.where(outside, 0.0, sample_weight)
            # Keep the custom grid so that it is not built again below.
            custom_axes[feature] = axis

    # Then, build the grid of each feature from the samples that are left.
    values = []
    for feature, (x, is_cat) in enumerate(zip(columns, is_categorical)):
        if feature in custom_axes:
            # Custom values: use them as they are.
            axis = custom_axes[feature]
        elif is_cat:
            # Categorical: the categories observed in the samples that are left.
            axis = np.unique(x[sample_weight > 0])
        else:
            # TODO: fix string features not listed in `categorical_features` fail
            # here with "could not convert string to float". Raise a clear error.
            # Continuous: `grid_resolution` weighted quantiles ("inverted_cdf"),
            # which are always observed values.
            axis = _weighted_percentile(
                x.astype(np.float64),
                sample_weight,
                np.linspace(0, 100, grid_resolution),
            )
            # Merge repeated quantiles and keep the dtype of the feature.
            axis = np.unique(axis).astype(x.dtype)
        values.append(axis)

    # TODO: fix a single custom bound ignores all samples, so this error is raised
    # instead of the "at least two distinct values" one below.
    if not np.any(sample_weight > 0):
        raise ValueError(
            "No sample with a positive weight is available for features "
            f"{[feature_names[idx] for idx in features_indices]!r}."
        )
    for feature_idx, axis in zip(features_indices, values):
        if axis.shape[0] < 2:
            raise ValueError(
                f"Feature {feature_names[feature_idx]!r} must have at least two "
                "distinct values to compute accumulated local effects. Got "
                f"{axis!r}."
            )

    if features_indices.size == 1:
        column = features_indices[0]
        x = np.asarray(_safe_indexing(X_subset, 0, axis=1))
        x = np.where(sample_weight > 0, x, values[0][0])
        if is_categorical[0]:
            effects, bin_sizes = _ale_categorical(
                estimator, X, column, x, sample_weight, values[0], response_method
            )
        else:
            effects, bin_sizes = _ale_continuous(
                estimator, X, column, x, sample_weight, values[0], response_method
            )
    else:
        # TODO: add support for 2 feature interaction
        raise NotImplementedError(
            "Second-order accumulated local effects of two features are not "
            "supported yet."
        )

    return Bunch(effects=effects, grid_values=values, bin_sizes=bin_sizes)
