"""Accumulated local effects (ALE) plots for regression and classification models."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections.abc import Iterable
import warnings
import numpy as np

from sklearn.base import is_classifier, is_regressor
from sklearn.inspection._pd_utils import _check_feature_names, _get_feature_index
from sklearn.utils import Bunch, _safe_indexing, check_array
from sklearn.utils._indexing import _get_column_indices
from sklearn.utils._param_validation import (
    HasMethods,
    Integral,
    Interval,
    StrOptions,
    validate_params,
)
from sklearn.utils._response import _get_response_values
from sklearn.utils.validation import _check_sample_weight, check_is_fitted

__all__ = [
    "accumulated_local_effects",
]


def _compute_1d_ale(
    estimator,
    X,
    feature_idx,
    response_method,
    grid_resolution,
    percentiles,
    sample_weight=None,
):
    """Compute 1D accumulated local effects for a single continuous feature."""
    x_feat = np.asarray(X[:, feature_idx])
    n_samples = X.shape[0]

    # Compute empirical quantiles for interval boundaries
    q_grid = np.linspace(percentiles[0], percentiles[1], grid_resolution + 1)
    breaks = np.unique(np.quantile(x_feat, q_grid))
    n_intervals = len(breaks) - 1

    if n_intervals < 1:
        warnings.warn(
            f"Feature index {feature_idx} has too few unique values to compute "
            "multiple intervals. Returning zero ALE effect.",
            UserWarning,
        )
        # Handle degenerate case: constant feature
        breaks = np.array([x_feat.min(), x_feat.max()]) if x_feat.size > 0 else np.array([0.0, 1.0])
        # Obtain target dimension by calling prediction once
        pred_sample, _ = _get_response_values(estimator, X[:1], response_method=response_method)
        n_targets = pred_sample.shape[1] if pred_sample.ndim > 1 else 1
        return np.zeros((n_targets, len(breaks))), breaks

    # Assign samples to intervals: bin k is [breaks[k], breaks[k+1])
    indices = np.digitize(x_feat, breaks, right=False) - 1
    indices = np.clip(indices, 0, n_intervals - 1)

    # Initial prediction to determine output dimension (n_targets)
    pred_init, _ = _get_response_values(estimator, X[:1], response_method=response_method)
    if pred_init.ndim == 1:
        n_targets = 1
    else:
        n_targets = pred_init.shape[1]

    deltas = np.zeros((n_targets, n_intervals), dtype=np.float64)

    for k in range(n_intervals):
        in_bin = (indices == k)
        if not np.any(in_bin):
            continue

        # Extract samples falling within the interval
        X_bin = X[in_bin].copy()
        
        # Evaluate predictions at left interval boundary
        X_bin[:, feature_idx] = breaks[k]
        pred_left, _ = _get_response_values(estimator, X_bin, response_method=response_method)
        if pred_left.ndim == 1:
            pred_left = pred_left[np.newaxis, :]
        else:
            pred_left = pred_left.T

        # Evaluate predictions at right interval boundary
        X_bin[:, feature_idx] = breaks[k + 1]
        pred_right, _ = _get_response_values(estimator, X_bin, response_method=response_method)
        if pred_right.ndim == 1:
            pred_right = pred_right[np.newaxis, :]
        else:
            pred_right = pred_right.T

        diff = pred_right - pred_left  # shape (n_targets, n_samples_in_bin)

        if sample_weight is not None:
            sw_bin = sample_weight[in_bin]
            sw_sum = np.sum(sw_bin)
            if sw_sum > 0:
                deltas[:, k] = np.sum(diff * sw_bin, axis=1) / sw_sum
            else:
                deltas[:, k] = np.mean(diff, axis=1)
        else:
            deltas[:, k] = np.mean(diff, axis=1)

    # Accumulate local differences
    ale_uncentered = np.zeros((n_targets, len(breaks)), dtype=np.float64)
    ale_uncentered[:, 1:] = np.cumsum(deltas, axis=1)

    # Center ALE so the weighted average across sample distribution is zero
    sample_ale = np.zeros((n_targets, n_samples), dtype=np.float64)
    for k in range(n_intervals):
        in_bin = (indices == k)
        if not np.any(in_bin):
            continue
        interval_width = breaks[k + 1] - breaks[k]
        if interval_width > 0:
            frac = (x_feat[in_bin] - breaks[k]) / interval_width
        else:
            frac = 0.0
        sample_ale[:, in_bin] = ale_uncentered[:, [k]] + frac * deltas[:, [k]]

    if sample_weight is not None and np.sum(sample_weight) > 0:
        mean_ale = np.sum(sample_ale * sample_weight, axis=1, keepdims=True) / np.sum(sample_weight)
    else:
        mean_ale = np.mean(sample_ale, axis=1, keepdims=True)

    ale_centered = ale_uncentered - mean_ale

    return ale_centered, breaks


@validate_params(
    {
        "estimator": [HasMethods(["fit", "predict"]), HasMethods(["fit", "predict_proba"]), HasMethods(["fit", "decision_function"])],
        "X": ["array-like"],
        "features": ["array-like", Integral, str],
        "sample_weight": ["array-like", None],
        "feature_names": ["array-like", None],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function", "predict"})],
        "grid_resolution": [Interval(Integral, 2, None, closed="left")],
        "percentiles": [tuple],
    },
    prefer_skip_nested_validation=True,
)
def accumulated_local_effects(
    estimator,
    X,
    features,
    *,
    sample_weight=None,
    feature_names=None,
    response_method="auto",
    grid_resolution=100,
    percentiles=(0.0, 1.0),
):
    """Calculate Accumulated Local Effects (ALE) for continuous features.

    Accumulated Local Effects describe how features influence the prediction of an
    estimator on average. Unlike Partial Dependence Plots (PDP), ALE computes
    differences in predictions over conditional local intervals rather than marginal
    averages. This prevents evaluation of unlikely or impossible data instances
    when features are strongly correlated.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted scikit-learn estimator implementing `predict`, `predict_proba`,
        or `decision_function`.

    X : {array-like, dataframe} of shape (n_samples, n_features)
        The training or evaluation data.

    features : array-like of {int, str} or int or str
        The feature(s) for which to compute the accumulated local effects.
        Can be an integer index, a feature name, or an iterable of indices/names.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights used to calculate the weighted average of local differences
        and the mean-centering offset.

    feature_names : array-like of shape (n_features,), default=None
        Name of each feature; used if `X` does not define `columns`.

    response_method : {'auto', 'predict_proba', 'decision_function', 'predict'}, default='auto'
        Specifies whether to use `predict_proba`, `decision_function`, or `predict`
        as the target response. For regressors, defaults to `'predict'`. For classifiers,
        `'auto'` prioritizes `predict_proba` over `decision_function`.

    grid_resolution : int, default=100
        The target number of intervals used to partition each continuous feature.

    percentiles : tuple of float, default=(0.0, 1.0)
        The lower and upper percentiles used to define the extreme values of the grid.
        Must be between 0.0 and 1.0 with percentiles[0] < percentiles[1].

    Returns
    -------
    result : Bunch
        Dictionary-like object with the following attributes:

        ale : list of ndarray of shape (n_targets, n_grid_points)
            The centered accumulated local effects for each requested feature.
        values : list of 1D ndarray
            The grid break points for each feature.
        feature_names : list of str
            The names corresponding to each requested feature.

    References
    ----------
    .. [1] Apley, D. W., & Zhu, J. (2020). "Visualizing the effects of predictor
           variables in black box supervised learning models". Journal of the Royal
           Statistical Society: Series B (Statistical Methodology), 82(4), 1059-1086.
    """
    check_is_fitted(estimator)

    if not isinstance(percentiles, Iterable) or len(percentiles) != 2:
        raise ValueError("'percentiles' must be a sequence of 2 elements.")
    if not (0.0 <= percentiles[0] < percentiles[1] <= 1.0):
        raise ValueError("'percentiles' values must be in [0, 1] with percentiles[0] < percentiles[1].")

    if sample_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, X)

    # Process feature names
    feature_names_ = _check_feature_names(X, feature_names)

    # Convert X to 2D numpy array for fast column indexing and copies
    if hasattr(X, "iloc"):
        # Pandas or Polars DataFrame
        X_arr = np.asarray(X.to_numpy(), dtype=np.float64)
    else:
        X_arr = check_array(X, dtype=np.float64, ensure_all_finite="allow-nan")

    # Normalize `features` parameter
    if isinstance(features, (int, str)):
        features_list = [features]
    else:
        features_list = list(features)

    # Resolve feature indices
    feature_indices = [
        _get_feature_index(f, feature_names=feature_names_) for f in features_list
    ]
    resolved_names = [
        feature_names_[idx] if feature_names_ is not None else f"Feature {idx}"
        for idx in feature_indices
    ]

    # Resolve response method
    if is_regressor(estimator):
        if response_method == "auto":
            resolved_response = "predict"
        else:
            resolved_response = response_method
    elif is_classifier(estimator):
        if response_method == "auto":
            resolved_response = ["predict_proba", "decision_function"]
        else:
            resolved_response = response_method
    else:
        resolved_response = "predict" if response_method == "auto" else response_method

    ale_list = []
    values_list = []

    for idx in feature_indices:
        ale_feat, values_feat = _compute_1d_ale(
            estimator=estimator,
            X=X_arr,
            feature_idx=idx,
            response_method=resolved_response,
            grid_resolution=grid_resolution,
            percentiles=percentiles,
            sample_weight=sample_weight,
        )
        ale_list.append(ale_feat)
        values_list.append(values_feat)

    return Bunch(
        ale=ale_list,
        values=values_list,
        grid_values=values_list,
        feature_names=resolved_names,
    )
