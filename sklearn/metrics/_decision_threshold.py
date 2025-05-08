"""Metric per threshold curve to assess binary classification performance.

Given threshold grid, one can undestand the behaviour of threshold-dependent
metrics when changing the threshold. In imbalanced scenarios or
cost-sensitive learning, a 0.5 threshold may not be optimal and tools like
this can help you visualize how the performance changes.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Integral, Real

from ..utils._param_validation import Interval, Options, validate_params


@validate_params(
    {
        "scoring_function": [callable],
        "y_true": ["array-like"],
        "y_score": ["array-like"],
        "thresholds": [
            Interval(Integral, 2, None, closed="left"),
            "array-like",
        ],
        "sign": Options(Real, {0, 1}),
        "labels": ["array-like", None],
        "pos_label": [Real, str, "boolean", None],
    },
    prefer_skip_nested_validation=True,
)
def decision_threshold_curve(
    scoring_function,
    y_true,
    y_score,
    # Should below 2 have a default value?
    thresholds=20,
    sign=1,
    labels=None,
    pos_label=None,
    **kwargs,
):
    """Compute threshold-dependent metric of interest per threshold.

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <metric_threshold_curve>`.

    .. versionadded:: 1.8

    Parameters
    ----------
    scoring_function : callable
        The score function to use. It will be called as
        `score_func(y_true, y_pred, **kwargs)`.
        TODO: decided on `scoring_function` as term also used in forest estimators

    y_true : array-like of shape (n_samples,)
        Ground truth (correct) target labels.

    y_score : array-like of shape (n_samples,)
        Continuous response scores.

    thresholds : int or array-like, default=20
        Specifies number of decision thresholds to compute score for. If an integer,
        it will be used to generate `thresholds` thresholds uniformly distributed
        between the minimum and maximum of `y_score`. If an array-like, it will be
        used as the thresholds.

    sign : int, default=1
        Either 1 or -1. Score is computed as `sign * score_func(estimator, X, y)`.
        Thus, `sign` defines whether higher scores are better or worse.

    labels: array-like, default=None
        Class labels. If `None`, inferred from `y_true`.
        TODO: used `labels` instead of `classes` to be consistent with other metrics.

    pos_label : int, float, bool or str, default=None
        The label of the positive class, used when thresholding `y_score`.
        If `score_func` also has a `pos_label` parameter, this value will also
        be passed `score_func`.
        If `None`, the default value of `score_func(pos_label)`, if present, is
        used. If not present, `1` is used.
        TODO: do we need to allow the user to set this even when `score_func`
        does not take `pos_label`? I think yes, so user can control
        output of `_threshold_scores_to_class_labels`.

    **kwargs : dict
        Parameters to pass to `score_func`.

    Returns
    -------
    score_thresholds : ndarray of shape (n_thresholds,)
        The scores associated with each threshold.

    thresholds : ndarray of shape (n_thresholds,)
        The thresholds used to compute the scores.

    See Also
    --------
    precision_recall_curve : Compute precision-recall pairs for different
        probability thresholds.
    det_curve : Compute error rates for different probability thresholds.
    roc_curve : Compute Receiver operating characteristic (ROC) curve.

    Examples #TODO(Carlo) change the example and fix threshold.
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import accuracy_score, decision_threshold_curve
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_score = np.array([0.1, 0.4, 0.35, 0.8])
    >>> accuracy_values, thresholds = decision_threshold_curve(
    ...     y_true, y_score, accuracy_score)
    >>> thresholds
    array([0.1 , 0.35, 0.4 , 0.8 ])
    >>> accuracy_values
    array([0.75, 0.5 , 0.75, 0.5 ])
    """
    # To prevent circular import
    from ._scorer import _CurveScorer

    curve_scorer = _CurveScorer(scoring_function, sign, {}, thresholds)
    return curve_scorer._scores_from_predictions(
        y_true,
        y_score,
        labels,
        pos_label,
        **kwargs,
    )
