"""Metric per threshold curve to assess binary classification performance.

Compute metric per threshold, over a range of threshold values to aid visualization
of threshold-dependent metric behavior.

Utilizes `_CurveScorer` methods to do all the computation.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._param_validation import validate_params


@validate_params(
    {
        "score_func": [callable],
        "y_true": ["array-like"],
        "y_score": ["array-like"],
        "greater_is_better": ["boolean"],
        "labels": ["array-like", None],
    },
    prefer_skip_nested_validation=True,
)
def decision_threshold_curve(
    score_func,
    y_true,
    y_score,
    greater_is_better=True,
    labels=None,
    **kwargs,
):
    """Compute threshold-dependent metric of interest per threshold.

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <threshold_tunning>`.

    .. versionadded:: 1.8

    Parameters
    ----------
    score_func : callable
        The score function to use. It will be called as
        `score_func(y_true, y_pred, **kwargs)`.

    y_true : array-like of shape (n_samples,)
        Ground truth (correct) target labels.

    y_score : array-like of shape (n_samples,)
        Continuous response scores.

    greater_is_better : bool, default=True
        Whether `score_func` is a score function (default), meaning high is
        good, or a loss function, meaning low is good. In the latter case, the
        the output of `score_func` will be sign-flipped.

    labels : array-like, default=None
        Class labels. If `None`, inferred from `y_true`.
        TODO: used `labels` instead of `classes` to be consistent with other metrics.

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

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import accuracy_score, decision_threshold_curve
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_score = np.array([0.1, 0.4, 0.35, 0.8])
    >>> score_thresholds, thresholds = decision_threshold_curve(
    ...     accuracy_score, y_true, y_score)
    >>> thresholds
    array([0.1, 0.33333333, 0.56666667, 0.8 ])
    >>> score_thresholds
    array([0.5, 0.75, 0.75, 0.75])
    """
    # To prevent circular import
    from sklearn.metrics._scorer import _CurveScorer

    sign = 1 if greater_is_better else -1
    curve_scorer = _CurveScorer(score_func, sign, {}, thresholds=None)
    return curve_scorer._scores_from_predictions(
        y_true,
        y_score,
        labels,
        **kwargs,
    )
