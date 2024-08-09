"""Metric per threshold curve to assess binary classification performance.

Given threshold grid, one can undestand the behaviour of threshold-dependent
metrics when changing the threshold. In imbalanced scenarios or
cost-sensitive learning, a 0.5 threshold may not be optimal and tools like
this can help you visualize how the performance changes.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Integral

from ..utils._param_validation import Interval, validate_params
from ._scorer import _CurveScorer


@validate_params(
    {
        "y_true": ["array-like"],
        "y_score": ["array-like"],
        "scoring": [callable],
        "thresholds": [
            Interval(Integral, 3, None, closed="left"),
            "array-like",
            None,
        ],
    },
    prefer_skip_nested_validation=True,
)
def decision_threshold_curve(
    y_true,
    y_score,
    scoring,
    thresholds=100,
):
    """Compute the threshold-dependent metric of interest per threshold.

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <metric_threshold_curve>`.

    .. versionadded:: 1.6

    Parameters
    ----------
    y_true : array-like of shape (n_samples,), default=None
        True targets of binary classification.

    y_score : array-like of shape (n_samples,), default=None
        Estimated probabilities or output of a decision function.

    scoring : callable, default=None
        The objective metric to be estimated. It should be a callable object created
        with :func:`~sklearn.metrics.make_scorer`.
        # TODO(Carlo): Change it to also just be a function callable. In this case,
        # transform it in a scorer inside the function.

    thresholds : int or array-like, default=100
        Related to the number of decision thresholds for which we want to compute the
        score. If an integer, it will be used to generate `thresholds` thresholds
        uniformly distributed between the minimum and maximum of `y_score`. If an
        array-like, it will be used as the thresholds.

    Returns
    -------
    metric_values : ndarray of shape (n_thresholds,)
        The scores associated to each threshold. At index i being the value of the
        theshold-dependent metric for predictions score >= thresholds[i].
        # TODO(Carlo) Check if > or >=

    thresholds : ndarray of shape (n_thresholds,)
        Ascending score values used as thresholds.

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
    # if scoring is function ... transform into scorer (do I need an estimator?)
    curve_scorer = _CurveScorer.from_scorer(scoring, thresholds)
    metric_values, thresholds = curve_scorer._score_given_prediction(y_score)

    return metric_values, thresholds
