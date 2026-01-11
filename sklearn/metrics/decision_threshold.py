"""Metric per threshold curve to assess binary classification performance.

Compute metric per threshold (unique `y_score`) to aid visualization
of threshold-dependent metric behavior.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.metrics._ranking import _validate_inputs_for_classification_thresholds
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
    metric_func,
    y_true,
    y_score,
    pos_label=None,
    **kwargs,
):
    """Compute `score_func` per threshold.

    Read more in the :ref:`User Guide <threshold_tunning>`.

    .. versionadded:: 1.8

    Parameters
    ----------
    metric_func : callable
        The metric function to use. It will be called as
        `score_func(y_true, y_pred, **kwargs)`.

    y_true : array-like of shape (n_samples,)
        Ground truth (correct) target labels.

    y_score : array-like of shape (n_samples,)
        Continuous response scores.

    pos_label : int, float, bool or str, default=None
        The label of the positive class.

    **kwargs : dict
        Parameters to pass to `score_func`.

    Returns
    -------
    metric_values : ndarray of shape (n_thresholds,)
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
    >>> scores, thresholds = decision_threshold_curve(
    ...     accuracy_score, y_true, y_score)
    >>> thresholds
    array([0.8 , 0.4 , 0.35, 0.1 ])
    >>> scores
    array([0.75, 0.5 , 0.75, 0.5 ])
    """
    # `y_true` gets converted into a boolean vector
    y_true, y_score, threshold_idxs, _ = _validate_inputs_for_classification_thresholds(
        y_true, y_score, pos_label=pos_label
    )
    y_true = y_true.astype(int)

    thresholds = y_score[threshold_idxs]
    metric_values = np.empty(thresholds.shape[0])
    for idx, threshold in enumerate(thresholds):
        thresholded_labels = (y_score >= threshold).astype(np.int32)
        metric_values[idx] = metric_func(y_true, thresholded_labels, **kwargs)

    return metric_values, thresholds
