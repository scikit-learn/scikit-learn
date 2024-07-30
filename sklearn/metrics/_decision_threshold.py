"""Metric per threshold curve to assess binary classification performance.

Given threshold grid, one can undestand the behaviour of threshold-dependent
metrics when changing the threshold. In imbalanced scenarios or
cost-sensitive learning, a 0.5 threshold may not be optimal and tools like
this can help you visualize how the performance changes.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


from numbers import Integral

import numpy as np

from ..utils import assert_all_finite, check_consistent_length, column_or_1d
from ..utils._param_validation import Interval, validate_params
from ..utils.multiclass import type_of_target
from ..utils.validation import (
    _check_pos_label_consistency,
    _check_sample_weight,
)


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
        "scoring_kwargs": [dict, None],
    },
    prefer_skip_nested_validation=True,
)
def decision_threshold_curve(
    y_true,
    y_score,
    scoring,
    *,
    thresholds=100,
    scoring_kwargs=None,
):
    """Compute the threshold-dependent metric of interest per threshold.

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <metric_threshold_curve>`.

    .. versionadded:: 1.3

    Parameters
    ----------
    y_true : array-like of shape (n_samples,), default=None
        True targets of binary classification.

    y_score : array-like of shape (n_samples,), default=None
        Estimated probabilities or output of a decision function.

    scoring : callable
        Threshold-dependent score function (or loss function) with signature
        `scoring(y, y_pred, **scoring_kwargs)`.

    thresholds : array-like or int, default=101
        Values of threhsold for each score calculation. If int then
        `thresholds` percentiles of `y_score` are selected. If int is lower
        then `len(set(y_score))` then all possible thresholds are selected.

    scoring_kwargs : dict, default=None
        Keyword arguments to pass to specified `scoring` function.

    Returns
    -------
    metric_values : ndarray of shape (n_thresholds,)
        Score value for each threshold. At index i being the value of the
        theshold-dependent metric for predictions score >= thresholds[i].

    thresholds : ndarray of shape (n_thresholds,)
        Ascending score values used as thresholds.

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
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> accuracy_values, thresholds = decision_threshold_curve(
    ...     y_true, y_scores, accuracy_score)
    >>> thresholds
    array([0.1 , 0.35, 0.4 , 0.8 ])
    >>> accuracy_values
    array([0.75, 0.5 , 0.75, 0.5 ])
    """
    if scoring_kwargs is None:
        scoring_kwargs = {}

    # Check to make sure y_true is valid.
    y_type = type_of_target(y_true, input_name="y_true")
    pos_label = scoring_kwargs.get("pos_label")
    if not (y_type == "binary" or (y_type == "multiclass" and pos_label is not None)):
        if y_type == "multiclass":
            raise ValueError(
                "In a multiclass scenario, you must pass a `pos_label` \
                to `scoring_kwargs`."
            )
        raise ValueError("{0} format is not supported".format(y_type))

    sample_weight = scoring_kwargs.get("sample_weight")
    check_consistent_length(y_true, y_score, sample_weight)
    y_true = column_or_1d(y_true)
    y_score = column_or_1d(y_score)
    assert_all_finite(y_true)
    assert_all_finite(y_score)

    # Filter out zero-weighted samples, as they should not impact the result.

    if sample_weight is not None:
        sample_weight = column_or_1d(sample_weight)
        sample_weight = _check_sample_weight(sample_weight, y_true)
        nonzero_weight_mask = sample_weight != 0
        y_true = y_true[nonzero_weight_mask]
        y_score = y_score[nonzero_weight_mask]
        sample_weight = sample_weight[nonzero_weight_mask]

    pos_label = _check_pos_label_consistency(pos_label, y_true)

    # Make y_true a boolean vector.
    y_true = y_true == pos_label

    # Sort scores and corresponding truth values.
    desc_score_indices = np.argsort(y_score, kind="mergesort")[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]
    if sample_weight is not None:
        sample_weight = sample_weight[desc_score_indices]

    if "sample_weight" in scoring_kwargs:
        scoring_kwargs["sample_weight"] = sample_weight

    # Logic to see if we need to use all possible thresholds (distinct values).
    all_thresholds = isinstance(thresholds, int) and len(set(y_score)) < thresholds

    if all_thresholds:
        # y_score typically has many tied values. Here we extract
        # the indices associated with the distinct values. We also
        # concatenate a value for the end of the curve.
        distinct_value_indices = np.where(np.diff(y_score))[0]
        threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]
        thresholds = y_score[threshold_idxs[::-1]]
    elif isinstance(thresholds, int):
        # It takes representative score points to calculate the metric
        # with these thresholds.
        thresholds = np.percentile(list(set(y_score)), np.linspace(0, 100, thresholds))
    else:
        # If thresholds is an array then run some checks and sort
        # it for consistency.
        thresholds = column_or_1d(thresholds)
        assert_all_finite(thresholds)
        thresholds = np.sort(thresholds)

    # For each threshold calculates the metric.
    metric_values = []
    for threshold in thresholds:
        preds_threshold = (y_score > threshold).astype(int)
        metric_values.append(scoring(y_true, preds_threshold, **scoring_kwargs))
    # TODO: should we multithread the metric calculations?

    return np.array(metric_values), thresholds
