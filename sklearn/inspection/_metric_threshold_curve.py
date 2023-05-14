"""Metrics per threshold curves are used to assess performance on binary
classification task given threshold grid. One can undestand the behaviour of
threshold-dependent metrics when changing the threshold.
"""

# Authors: ########
# License: BSD 3 clause

from numbers import Real, Integral

import numpy as np

from ..utils import assert_all_finite
from ..utils import check_consistent_length
from ..utils.validation import (
    _check_pos_label_consistency,
    _check_sample_weight,
)
from ..utils import column_or_1d
from ..utils.multiclass import type_of_target
from ..utils._param_validation import validate_params, Interval


@validate_params(
    {
        "y_true": ["array-like"],
        "y_score": ["array-like"],
        "score_func": [callable],
        "threshold_grid": [
            Interval(Integral, 3, None, closed="left"),
            "array-like",
            None,
        ],
        "pos_label": [Real, str, "boolean", None],
        "sample_weight": ["array-like", None],
    }
)
def metric_threshold_curve(
    y_true,
    y_score,
    score_func,
    *,
    threshold_grid=101,
    pos_label=None,
    sample_weight=None,
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

    score_func : callable
        Threshold dependent score function (or loss function) with signature
        `score_func(y, y_pred, sample_weight, **kwargs)`.

    threshold_grid : array-like, int or None, default=101
        Values of threhsold for each score calculation. If int then
        `threshold_grid` percentiles of `y_score` are selected. If `None` then
        all possible thresholds are selected. If int is lower then
        `len(set(y_score))` then all possible thresholds are selected.

    pos_label : int, float, bool or str, default=None
        The label of the positive class.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

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
    """
    # Check to make sure y_true is valid
    y_type = type_of_target(y_true, input_name="y_true")
    if not (y_type == "binary" or (y_type == "multiclass" and pos_label is not None)):
        raise ValueError("{0} format is not supported".format(y_type))

    check_consistent_length(y_true, y_score, sample_weight)
    y_true = column_or_1d(y_true)
    y_score = column_or_1d(y_score)
    assert_all_finite(y_true)
    assert_all_finite(y_score)

    # Filter out zero-weighted samples, as they should not impact the result
    if sample_weight is not None:
        sample_weight = column_or_1d(sample_weight)
        sample_weight = _check_sample_weight(sample_weight, y_true)
        nonzero_weight_mask = sample_weight != 0
        y_true = y_true[nonzero_weight_mask]
        y_score = y_score[nonzero_weight_mask]
        sample_weight = sample_weight[nonzero_weight_mask]

    pos_label = _check_pos_label_consistency(pos_label, y_true)

    # make y_true a boolean vector
    y_true = y_true == pos_label

    # sort scores and corresponding truth values
    desc_score_indices = np.argsort(y_score, kind="mergesort")[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]
    if sample_weight is not None:
        sample_weight = sample_weight[desc_score_indices]

    # logic to see if we need to use all possible thresholds (distinct values)
    all_thresholds = False
    if threshold_grid is None:
        all_thresholds = True
    elif isinstance(threshold_grid, int):
        if len(set(y_score)) < threshold_grid:
            all_thresholds = True

    if all_thresholds:
        # y_score typically has many tied values. Here we extract
        # the indices associated with the distinct values. We also
        # concatenate a value for the end of the curve.
        distinct_value_indices = np.where(np.diff(y_score))[0]
        threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]
        thresholds = y_score[threshold_idxs[::-1]]
    elif isinstance(threshold_grid, int):
        # takes representative score points to calculate the metric
        # with these thresholds
        thresholds = np.percentile(
            list(set(y_score)), np.linspace(0, 100, threshold_grid)
        )
    else:
        # if threshold_grid is an array then run some checks and sort
        # it for consistency
        threshold_grid = column_or_1d(threshold_grid)
        assert_all_finite(threshold_grid)
        thresholds = np.sort(threshold_grid)

    # for each threshold calculates the metric
    metric_values = []
    for threshold in thresholds:
        preds_threshold = (y_score > threshold).astype(int)
        metric_values.append(
            score_func(y_true, preds_threshold, sample_weight=sample_weight)
        )
    # TODO: should we multithread the metric calculations?

    return np.array(metric_values), thresholds
