"""Metrics to assess performance on classification task given scores

Functions named as ``*_score`` return a scalar value to maximize: the higher
the better

Function named as ``*_error`` or ``*_loss`` return a scalar value to minimize:
the lower the better
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Arnaud Joly <a.joly@ulg.ac.be>
#          Jochen Wersdorfer <jochen@wersdoerfer.de>
#          Lars Buitinck
#          Joel Nothman <joel.nothman@gmail.com>
#          Noel Dawe <noel@dawe.me>
#          Fabian Falck <fabian.falck@web.de>
# License: BSD 3 clause


from __future__ import division

import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from math import floor, log
from scipy.sparse import csr_matrix
from scipy.stats import rankdata

from ..utils import assert_all_finite
from ..utils import check_consistent_length
from ..utils import column_or_1d, check_array
from ..utils.multiclass import type_of_target
from ..utils.extmath import stable_cumsum
from ..utils.sparsefuncs import count_nonzero
from ..exceptions import UndefinedMetricWarning
from ..preprocessing import label_binarize

from .base import _average_binary_score


def auc(x, y, reorder='deprecated'):
    """Compute Area Under the Curve (AUC) using the trapezoidal rule

    This is a general function, given points on a curve.  For computing the
    area under the ROC-curve, see :func:`roc_auc_score`.  For an alternative
    way to summarize a precision-recall curve, see
    :func:`average_precision_score`.

    Parameters
    ----------
    x : array, shape = [n]
        x coordinates. These must be either monotonic increasing or monotonic
        decreasing.
    y : array, shape = [n]
        y coordinates.
    reorder : boolean, optional (default='deprecated')
        Whether to sort x before computing. If False, assume that x must be
        either monotonic increasing or monotonic decreasing. If True, y is
        used to break ties when sorting x. Make sure that y has a monotonic
        relation to x when setting reorder to True.

        .. deprecated:: 0.20
           Parameter ``reorder`` has been deprecated in version 0.20 and will
           be removed in 0.22. It's introduced for roc_auc_score (not for
           general use) and is no longer used there. What's more, the result
           from auc will be significantly influenced if x is sorted
           unexpectedly due to slight floating point error (See issue #9786).
           Future (and default) behavior is equivalent to ``reorder=False``.

    Returns
    -------
    auc : float

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([1, 1, 2, 2])
    >>> pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=2)
    >>> metrics.auc(fpr, tpr)
    0.75

    See also
    --------
    roc_auc_score : Compute the area under the ROC curve
    average_precision_score : Compute average precision from prediction scores
    precision_recall_curve :
        Compute precision-recall pairs for different probability thresholds
    """
    check_consistent_length(x, y)
    x = column_or_1d(x)
    y = column_or_1d(y)

    if x.shape[0] < 2:
        raise ValueError('At least 2 points are needed to compute'
                         ' area under curve, but x.shape = %s' % x.shape)

    if reorder != 'deprecated':
        warnings.warn("The 'reorder' parameter has been deprecated in "
                      "version 0.20 and will be removed in 0.22. It is "
                      "recommended not to set 'reorder' and ensure that x "
                      "is monotonic increasing or monotonic decreasing.",
                      DeprecationWarning)

    direction = 1
    if reorder is True:
        # reorder the data points according to the x axis and using y to
        # break ties
        order = np.lexsort((y, x))
        x, y = x[order], y[order]
    else:
        dx = np.diff(x)
        if np.any(dx < 0):
            if np.all(dx <= 0):
                direction = -1
            else:
                raise ValueError("x is neither increasing nor decreasing "
                                 ": {}.".format(x))

    area = direction * np.trapz(y, x)
    if isinstance(area, np.memmap):
        # Reductions such as .sum used internally in np.trapz do not return a
        # scalar by default for numpy.memmap instances contrary to
        # regular numpy.ndarray instances.
        area = area.dtype.type(area)
    return area


def average_precision_score(y_true, y_score, average="macro",
                            sample_weight=None):
    """Compute average precision (AP) from prediction scores

    AP summarizes a precision-recall curve as the weighted mean of precisions
    achieved at each threshold, with the increase in recall from the previous
    threshold used as the weight:

    .. math::
        \\text{AP} = \\sum_n (R_n - R_{n-1}) P_n

    where :math:`P_n` and :math:`R_n` are the precision and recall at the nth
    threshold [1]_. This implementation is not interpolated and is different
    from computing the area under the precision-recall curve with the
    trapezoidal rule, which uses linear interpolation and can be too
    optimistic.

    Note: this implementation is restricted to the binary classification task
    or multilabel classification task.

    Read more in the :ref:`User Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    y_true : array, shape = [n_samples] or [n_samples, n_classes]
        True binary labels (either {0, 1} or {-1, 1}).

    y_score : array, shape = [n_samples] or [n_samples, n_classes]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers).

    average : string, [None, 'micro', 'macro' (default), 'samples', 'weighted']
        If ``None``, the scores for each class are returned. Otherwise,
        this determines the type of averaging performed on the data:

        ``'micro'``:
            Calculate metrics globally by considering each element of the label
            indicator matrix as a label.
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean.  This does not take label imbalance into account.
        ``'weighted'``:
            Calculate metrics for each label, and find their average, weighted
            by support (the number of true instances for each label).
        ``'samples'``:
            Calculate metrics for each instance, and find their average.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    average_precision : float

    References
    ----------
    .. [1] `Wikipedia entry for the Average precision
           <http://en.wikipedia.org/w/index.php?title=Information_retrieval&
           oldid=793358396#Average_precision>`_

    See also
    --------
    roc_auc_score : Compute the area under the ROC curve

    precision_recall_curve :
        Compute precision-recall pairs for different probability thresholds

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import average_precision_score
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> average_precision_score(y_true, y_scores)  # doctest: +ELLIPSIS
    0.83...

    Notes
    -----
    .. versionchanged:: 0.19
      Instead of linearly interpolating between operating points, precisions
      are weighted by the change in recall since the last operating point.
    """
    def _binary_uninterpolated_average_precision(
            y_true, y_score, sample_weight=None):
        precision, recall, _ = precision_recall_curve(
            y_true, y_score, sample_weight=sample_weight)
        # Return the step function integral
        # The following works because the last entry of precision is
        # guaranteed to be 1, as returned by precision_recall_curve
        return -np.sum(np.diff(recall) * np.array(precision)[:-1])

    return _average_binary_score(_binary_uninterpolated_average_precision,
                                 y_true, y_score, average,
                                 sample_weight=sample_weight)


def roc_auc_score(y_true, y_score, average="macro", sample_weight=None,
                  max_fpr=None):
    """Compute Area Under the Receiver Operating Characteristic Curve (ROC AUC)
    from prediction scores.

    Note: this implementation is restricted to the binary classification task
    or multilabel classification task in label indicator format.

    Read more in the :ref:`User Guide <roc_metrics>`.

    Parameters
    ----------
    y_true : array, shape = [n_samples] or [n_samples, n_classes]
        True binary labels or binary label indicators.

    y_score : array, shape = [n_samples] or [n_samples, n_classes]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers). For binary
        y_true, y_score is supposed to be the score of the class with greater
        label.

    average : string, [None, 'micro', 'macro' (default), 'samples', 'weighted']
        If ``None``, the scores for each class are returned. Otherwise,
        this determines the type of averaging performed on the data:

        ``'micro'``:
            Calculate metrics globally by considering each element of the label
            indicator matrix as a label.
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean.  This does not take label imbalance into account.
        ``'weighted'``:
            Calculate metrics for each label, and find their average, weighted
            by support (the number of true instances for each label).
        ``'samples'``:
            Calculate metrics for each instance, and find their average.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    max_fpr : float > 0 and <= 1, optional
        If not ``None``, the standardized partial AUC [3]_ over the range
        [0, max_fpr] is returned.

    Returns
    -------
    auc : float

    References
    ----------
    .. [1] `Wikipedia entry for the Receiver operating characteristic
            <https://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_

    .. [2] Fawcett T. An introduction to ROC analysis[J]. Pattern Recognition
           Letters, 2006, 27(8):861-874.

    .. [3] `Analyzing a portion of the ROC curve. McClish, 1989
            <http://www.ncbi.nlm.nih.gov/pubmed/2668680>`_

    See also
    --------
    average_precision_score : Area under the precision-recall curve

    roc_curve : Compute Receiver operating characteristic (ROC) curve

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import roc_auc_score
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> roc_auc_score(y_true, y_scores)
    0.75

    """
    def _binary_roc_auc_score(y_true, y_score, sample_weight=None):
        if len(np.unique(y_true)) != 2:
            raise ValueError("Only one class present in y_true. ROC AUC score "
                             "is not defined in that case.")

        fpr, tpr, _ = roc_curve(y_true, y_score,
                                sample_weight=sample_weight)
        if max_fpr is None or max_fpr == 1:
            return auc(fpr, tpr)
        if max_fpr <= 0 or max_fpr > 1:
            raise ValueError("Expected max_frp in range ]0, 1], got: %r"
                             % max_fpr)

        # Add a single point at max_fpr by linear interpolation
        stop = np.searchsorted(fpr, max_fpr, 'right')
        x_interp = [fpr[stop - 1], fpr[stop]]
        y_interp = [tpr[stop - 1], tpr[stop]]
        tpr = np.append(tpr[:stop], np.interp(max_fpr, x_interp, y_interp))
        fpr = np.append(fpr[:stop], max_fpr)
        partial_auc = auc(fpr, tpr)

        # McClish correction: standardize result to be 0.5 if non-discriminant
        # and 1 if maximal
        min_area = 0.5 * max_fpr**2
        max_area = max_fpr
        return 0.5 * (1 + (partial_auc - min_area) / (max_area - min_area))

    y_type = type_of_target(y_true)
    if y_type == "binary":
        labels = np.unique(y_true)
        y_true = label_binarize(y_true, labels)[:, 0]

    return _average_binary_score(
        _binary_roc_auc_score, y_true, y_score, average,
        sample_weight=sample_weight)


def _binary_clf_curve(y_true, y_score, pos_label=None, sample_weight=None):
    """Calculate true and false positives per binary classification threshold.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets of binary classification

    y_score : array, shape = [n_samples]
        Estimated probabilities or decision function

    pos_label : int or str, default=None
        The label of the positive class

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    fps : array, shape = [n_thresholds]
        A count of false positives, at index i being the number of negative
        samples assigned a score >= thresholds[i]. The total number of
        negative samples is equal to fps[-1] (thus true negatives are given by
        fps[-1] - fps).

    tps : array, shape = [n_thresholds <= len(np.unique(y_score))]
        An increasing count of true positives, at index i being the number
        of positive samples assigned a score >= thresholds[i]. The total
        number of positive samples is equal to tps[-1] (thus false negatives
        are given by tps[-1] - tps).

    thresholds : array, shape = [n_thresholds]
        Decreasing score values.
    """
    # Check to make sure y_true is valid
    y_type = type_of_target(y_true)
    if not (y_type == "binary" or
            (y_type == "multiclass" and pos_label is not None)):
        raise ValueError("{0} format is not supported".format(y_type))

    check_consistent_length(y_true, y_score, sample_weight)
    y_true = column_or_1d(y_true)
    y_score = column_or_1d(y_score)
    assert_all_finite(y_true)
    assert_all_finite(y_score)

    if sample_weight is not None:
        sample_weight = column_or_1d(sample_weight)

    # ensure binary classification if pos_label is not specified
    classes = np.unique(y_true)
    if (pos_label is None and
        not (np.array_equal(classes, [0, 1]) or
             np.array_equal(classes, [-1, 1]) or
             np.array_equal(classes, [0]) or
             np.array_equal(classes, [-1]) or
             np.array_equal(classes, [1]))):
        raise ValueError("Data is not binary and pos_label is not specified")
    elif pos_label is None:
        pos_label = 1.

    # make y_true a boolean vector
    y_true = (y_true == pos_label)

    # sort scores and corresponding truth values
    desc_score_indices = np.argsort(y_score, kind="mergesort")[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]
    if sample_weight is not None:
        weight = sample_weight[desc_score_indices]
    else:
        weight = 1.

    # y_score typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve.
    distinct_value_indices = np.where(np.diff(y_score))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]

    # accumulate the true positives with decreasing threshold
    tps = stable_cumsum(y_true * weight)[threshold_idxs]
    if sample_weight is not None:
        # express fps as a cumsum to ensure fps is increasing even in
        # the presence of floating point errors
        fps = stable_cumsum((1 - y_true) * weight)[threshold_idxs]
    else:
        fps = 1 + threshold_idxs - tps
    return fps, tps, y_score[threshold_idxs]


def precision_recall_curve(y_true, probas_pred, pos_label=None,
                           sample_weight=None):
    """Compute precision-recall pairs for different probability thresholds

    Note: this implementation is restricted to the binary classification task.

    The precision is the ratio ``tp / (tp + fp)`` where ``tp`` is the number of
    true positives and ``fp`` the number of false positives. The precision is
    intuitively the ability of the classifier not to label as positive a sample
    that is negative.

    The recall is the ratio ``tp / (tp + fn)`` where ``tp`` is the number of
    true positives and ``fn`` the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The last precision and recall values are 1. and 0. respectively and do not
    have a corresponding threshold.  This ensures that the graph starts on the
    y axis.

    Read more in the :ref:`User Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets of binary classification in range {-1, 1} or {0, 1}.

    probas_pred : array, shape = [n_samples]
        Estimated probabilities or decision function.

    pos_label : int or str, default=None
        The label of the positive class

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    precision : array, shape = [n_thresholds + 1]
        Precision values such that element i is the precision of
        predictions with score >= thresholds[i] and the last element is 1.

    recall : array, shape = [n_thresholds + 1]
        Decreasing recall values such that element i is the recall of
        predictions with score >= thresholds[i] and the last element is 0.

    thresholds : array, shape = [n_thresholds <= len(np.unique(probas_pred))]
        Increasing thresholds on the decision function used to compute
        precision and recall.

    See also
    --------
    average_precision_score : Compute average precision from prediction scores

    roc_curve : Compute Receiver operating characteristic (ROC) curve

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import precision_recall_curve
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> precision, recall, thresholds = precision_recall_curve(
    ...     y_true, y_scores)
    >>> precision  # doctest: +ELLIPSIS
    array([0.66666667, 0.5       , 1.        , 1.        ])
    >>> recall
    array([1. , 0.5, 0.5, 0. ])
    >>> thresholds
    array([0.35, 0.4 , 0.8 ])

    """
    fps, tps, thresholds = _binary_clf_curve(y_true, probas_pred,
                                             pos_label=pos_label,
                                             sample_weight=sample_weight)

    precision = tps / (tps + fps)
    recall = tps / tps[-1]

    # stop when full recall attained
    # and reverse the outputs so recall is decreasing
    last_ind = tps.searchsorted(tps[-1])
    sl = slice(last_ind, None, -1)
    return np.r_[precision[sl], 1], np.r_[recall[sl], 0], thresholds[sl]


def roc_curve(y_true, y_score, pos_label=None, sample_weight=None,
              drop_intermediate=True):
    """Compute Receiver operating characteristic (ROC)

    Note: this implementation is restricted to the binary classification task.

    Read more in the :ref:`User Guide <roc_metrics>`.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        True binary labels. If labels are not either {-1, 1} or {0, 1}, then
        pos_label should be explicitly given.

    y_score : array, shape = [n_samples]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers).

    pos_label : int or str, default=None
        Label considered as positive and others are considered negative.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    drop_intermediate : boolean, optional (default=True)
        Whether to drop some suboptimal thresholds which would not appear
        on a plotted ROC curve. This is useful in order to create lighter
        ROC curves.

        .. versionadded:: 0.17
           parameter *drop_intermediate*.

    Returns
    -------
    fpr : array, shape = [>2]
        Increasing false positive rates such that element i is the false
        positive rate of predictions with score >= thresholds[i].

    tpr : array, shape = [>2]
        Increasing true positive rates such that element i is the true
        positive rate of predictions with score >= thresholds[i].

    thresholds : array, shape = [n_thresholds]
        Decreasing thresholds on the decision function used to compute
        fpr and tpr. `thresholds[0]` represents no instances being predicted
        and is arbitrarily set to `max(y_score) + 1`.

    See also
    --------
    roc_auc_score : Compute the area under the ROC curve

    Notes
    -----
    Since the thresholds are sorted from low to high values, they
    are reversed upon returning them to ensure they correspond to both ``fpr``
    and ``tpr``, which are sorted in reversed order during their calculation.

    References
    ----------
    .. [1] `Wikipedia entry for the Receiver operating characteristic
            <https://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_

    .. [2] Fawcett T. An introduction to ROC analysis[J]. Pattern Recognition
           Letters, 2006, 27(8):861-874.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([1, 1, 2, 2])
    >>> scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=2)
    >>> fpr
    array([0. , 0. , 0.5, 0.5, 1. ])
    >>> tpr
    array([0. , 0.5, 0.5, 1. , 1. ])
    >>> thresholds
    array([1.8 , 0.8 , 0.4 , 0.35, 0.1 ])

    """
    fps, tps, thresholds = _binary_clf_curve(
        y_true, y_score, pos_label=pos_label, sample_weight=sample_weight)

    # Attempt to drop thresholds corresponding to points in between and
    # collinear with other points. These are always suboptimal and do not
    # appear on a plotted ROC curve (and thus do not affect the AUC).
    # Here np.diff(_, 2) is used as a "second derivative" to tell if there
    # is a corner at the point. Both fps and tps must be tested to handle
    # thresholds with multiple data points (which are combined in
    # _binary_clf_curve). This keeps all cases where the point should be kept,
    # but does not drop more complicated cases like fps = [1, 3, 7],
    # tps = [1, 2, 4]; there is no harm in keeping too many thresholds.
    if drop_intermediate and len(fps) > 2:
        optimal_idxs = np.where(np.r_[True,
                                      np.logical_or(np.diff(fps, 2),
                                                    np.diff(tps, 2)),
                                      True])[0]
        fps = fps[optimal_idxs]
        tps = tps[optimal_idxs]
        thresholds = thresholds[optimal_idxs]

    if tps.size == 0 or fps[0] != 0 or tps[0] != 0:
        # Add an extra threshold position if necessary
        # to make sure that the curve starts at (0, 0)
        tps = np.r_[0, tps]
        fps = np.r_[0, fps]
        thresholds = np.r_[thresholds[0] + 1, thresholds]

    if fps[-1] <= 0:
        warnings.warn("No negative samples in y_true, "
                      "false positive value should be meaningless",
                      UndefinedMetricWarning)
        fpr = np.repeat(np.nan, fps.shape)
    else:
        fpr = fps / fps[-1]

    if tps[-1] <= 0:
        warnings.warn("No positive samples in y_true, "
                      "true positive value should be meaningless",
                      UndefinedMetricWarning)
        tpr = np.repeat(np.nan, tps.shape)
    else:
        tpr = tps / tps[-1]

    return fpr, tpr, thresholds


def label_ranking_average_precision_score(y_true, y_score, sample_weight=None):
    """Compute ranking-based average precision

    Label ranking average precision (LRAP) is the average over each ground
    truth label assigned to each sample, of the ratio of true vs. total
    labels with lower score.

    This metric is used in multilabel ranking problem, where the goal
    is to give better rank to the labels associated to each sample.

    The obtained score is always strictly greater than 0 and
    the best value is 1.

    Read more in the :ref:`User Guide <label_ranking_average_precision>`.

    Parameters
    ----------
    y_true : array or sparse matrix, shape = [n_samples, n_labels]
        True binary labels in binary indicator format.

    y_score : array, shape = [n_samples, n_labels]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers).

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    score : float

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import label_ranking_average_precision_score
    >>> y_true = np.array([[1, 0, 0], [0, 0, 1]])
    >>> y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
    >>> label_ranking_average_precision_score(y_true, y_score) \
        # doctest: +ELLIPSIS
    0.416...

    """
    check_consistent_length(y_true, y_score, sample_weight)
    y_true = check_array(y_true, ensure_2d=False)
    y_score = check_array(y_score, ensure_2d=False)

    if y_true.shape != y_score.shape:
        raise ValueError("y_true and y_score have different shape")

    # Handle badly formatted array and the degenerate case with one label
    y_type = type_of_target(y_true)
    if (y_type != "multilabel-indicator" and
            not (y_type == "binary" and y_true.ndim == 2)):
        raise ValueError("{0} format is not supported".format(y_type))

    y_true = csr_matrix(y_true)
    y_score = -y_score

    n_samples, n_labels = y_true.shape

    out = 0.
    for i, (start, stop) in enumerate(zip(y_true.indptr, y_true.indptr[1:])):
        relevant = y_true.indices[start:stop]

        if (relevant.size == 0 or relevant.size == n_labels):
            # If all labels are relevant or unrelevant, the score is also
            # equal to 1. The label ranking has no meaning.
            out += 1.
            continue

        scores_i = y_score[i]
        rank = rankdata(scores_i, 'max')[relevant]
        L = rankdata(scores_i[relevant], 'max')
        aux = (L / rank).mean()
        if sample_weight is not None:
            aux = aux * sample_weight[i]
        out += aux

    if sample_weight is None:
        out /= n_samples
    else:
        out /= np.sum(sample_weight)

    return out


def coverage_error(y_true, y_score, sample_weight=None):
    """Coverage error measure

    Compute how far we need to go through the ranked scores to cover all
    true labels. The best value is equal to the average number
    of labels in ``y_true`` per sample.

    Ties in ``y_scores`` are broken by giving maximal rank that would have
    been assigned to all tied values.

    Note: Our implementation's score is 1 greater than the one given in
    Tsoumakas et al., 2010. This extends it to handle the degenerate case
    in which an instance has 0 true labels.

    Read more in the :ref:`User Guide <coverage_error>`.

    Parameters
    ----------
    y_true : array, shape = [n_samples, n_labels]
        True binary labels in binary indicator format.

    y_score : array, shape = [n_samples, n_labels]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers).

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    coverage_error : float

    References
    ----------
    .. [1] Tsoumakas, G., Katakis, I., & Vlahavas, I. (2010).
           Mining multi-label data. In Data mining and knowledge discovery
           handbook (pp. 667-685). Springer US.

    """
    y_true = check_array(y_true, ensure_2d=False)
    y_score = check_array(y_score, ensure_2d=False)
    check_consistent_length(y_true, y_score, sample_weight)

    y_type = type_of_target(y_true)
    if y_type != "multilabel-indicator":
        raise ValueError("{0} format is not supported".format(y_type))

    if y_true.shape != y_score.shape:
        raise ValueError("y_true and y_score have different shape")

    y_score_mask = np.ma.masked_array(y_score, mask=np.logical_not(y_true))
    y_min_relevant = y_score_mask.min(axis=1).reshape((-1, 1))
    coverage = (y_score >= y_min_relevant).sum(axis=1)
    coverage = coverage.filled(0)

    return np.average(coverage, weights=sample_weight)


def label_ranking_loss(y_true, y_score, sample_weight=None):
    """Compute Ranking loss measure

    Compute the average number of label pairs that are incorrectly ordered
    given y_score weighted by the size of the label set and the number of
    labels not in the label set.

    This is similar to the error set size, but weighted by the number of
    relevant and irrelevant labels. The best performance is achieved with
    a ranking loss of zero.

    Read more in the :ref:`User Guide <label_ranking_loss>`.

    .. versionadded:: 0.17
       A function *label_ranking_loss*

    Parameters
    ----------
    y_true : array or sparse matrix, shape = [n_samples, n_labels]
        True binary labels in binary indicator format.

    y_score : array, shape = [n_samples, n_labels]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or non-thresholded measure of decisions
        (as returned by "decision_function" on some classifiers).

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    loss : float

    References
    ----------
    .. [1] Tsoumakas, G., Katakis, I., & Vlahavas, I. (2010).
           Mining multi-label data. In Data mining and knowledge discovery
           handbook (pp. 667-685). Springer US.

    """
    y_true = check_array(y_true, ensure_2d=False, accept_sparse='csr')
    y_score = check_array(y_score, ensure_2d=False)
    check_consistent_length(y_true, y_score, sample_weight)

    y_type = type_of_target(y_true)
    if y_type not in ("multilabel-indicator",):
        raise ValueError("{0} format is not supported".format(y_type))

    if y_true.shape != y_score.shape:
        raise ValueError("y_true and y_score have different shape")

    n_samples, n_labels = y_true.shape

    y_true = csr_matrix(y_true)

    loss = np.zeros(n_samples)
    for i, (start, stop) in enumerate(zip(y_true.indptr, y_true.indptr[1:])):
        # Sort and bin the label scores
        unique_scores, unique_inverse = np.unique(y_score[i],
                                                  return_inverse=True)
        true_at_reversed_rank = np.bincount(
            unique_inverse[y_true.indices[start:stop]],
            minlength=len(unique_scores))
        all_at_reversed_rank = np.bincount(unique_inverse,
                                        minlength=len(unique_scores))
        false_at_reversed_rank = all_at_reversed_rank - true_at_reversed_rank

        # if the scores are ordered, it's possible to count the number of
        # incorrectly ordered paires in linear time by cumulatively counting
        # how many false labels of a given score have a score higher than the
        # accumulated true labels with lower score.
        loss[i] = np.dot(true_at_reversed_rank.cumsum(),
                         false_at_reversed_rank)

    n_positives = count_nonzero(y_true, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        loss /= ((n_labels - n_positives) * n_positives)

    # When there is no positive or no negative labels, those values should
    # be consider as correct, i.e. the ranking doesn't matter.
    loss[np.logical_or(n_positives == 0, n_positives == n_labels)] = 0.

    return np.average(loss, weights=sample_weight)


def amoc_curve(time, y_score, ts_ids, k_fold_label=None):
    """ Compute Activity Monitoring Operating Characteristic (AMOC) [1]

    The Activity Monitoring Operating Characteristic (AMOC) plots "false-alarms", measured in
    the false-positive-rate (FPR) on the y-axis versus the time to detection (defined in a moment) on the x-axis.
    The *time to detection* is defined as the first point in time after an event, by assumption at time=0, that the
    classifier predicts a positive label.

    The AMOC curve can be used in time series, binary classification settings where a label changes from negative to positive
    at a certain time, i.e. the label is negative before and positive after the event. The curve is of particular interest
    in settings where false alarms are costly, such as in medical applications.

    The intuition of the curve is the trade-off between the two metrics plotted: If the threshold of predicting a positive
    label is lowered, positive predictions are more likely, so that the time to detection is small, but the false-alarm
    rate ("wrong positive predictions") is increases.

    The structure of the input to this function is assumed to be the following: One must provide 4 arrays:
    1) time, the time at which the label is predicted. The event at which the label switches from 0 to 1 must be at time=0,
        i.e. the inputs to the function must be normalized accordingly beforehand.
    2) y_score, the prediction of a classifier for that label
    3) ts_ids, the ID of the time series. Every time series is constituted by multiple observations and their corresponding
        prediction and label. Example: In a medical setting, one patient might be observed over 3 hours, once per hour.
        Then, patient A has e.g. time series ID 0, where the part of the input for this patient looks like this:
        ts_ids = [..., 0, 0, 0, ...]
        y_true = [..., 0, 1, 1, ...] (inferred from input: time)
        y_pred = [..., 0.1, 0.7, 0.9, ...]
    Note that value i in any of the 4 lists corresponds to value i in any of the other 3 lists.

    Parameters
    ---------

    time: array, shape = [n_samples]
        Time stamp of the prediction.
        The event at which the inferred label changes from 0 (negative) to 1 (positive) must for all time series
        be at time = 0.0.
        The time stamps are assumed to be uniformly distributed.
        time can be given in any unit (e.g. in minutes).

    y_score: array, shape = [n_samples]
        Target scores, can, but don't have to be probabilities (i.e. in range [0, 1]).
        The higher a y_score value, the more likely the positive class (as of the prediction).

    ts_ids: array, shape = [n_samples]
        Time series ID corresponding to the time and y_score value.
        Each time series corresponds to observations and labels over time.
        Must be unique for every time series, also across k_folds.

    k_fold_label: array, shape = [n_samples], default=None
        The k-th fold the (test) time series corresponds to.
        Allows to compute confidence bounds over the folds.

    Returns
    -------

    time_to_detection: array
        Time of the first detection after the event occurs, averaged over all time series.

    FPR_mean: array
        False positive rate averaged over all time series

    FPR_std: array
        Standard deviation of false positive rate over the k folds.

    thresholds: array
        Thresholds used to compute the AMOC curve (cut-offs to decide on positive
        or negative prediction).

    Examples
    --------

    -> medical setting
    should actually work!!

    Potential improvements
    ----------------------

    * The Generalized-AMOC (G-AMOC) extends the standard AMOC by a protocol function that allows
        a more complex definition of a time to detection (details can be found in [1]). However, this
        implementation is limited to the standard AMOC implementation.
    * Make small toy example
    * In the current version, thresholds are drawn from a uniform distribution between the minimum
        and maximum prediction score. This can cause issues, if the predictions scores are not
        uniformly distributed.

    References
    ----------

    .. Implementation is adapted version as first published in:
    .. [1] `Fawcett, Tom, and Foster Provost. Activity monitoring: Noticing interesting changes in behavior.
           Proceedings of the fifth ACM SIGKDD international conference on Knowledge discovery and data mining.
           ACM, 1999. <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3654&rep=rep1&type=pdf>`_
    .. [2] `Jiang, Xia, Gregory F. Cooper, and Daniel B. Neill.
            Generalized AMOC curves for evaluation and improvement of event surveillance.
            AMIA Annual Symposium Proceedings. Vol. 2009. American Medical Informatics Association, 2009.`_

    Author
    ------

    Fabian Falck

    """

    # preliminary checks of inputs
    if k_fold_label is None:
        check_consistent_length(time, y_score, ts_ids)
        time, y_score, ts_ids = column_or_1d(time), column_or_1d(y_score), column_or_1d(ts_ids)
    else:
        check_consistent_length(time, y_score, ts_ids, k_fold_label)
        time, y_score, ts_ids, k_fold_label = column_or_1d(time), column_or_1d(y_score), column_or_1d(ts_ids), column_or_1d(k_fold_label)
        assert_all_finite(k_fold_label)
    assert_all_finite(time)
    assert_all_finite(y_score)
    assert_all_finite(ts_ids)

    # infer y_true from time stamps
    y_true = np.ones(y_score.shape[0])
    y_true[time < 0.0] = 0.0

    # if no k_folds given: only one fold infered
    if k_fold_label is None:
        k_fold_label = np.ones(y_score.shape[0])

    if k_fold_label is None:
        # pseudo k_fold_label with all 0s
        k_fold_label = np.zeros(y_true.shape[0])
    # find the minimum and maximum prediction score for generating thresholds later
    min_treshold, max_threshold = y_score.min()-0.01, y_score.max()+0.01
    # fold labels
    unique_folds = np.unique(k_fold_label)
    det_time_folds, FPR_folds = [], []
    # number of thresholds computed (the higher, the smoother the AMOC, but the longer the computation)
    resolution = 200

    for k, k_fold in enumerate(unique_folds):
        # labels, predictions and ts_ids corresponding to one fold
        y_true_kth, y_score_kth, time_kth, ts_ids_kth = y_true[k_fold_label == k_fold], y_score[k_fold_label == k_fold], time[k_fold_label == k_fold], ts_ids[k_fold_label == k_fold]
        # ts ids within one fold
        unique_ts_ids = np.unique(ts_ids_kth)

        # dictionary to store data filtered by time series ID
        data_by_ts = {}
        # number of negative observations
        N = 0.
        for i, ts_id in enumerate(unique_ts_ids):
            # labels and predictions corresponding to one time series
            y_true_ts = y_true_kth[ts_ids_kth == ts_id]
            y_score_ts = y_score_kth[ts_ids_kth == ts_id]
            time_ts = time_kth[ts_ids_kth == ts_id]
            # order them according to time
            order = np.argsort(time_ts)
            y_true_ts = y_true_ts[order]
            y_score_ts = y_score_ts[order]
            time_ts = time_ts[order]
            # store them in a per-time series dictionary
            data_by_ts[ts_id] = (y_true_ts, y_score_ts, time_ts)
            # add the negative observations of that time series
            N = N + float(np.count_nonzero(y_true_ts == 0.0))

        # storing the data points for each threshold
        det_time_thres, FPR_thres = [], []
        thresholds = np.linspace(min_treshold, max_threshold, num=resolution)
        for j, threshold in enumerate(thresholds):
            # counting the false positives over all time series
            FP = 0.0
            # storing time to detection (earliest time after the event with a positive prediction)
            # for each time series
            det_time = []
            for i, ts_id in enumerate(unique_ts_ids):
                # whether the event was detected at all or not
                detected = False
                y_true_ts, y_score_ts, time_ts = data_by_ts[ts_id]

                # loop over the time series
                for j, t in enumerate(time_ts):
                    if y_score_ts[j] > threshold and t < 0:  # FP
                        FP += 1.
                    elif y_score_ts[j] > threshold and t >= 0:  # TP
                        detected = True
                        det_time.append(t)
                        break
                # event was not detected at all
                if not detected:
                    # assumption: adding the maximum observation time
                    det_time.append(time_ts.max())
            det_time_thres.append(np.mean(det_time))
            FPR_thres.append(FP / N)
        det_time_folds.append(np.array(det_time_thres))
        FPR_folds.append(np.array(FPR_thres))

    # find the maximum time to detection
    det_time_max = 0.0
    for det_time_thres in det_time_folds:
        det_time_max = max(det_time_thres.max(), det_time_max)  # maximum of already found max and current max

    #interpolate the k fold curves
    time_to_detection = np.linspace(0, det_time_max, resolution)  #
    FPR_interps = []
    for det_time_thres, FPR_thres in zip(det_time_folds, FPR_folds):
        # sort by det_time
        order = np.argsort(det_time_thres)
        det_time_thres = det_time_thres[order]
        FPR_thres = FPR_thres[order]
        FPR_interp = np.interp(time_to_detection, det_time_thres, FPR_thres)
        FPR_interps.append(FPR_interp)
    # compute mean and standard deviation of FPR
    FPR_interps = np.array(FPR_interps)
    FPR_mean = np.mean(FPR_interps, axis=0)
    FPR_std = np.std(FPR_interps, axis=0)

    return time_to_detection, FPR_mean, FPR_std, thresholds


def plot_amoc(time_list, y_score_list, ts_id_list, k_fold_label_list, model_labels, evaluation_path, conf_stds=1.0):
    """
    Example of how to use the amoc_curve() function with multiple test sets to compare.

    :param time_list: list of time arrays, where one element is an input to amoc_curve()
    :param y_score_list: list of prediction arrays, where one element is an input to amoc_curve()
    :param ts_id_list: list of time series ID arrays, where one element is an input to amoc_curve()
    :param k_fold_label_list: list of k_fold label arrays, where one element is an input to amoc_curve()
    :param model_labels: list of strings. Describing the name of the model, displayed in the legend of the figure.
    :param evaluation_path: a string of the path where to save the AMOC figures to.
    :param conf_stds: controlling how thick the confidence bounds are in temrs of multiples of the standard deviation.
            Default: 1 standard deviation
    """

    number_of_curves = len(y_score_list)
    time_to_detection_list, FPR_mean_list, FPR_std_list, thresholds_list = [], [], [], []

    # compute the amoc curves for mutiple models that are compared
    for c in range(number_of_curves):
        # /60. in order to scale to minutes (original data was given in seconds)
        time_to_detection, FPR_mean, FPR_std, thresholds = amoc_curve(time_list[c], y_score_list[c], ts_id_list[c], k_fold_label_list[c])
        time_to_detection_list.append(time_to_detection)
        FPR_mean_list.append(FPR_mean)
        FPR_std_list.append(FPR_std)
        thresholds_list.append(thresholds)

    # plot the amoc curves all in one figure
    for plot_nr in range(3):  # plot both unit-unit and unit-log version
        fig = plt.figure(figsize=(8.0, 5.0))
        ax = fig.add_subplot(111)
        colors = ['blue', 'orange', 'green', 'purple', 'grey']
        x_axis_max = 5  # upper bound of x axis
        global_min = 1.0
        for c in range(number_of_curves):
            time_to_detection, FPR_mean, FPR_std, thresholds = time_to_detection_list[c], FPR_mean_list[c], FPR_std_list[c], thresholds_list[c]
            plt.plot(time_to_detection, FPR_mean, color=colors[c], label=model_labels[c])
            # compute lower and upper confidence bounds with one standard deviation
            lower_conf = np.maximum(FPR_mean - conf_stds * FPR_std, 0)  # lower bound is 0
            upper_conf = np.minimum(FPR_mean + conf_stds * FPR_std, 1)  # upper bound is 1
            # plot the confidence bounds
            plt.fill_between(time_to_detection, lower_conf, upper_conf, color=colors[c], alpha=0.15)
            # find global minimum in plot for y axis in log plot
            lower_conf_plotted = lower_conf[time_to_detection <= x_axis_max]
            global_min = min(lower_conf_plotted.min(), global_min)
        # cosmetics
        plt.title('Activity Monitoring Operating Characteristic (AMOC)')
        plt.xlabel('Time to detection (minutes)', fontsize='large')
        plt.xlim(xmin=0, xmax=x_axis_max)
        plt.legend(loc='upper right', prop={'size': 9})

        # reg scale
        if plot_nr == 1:
            plt.ylim(ymin=0, ymax=1)
            plt.ylabel(r'$FPR = \frac{FP_{obs}}{N_{obs}}$', fontsize='large')
            plt.savefig(os.path.join(evaluation_path, 'AMOC_det-FPR_u-u' + '.pdf'), format='pdf', dpi=1200)
        # log scale
        elif plot_nr == 2:
            # cut log scale off at clostest power of 10 smaller than the global minimum of "anything plotted"
            # log must be defined
            if not global_min <= 0.0:
                closest_10base = 10 ** floor(log(global_min, 10))
                plt.ylim(ymin=closest_10base, ymax=1)
            plt.ylabel(r'$FPR = \frac{FP_{obs}}{N_{obs}}$ (log scale)', fontsize='large')
            ax.set_yscale('log')
            plt.savefig(os.path.join(evaluation_path, 'AMOC_det-FPR_u-log' + '.pdf'), format='pdf', dpi=1200)

