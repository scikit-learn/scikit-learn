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
# License: BSD 3 clause


import warnings
from functools import partial

import numpy as np
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
from ..preprocessing._label import _encode

from ._base import _average_binary_score, _average_multiclass_ovo_score


def auc(x, y):
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

    direction = 1
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


def average_precision_score(y_true, y_score, average="macro", pos_label=1,
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
        True binary labels or binary label indicators.

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

        Will be ignored when ``y_true`` is binary.

    pos_label : int or str (default=1)
        The label of the positive class. Only applied to binary ``y_true``.
        For multilabel-indicator ``y_true``, ``pos_label`` is fixed to 1.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    Returns
    -------
    average_precision : float

    References
    ----------
    .. [1] `Wikipedia entry for the Average precision
           <https://en.wikipedia.org/w/index.php?title=Information_retrieval&
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
    >>> average_precision_score(y_true, y_scores)
    0.83...

    Notes
    -----
    .. versionchanged:: 0.19
      Instead of linearly interpolating between operating points, precisions
      are weighted by the change in recall since the last operating point.
    """
    def _binary_uninterpolated_average_precision(
            y_true, y_score, pos_label=1, sample_weight=None):
        precision, recall, _ = precision_recall_curve(
            y_true, y_score, pos_label=pos_label, sample_weight=sample_weight)
        # Return the step function integral
        # The following works because the last entry of precision is
        # guaranteed to be 1, as returned by precision_recall_curve
        return -np.sum(np.diff(recall) * np.array(precision)[:-1])

    y_type = type_of_target(y_true)
    if y_type == "multilabel-indicator" and pos_label != 1:
        raise ValueError("Parameter pos_label is fixed to 1 for "
                         "multilabel-indicator y_true. Do not set "
                         "pos_label or set pos_label to 1.")
    elif y_type == "binary":
        present_labels = np.unique(y_true)
        if len(present_labels) == 2 and pos_label not in present_labels:
            raise ValueError("pos_label=%r is invalid. Set it to a label in "
                             "y_true." % pos_label)
    average_precision = partial(_binary_uninterpolated_average_precision,
                                pos_label=pos_label)
    return _average_binary_score(average_precision, y_true, y_score,
                                 average, sample_weight=sample_weight)


def _binary_roc_auc_score(y_true, y_score, sample_weight=None, max_fpr=None):
    """Binary roc auc score"""
    if len(np.unique(y_true)) != 2:
        raise ValueError("Only one class present in y_true. ROC AUC score "
                         "is not defined in that case.")

    fpr, tpr, _ = roc_curve(y_true, y_score,
                            sample_weight=sample_weight)
    if max_fpr is None or max_fpr == 1:
        return auc(fpr, tpr)
    if max_fpr <= 0 or max_fpr > 1:
        raise ValueError("Expected max_fpr in range (0, 1], got: %r" % max_fpr)

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


def roc_auc_score(y_true, y_score, average="macro", sample_weight=None,
                  max_fpr=None, multi_class="raise", labels=None):
    """Compute Area Under the Receiver Operating Characteristic Curve (ROC AUC)
    from prediction scores.

    Note: this implementation can be used with binary, multiclass and
    multilabel classification, but some restrictions apply (see Parameters).

    Read more in the :ref:`User Guide <roc_metrics>`.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,) or (n_samples, n_classes)
        True labels or binary label indicators. The binary and multiclass cases
        expect labels with shape (n_samples,) while the multilabel case expects
        binary label indicators with shape (n_samples, n_classes).

    y_score : array-like of shape (n_samples,) or (n_samples, n_classes)
        Target scores. In the binary and multilabel cases, these can be either
        probability estimates or non-thresholded decision values (as returned
        by `decision_function` on some classifiers). In the multiclass case,
        these must be probability estimates which sum to 1. The binary
        case expects a shape (n_samples,), and the scores must be the scores of
        the class with the greater label. The multiclass and multilabel
        cases expect a shape (n_samples, n_classes). In the multiclass case,
        the order of the class scores must correspond to the order of
        ``labels``, if provided, or else to the numerical or lexicographical
        order of the labels in ``y_true``.

    average : {'micro', 'macro', 'samples', 'weighted'} or None, \
            default='macro'
        If ``None``, the scores for each class are returned. Otherwise,
        this determines the type of averaging performed on the data:
        Note: multiclass ROC AUC currently only handles the 'macro' and
        'weighted' averages.

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

        Will be ignored when ``y_true`` is binary.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    max_fpr : float > 0 and <= 1, default=None
        If not ``None``, the standardized partial AUC [2]_ over the range
        [0, max_fpr] is returned. For the multiclass case, ``max_fpr``,
        should be either equal to ``None`` or ``1.0`` as AUC ROC partial
        computation currently is not supported for multiclass.

    multi_class : {'raise', 'ovr', 'ovo'}, default='raise'
        Multiclass only. Determines the type of configuration to use. The
        default value raises an error, so either ``'ovr'`` or ``'ovo'`` must be
        passed explicitly.

        ``'ovr'``:
            Computes the AUC of each class against the rest [3]_ [4]_. This
            treats the multiclass case in the same way as the multilabel case.
            Sensitive to class imbalance even when ``average == 'macro'``,
            because class imbalance affects the composition of each of the
            'rest' groupings.
        ``'ovo'``:
            Computes the average AUC of all possible pairwise combinations of
            classes [5]_. Insensitive to class imbalance when
            ``average == 'macro'``.

    labels : array-like of shape (n_classes,), default=None
        Multiclass only. List of labels that index the classes in ``y_score``.
        If ``None``, the numerical or lexicographical order of the labels in
        ``y_true`` is used.

    Returns
    -------
    auc : float

    References
    ----------
    .. [1] `Wikipedia entry for the Receiver operating characteristic
            <https://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_

    .. [2] `Analyzing a portion of the ROC curve. McClish, 1989
            <https://www.ncbi.nlm.nih.gov/pubmed/2668680>`_

    .. [3] Provost, F., Domingos, P. (2000). Well-trained PETs: Improving
           probability estimation trees (Section 6.2), CeDER Working Paper
           #IS-00-04, Stern School of Business, New York University.

    .. [4] `Fawcett, T. (2006). An introduction to ROC analysis. Pattern
            Recognition Letters, 27(8), 861-874.
            <https://www.sciencedirect.com/science/article/pii/S016786550500303X>`_

    .. [5] `Hand, D.J., Till, R.J. (2001). A Simple Generalisation of the Area
            Under the ROC Curve for Multiple Class Classification Problems.
            Machine Learning, 45(2), 171-186.
            <http://link.springer.com/article/10.1023/A:1010920819831>`_

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

    y_type = type_of_target(y_true)
    y_true = check_array(y_true, ensure_2d=False, dtype=None)
    y_score = check_array(y_score, ensure_2d=False)

    if y_type == "multiclass" or (y_type == "binary" and
                                  y_score.ndim == 2 and
                                  y_score.shape[1] > 2):
        # do not support partial ROC computation for multiclass
        if max_fpr is not None and max_fpr != 1.:
            raise ValueError("Partial AUC computation not available in "
                             "multiclass setting, 'max_fpr' must be"
                             " set to `None`, received `max_fpr={0}` "
                             "instead".format(max_fpr))
        if multi_class == 'raise':
            raise ValueError("multi_class must be in ('ovo', 'ovr')")
        return _multiclass_roc_auc_score(y_true, y_score, labels,
                                         multi_class, average, sample_weight)
    elif y_type == "binary":
        labels = np.unique(y_true)
        y_true = label_binarize(y_true, labels)[:, 0]
        return _average_binary_score(partial(_binary_roc_auc_score,
                                             max_fpr=max_fpr),
                                     y_true, y_score, average,
                                     sample_weight=sample_weight)
    else:  # multilabel-indicator
        return _average_binary_score(partial(_binary_roc_auc_score,
                                             max_fpr=max_fpr),
                                     y_true, y_score, average,
                                     sample_weight=sample_weight)


def _multiclass_roc_auc_score(y_true, y_score, labels,
                              multi_class, average, sample_weight):
    """Multiclass roc auc score

    Parameters
    ----------
    y_true : array-like of shape (n_samples,)
        True multiclass labels.

    y_score : array-like of shape (n_samples, n_classes)
        Target scores corresponding to probability estimates of a sample
        belonging to a particular class

    labels : array, shape = [n_classes] or None, optional (default=None)
        List of labels to index ``y_score`` used for multiclass. If ``None``,
        the lexical order of ``y_true`` is used to index ``y_score``.

    multi_class : string, 'ovr' or 'ovo'
        Determines the type of multiclass configuration to use.
        ``'ovr'``:
            Calculate metrics for the multiclass case using the one-vs-rest
            approach.
        ``'ovo'``:
            Calculate metrics for the multiclass case using the one-vs-one
            approach.

    average : 'macro' or 'weighted', optional (default='macro')
        Determines the type of averaging performed on the pairwise binary
        metric scores
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean. This does not take label imbalance into account. Classes
            are assumed to be uniformly distributed.
        ``'weighted'``:
            Calculate metrics for each label, taking into account the
            prevalence of the classes.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.

    """
    # validation of the input y_score
    if not np.allclose(1, y_score.sum(axis=1)):
        raise ValueError(
            "Target scores need to be probabilities for multiclass "
            "roc_auc, i.e. they should sum up to 1.0 over classes")

    # validation for multiclass parameter specifications
    average_options = ("macro", "weighted")
    if average not in average_options:
        raise ValueError("average must be one of {0} for "
                         "multiclass problems".format(average_options))

    multiclass_options = ("ovo", "ovr")
    if multi_class not in multiclass_options:
        raise ValueError("multi_class='{0}' is not supported "
                         "for multiclass ROC AUC, multi_class must be "
                         "in {1}".format(
                                multi_class, multiclass_options))

    if labels is not None:
        labels = column_or_1d(labels)
        classes = _encode(labels)
        if len(classes) != len(labels):
            raise ValueError("Parameter 'labels' must be unique")
        if not np.array_equal(classes, labels):
            raise ValueError("Parameter 'labels' must be ordered")
        if len(classes) != y_score.shape[1]:
            raise ValueError(
                "Number of given labels, {0}, not equal to the number "
                "of columns in 'y_score', {1}".format(
                    len(classes), y_score.shape[1]))
        if len(np.setdiff1d(y_true, classes)):
            raise ValueError(
                "'y_true' contains labels not in parameter 'labels'")
    else:
        classes = _encode(y_true)
        if len(classes) != y_score.shape[1]:
            raise ValueError(
                "Number of classes in y_true not equal to the number of "
                "columns in 'y_score'")

    if multi_class == "ovo":
        if sample_weight is not None:
            raise ValueError("sample_weight is not supported "
                             "for multiclass one-vs-one ROC AUC, "
                             "'sample_weight' must be None in this case.")
        _, y_true_encoded = _encode(y_true, uniques=classes, encode=True)
        # Hand & Till (2001) implementation (ovo)
        return _average_multiclass_ovo_score(_binary_roc_auc_score,
                                             y_true_encoded,
                                             y_score, average=average)
    else:
        # ovr is same as multi-label
        y_true_multilabel = label_binarize(y_true, classes)
        return _average_binary_score(_binary_roc_auc_score, y_true_multilabel,
                                     y_score, average,
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

    sample_weight : array-like of shape (n_samples,), default=None
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
    # classes.dtype.kind in ('O', 'U', 'S') is required to avoid
    # triggering a FutureWarning by calling np.array_equal(a, b)
    # when elements in the two arrays are not comparable.
    classes = np.unique(y_true)
    if (pos_label is None and (
            classes.dtype.kind in ('O', 'U', 'S') or
            not (np.array_equal(classes, [0, 1]) or
                 np.array_equal(classes, [-1, 1]) or
                 np.array_equal(classes, [0]) or
                 np.array_equal(classes, [-1]) or
                 np.array_equal(classes, [1])))):
        classes_repr = ", ".join(repr(c) for c in classes)
        raise ValueError("y_true takes value in {{{classes_repr}}} and "
                         "pos_label is not specified: either make y_true "
                         "take value in {{0, 1}} or {{-1, 1}} or "
                         "pass pos_label explicitly.".format(
                             classes_repr=classes_repr))
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
        True binary labels. If labels are not either {-1, 1} or {0, 1}, then
        pos_label should be explicitly given.

    probas_pred : array, shape = [n_samples]
        Estimated probabilities or decision function.

    pos_label : int or str, default=None
        The label of the positive class.
        When ``pos_label=None``, if y_true is in {-1, 1} or {0, 1},
        ``pos_label`` is set to 1, otherwise an error will be raised.

    sample_weight : array-like of shape (n_samples,), default=None
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
    >>> precision
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
    precision[np.isnan(precision)] = 0
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
        The label of the positive class.
        When ``pos_label=None``, if y_true is in {-1, 1} or {0, 1},
        ``pos_label`` is set to 1, otherwise an error will be raised.

    sample_weight : array-like of shape (n_samples,), default=None
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

    # Add an extra threshold position
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

    sample_weight : array-like of shape (n_samples,), default=None
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
    >>> label_ranking_average_precision_score(y_true, y_score)
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
            aux = 1.
        else:
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

    sample_weight : array-like of shape (n_samples,), default=None
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

    sample_weight : array-like of shape (n_samples,), default=None
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


def _dcg_sample_scores(y_true, y_score, k=None,
                       log_base=2, ignore_ties=False):
    """Compute Discounted Cumulative Gain.

    Sum the true scores ranked in the order induced by the predicted scores,
    after applying a logarithmic discount.

    This ranking metric yields a high value if true labels are ranked high by
    ``y_score``.

    Parameters
    ----------
    y_true : ndarray, shape (n_samples, n_labels)
        True targets of multilabel classification, or true scores of entities
        to be ranked.

    y_score : ndarray, shape (n_samples, n_labels)
        Target scores, can either be probability estimates, confidence values,
        or non-thresholded measure of decisions (as returned by
        "decision_function" on some classifiers).

    k : int, optional (default=None)
        Only consider the highest k scores in the ranking. If None, use all
        outputs.

    log_base : float, optional (default=2)
        Base of the logarithm used for the discount. A low value means a
        sharper discount (top results are more important).

    ignore_ties : bool, optional (default=False)
        Assume that there are no ties in y_score (which is likely to be the
        case if y_score is continuous) for efficiency gains.

    Returns
    -------
    discounted_cumulative_gain : ndarray, shape (n_samples,)
        The DCG score for each sample.

    See also
    --------
    ndcg_score :
        The Discounted Cumulative Gain divided by the Ideal Discounted
        Cumulative Gain (the DCG obtained for a perfect ranking), in order to
        have a score between 0 and 1.

    """
    discount = 1 / (np.log(np.arange(y_true.shape[1]) + 2) / np.log(log_base))
    if k is not None:
        discount[k:] = 0
    if ignore_ties:
        ranking = np.argsort(y_score)[:, ::-1]
        ranked = y_true[np.arange(ranking.shape[0])[:, np.newaxis], ranking]
        cumulative_gains = discount.dot(ranked.T)
    else:
        discount_cumsum = np.cumsum(discount)
        cumulative_gains = [_tie_averaged_dcg(y_t, y_s, discount_cumsum)
                            for y_t, y_s in zip(y_true, y_score)]
        cumulative_gains = np.asarray(cumulative_gains)
    return cumulative_gains


def _tie_averaged_dcg(y_true, y_score, discount_cumsum):
    """
    Compute DCG by averaging over possible permutations of ties.

    The gain (`y_true`) of an index falling inside a tied group (in the order
    induced by `y_score`) is replaced by the average gain within this group.
    The discounted gain for a tied group is then the average `y_true` within
    this group times the sum of discounts of the corresponding ranks.

    This amounts to averaging scores for all possible orderings of the tied
    groups.

    (note in the case of dcg@k the discount is 0 after index k)

    Parameters
    ----------
    y_true : ndarray
        The true relevance scores

    y_score : ndarray
        Predicted scores

    discount_cumsum : ndarray
        Precomputed cumulative sum of the discounts.

    Returns
    -------
    The discounted cumulative gain.

    References
    ----------
    McSherry, F., & Najork, M. (2008, March). Computing information retrieval
    performance measures efficiently in the presence of tied scores. In
    European conference on information retrieval (pp. 414-421). Springer,
    Berlin, Heidelberg.

    """
    _, inv, counts = np.unique(
        - y_score, return_inverse=True, return_counts=True)
    ranked = np.zeros(len(counts))
    np.add.at(ranked, inv, y_true)
    ranked /= counts
    groups = np.cumsum(counts) - 1
    discount_sums = np.empty(len(counts))
    discount_sums[0] = discount_cumsum[groups[0]]
    discount_sums[1:] = np.diff(discount_cumsum[groups])
    return (ranked * discount_sums).sum()


def _check_dcg_target_type(y_true):
    y_type = type_of_target(y_true)
    supported_fmt = ("multilabel-indicator", "continuous-multioutput",
                     "multiclass-multioutput")
    if y_type not in supported_fmt:
        raise ValueError(
            "Only {} formats are supported. Got {} instead".format(
                supported_fmt, y_type))


def dcg_score(y_true, y_score, k=None,
              log_base=2, sample_weight=None, ignore_ties=False):
    """Compute Discounted Cumulative Gain.

    Sum the true scores ranked in the order induced by the predicted scores,
    after applying a logarithmic discount.

    This ranking metric yields a high value if true labels are ranked high by
    ``y_score``.

    Usually the Normalized Discounted Cumulative Gain (NDCG, computed by
    ndcg_score) is preferred.

    Parameters
    ----------
    y_true : ndarray, shape (n_samples, n_labels)
        True targets of multilabel classification, or true scores of entities
        to be ranked.

    y_score : ndarray, shape (n_samples, n_labels)
        Target scores, can either be probability estimates, confidence values,
        or non-thresholded measure of decisions (as returned by
        "decision_function" on some classifiers).

    k : int, optional (default=None)
        Only consider the highest k scores in the ranking. If None, use all
        outputs.

    log_base : float, optional (default=2)
        Base of the logarithm used for the discount. A low value means a
        sharper discount (top results are more important).

    sample_weight : ndarray, shape (n_samples,), optional (default=None)
        Sample weights. If None, all samples are given the same weight.

    ignore_ties : bool, optional (default=False)
        Assume that there are no ties in y_score (which is likely to be the
        case if y_score is continuous) for efficiency gains.

    Returns
    -------
    discounted_cumulative_gain : float
        The averaged sample DCG scores.

    See also
    --------
    ndcg_score :
        The Discounted Cumulative Gain divided by the Ideal Discounted
        Cumulative Gain (the DCG obtained for a perfect ranking), in order to
        have a score between 0 and 1.

    References
    ----------
    `Wikipedia entry for Discounted Cumulative Gain
    <https://en.wikipedia.org/wiki/Discounted_cumulative_gain>`_

    Jarvelin, K., & Kekalainen, J. (2002).
    Cumulated gain-based evaluation of IR techniques. ACM Transactions on
    Information Systems (TOIS), 20(4), 422-446.

    Wang, Y., Wang, L., Li, Y., He, D., Chen, W., & Liu, T. Y. (2013, May).
    A theoretical analysis of NDCG ranking measures. In Proceedings of the 26th
    Annual Conference on Learning Theory (COLT 2013)

    McSherry, F., & Najork, M. (2008, March). Computing information retrieval
    performance measures efficiently in the presence of tied scores. In
    European conference on information retrieval (pp. 414-421). Springer,
    Berlin, Heidelberg.

    Examples
    --------
    >>> from sklearn.metrics import dcg_score
    >>> # we have groud-truth relevance of some answers to a query:
    >>> true_relevance = np.asarray([[10, 0, 0, 1, 5]])
    >>> # we predict scores for the answers
    >>> scores = np.asarray([[.1, .2, .3, 4, 70]])
    >>> dcg_score(true_relevance, scores) # doctest: +ELLIPSIS
    9.49...
    >>> # we can set k to truncate the sum; only top k answers contribute
    >>> dcg_score(true_relevance, scores, k=2) # doctest: +ELLIPSIS
    5.63...
    >>> # now we have some ties in our prediction
    >>> scores = np.asarray([[1, 0, 0, 0, 1]])
    >>> # by default ties are averaged, so here we get the average true
    >>> # relevance of our top predictions: (10 + 5) / 2 = 7.5
    >>> dcg_score(true_relevance, scores, k=1) # doctest: +ELLIPSIS
    7.5
    >>> # we can choose to ignore ties for faster results, but only
    >>> # if we know there aren't ties in our scores, otherwise we get
    >>> # wrong results:
    >>> dcg_score(true_relevance,
    ...           scores, k=1, ignore_ties=True) # doctest: +ELLIPSIS
    5.0

    """
    y_true = check_array(y_true, ensure_2d=False)
    y_score = check_array(y_score, ensure_2d=False)
    check_consistent_length(y_true, y_score, sample_weight)
    _check_dcg_target_type(y_true)
    return np.average(
        _dcg_sample_scores(
            y_true, y_score, k=k, log_base=log_base,
            ignore_ties=ignore_ties),
        weights=sample_weight)


def _ndcg_sample_scores(y_true, y_score, k=None, ignore_ties=False):
    """Compute Normalized Discounted Cumulative Gain.

    Sum the true scores ranked in the order induced by the predicted scores,
    after applying a logarithmic discount. Then divide by the best possible
    score (Ideal DCG, obtained for a perfect ranking) to obtain a score between
    0 and 1.

    This ranking metric yields a high value if true labels are ranked high by
    ``y_score``.

    Parameters
    ----------
    y_true : ndarray, shape (n_samples, n_labels)
        True targets of multilabel classification, or true scores of entities
        to be ranked.

    y_score : ndarray, shape (n_samples, n_labels)
        Target scores, can either be probability estimates, confidence values,
        or non-thresholded measure of decisions (as returned by
        "decision_function" on some classifiers).

    k : int, optional (default=None)
        Only consider the highest k scores in the ranking. If None, use all
        outputs.

    ignore_ties : bool, optional (default=False)
        Assume that there are no ties in y_score (which is likely to be the
        case if y_score is continuous) for efficiency gains.

    Returns
    -------
    normalized_discounted_cumulative_gain : ndarray, shape (n_samples,)
        The NDCG score for each sample (float in [0., 1.]).

    See also
    --------
    dcg_score : Discounted Cumulative Gain (not normalized).

    """
    gain = _dcg_sample_scores(y_true, y_score, k, ignore_ties=ignore_ties)
    # Here we use the order induced by y_true so we can ignore ties since
    # the gain associated to tied indices is the same (permuting ties doesn't
    # change the value of the re-ordered y_true)
    normalizing_gain = _dcg_sample_scores(y_true, y_true, k, ignore_ties=True)
    all_irrelevant = normalizing_gain == 0
    gain[all_irrelevant] = 0
    gain[~all_irrelevant] /= normalizing_gain[~all_irrelevant]
    return gain


def ndcg_score(y_true, y_score, k=None, sample_weight=None, ignore_ties=False):
    """Compute Normalized Discounted Cumulative Gain.

    Sum the true scores ranked in the order induced by the predicted scores,
    after applying a logarithmic discount. Then divide by the best possible
    score (Ideal DCG, obtained for a perfect ranking) to obtain a score between
    0 and 1.

    This ranking metric yields a high value if true labels are ranked high by
    ``y_score``.

    Parameters
    ----------
    y_true : ndarray, shape (n_samples, n_labels)
        True targets of multilabel classification, or true scores of entities
        to be ranked.

    y_score : ndarray, shape (n_samples, n_labels)
        Target scores, can either be probability estimates, confidence values,
        or non-thresholded measure of decisions (as returned by
        "decision_function" on some classifiers).

    k : int, optional (default=None)
        Only consider the highest k scores in the ranking. If None, use all
        outputs.

    sample_weight : ndarray, shape (n_samples,), optional (default=None)
        Sample weights. If None, all samples are given the same weight.

    ignore_ties : bool, optional (default=False)
        Assume that there are no ties in y_score (which is likely to be the
        case if y_score is continuous) for efficiency gains.

    Returns
    -------
    normalized_discounted_cumulative_gain : float in [0., 1.]
        The averaged NDCG scores for all samples.

    See also
    --------
    dcg_score : Discounted Cumulative Gain (not normalized).

    References
    ----------
    `Wikipedia entry for Discounted Cumulative Gain
    <https://en.wikipedia.org/wiki/Discounted_cumulative_gain>`_

    Jarvelin, K., & Kekalainen, J. (2002).
    Cumulated gain-based evaluation of IR techniques. ACM Transactions on
    Information Systems (TOIS), 20(4), 422-446.

    Wang, Y., Wang, L., Li, Y., He, D., Chen, W., & Liu, T. Y. (2013, May).
    A theoretical analysis of NDCG ranking measures. In Proceedings of the 26th
    Annual Conference on Learning Theory (COLT 2013)

    McSherry, F., & Najork, M. (2008, March). Computing information retrieval
    performance measures efficiently in the presence of tied scores. In
    European conference on information retrieval (pp. 414-421). Springer,
    Berlin, Heidelberg.

    Examples
    --------
    >>> from sklearn.metrics import ndcg_score
    >>> # we have groud-truth relevance of some answers to a query:
    >>> true_relevance = np.asarray([[10, 0, 0, 1, 5]])
    >>> # we predict some scores (relevance) for the answers
    >>> scores = np.asarray([[.1, .2, .3, 4, 70]])
    >>> ndcg_score(true_relevance, scores) # doctest: +ELLIPSIS
    0.69...
    >>> scores = np.asarray([[.05, 1.1, 1., .5, .0]])
    >>> ndcg_score(true_relevance, scores) # doctest: +ELLIPSIS
    0.49...
    >>> # we can set k to truncate the sum; only top k answers contribute.
    >>> ndcg_score(true_relevance, scores, k=4) # doctest: +ELLIPSIS
    0.35...
    >>> # the normalization takes k into account so a perfect answer
    >>> # would still get 1.0
    >>> ndcg_score(true_relevance, true_relevance, k=4) # doctest: +ELLIPSIS
    1.0
    >>> # now we have some ties in our prediction
    >>> scores = np.asarray([[1, 0, 0, 0, 1]])
    >>> # by default ties are averaged, so here we get the average (normalized)
    >>> # true relevance of our top predictions: (10 / 10 + 5 / 10) / 2 = .75
    >>> ndcg_score(true_relevance, scores, k=1) # doctest: +ELLIPSIS
    0.75
    >>> # we can choose to ignore ties for faster results, but only
    >>> # if we know there aren't ties in our scores, otherwise we get
    >>> # wrong results:
    >>> ndcg_score(true_relevance,
    ...           scores, k=1, ignore_ties=True) # doctest: +ELLIPSIS
    0.5

    """
    y_true = check_array(y_true, ensure_2d=False)
    y_score = check_array(y_score, ensure_2d=False)
    check_consistent_length(y_true, y_score, sample_weight)
    _check_dcg_target_type(y_true)
    gain = _ndcg_sample_scores(y_true, y_score, k=k, ignore_ties=ignore_ties)
    return np.average(gain, weights=sample_weight)


def top_k_accuracy_score(y_true, y_score, k=3, normalize=True):
    """Top k Accuracy multiclass classification score.

    This metric computes the number of times where the correct label is among
    the top `k` labels predicted (ranked by predicted scores). Note that
    multilabel classification case isn't covered here.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,)
        True labels. Expected be a sequence of int in the range
        (0, n_classes - 1).

    y_score : array-like of shape (n_samples, n_classes)
        Target scores. These can be either probability estimates or
        non-thresholded decision values as returned by :term:`predict_proba` or
        :term:`decision_function`, respectively.

    k : int, default=3
        Number of guesses allowed to find the correct label.

    normalize : bool, default=True
        If `True`, return the fraction of correctly classified samples.
        Otherwise, return the number of correctly classified samples.

    Returns
    -------
    score : float
        The top-k accuracy score. The best performance is 1 with
        `normalize == True` and the number of samples with
        `normalize == False`.

    See also
    --------
    accuracy_score

    Notes
    -----
    In cases where two or more labels are assigned equal predicted scores,
    the labels with the smallest indices will be chosen first. Then, the result
    might be incorrect if the correct label falls after the threshold because
    of that.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import top_k_accuracy_score
    >>> y_true = np.array([0, 1, 2, 2])
    >>> y_score = np.array([[0.5, 0.2, 0.2],
    ...                     [0.3, 0.4, 0.2],
    ...                     [0.2, 0.4, 0.3],
    ...                     [0.7, 0.2, 0.1]])
    >>> top_k_accuracy_score(y_true, y_score, k=2)
    0.75
    >>> top_k_accuracy_score(y_true, y_score, k=2, normalize=False)
    3
    >>> top_k_accuracy_score([0]*4, y_score, k=2)
    0.75

    """
    if k == 1:
        raise ValueError(
            "'k'=1 is equivalent to 'metrics.accuracy_score'. Please, use "
            "that method instead."
        )

    check_consistent_length(y_true, y_score)
    y_true_type = type_of_target(y_true)

    if y_true_type not in {'binary', 'multiclass'}:
        raise ValueError(
            f"'y_true' type must be 'binary' or 'multiclass', got "
            f"'{y_true_type}' instead."
        )

    y_true = column_or_1d(y_true)
    y_score = check_array(y_score, ensure_min_features=3)

    n_classes = y_score.shape[1]
    if k >= n_classes:
        raise ValueError(
            f"'k' ({k}) can't be greater than or equal to 'n_classes' "
            f"({n_classes}). Will result in a perfect score and is therefore "
            "meaningless."
        )

    infered_n_classes = y_true.max() + 1
    if infered_n_classes > n_classes:
        raise ValueError(
            f"Number of classes in 'y_true' {infered_n_classes} is greater "
            f"than the number of columns in 'y_score' {n_classes}."
        )

    sorted_pred = np.argsort(-y_score, axis=1)
    score = sum(y in pred[:k] for y, pred in zip(y_true, sorted_pred))
    score = score / len(y_true) if normalize else score
    return score
