"""Utilities to evaluate the predictive performance of models

Functions named as ``*_score`` return a scalar value to maximize: the higher
the better

Function named as ``*_error`` or ``*_loss`` return a scalar value to minimize:
the lower the better
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Arnaud Joly <a.joly@ulg.ac.be>
# License: BSD Style.

from itertools import izip
import warnings
import numpy as np
from scipy.sparse import coo_matrix

from ..utils import check_arrays, deprecated


###############################################################################
# General utilities
###############################################################################
def auc(x, y, reorder=False):
    """Compute Area Under the Curve (AUC) using the trapezoidal rule

    This is a general fuction, given points on a curve.  For computing the area
    under the ROC-curve, see :func:`auc_score`.

    Parameters
    ----------
    x : array, shape = [n]
        x coordinates.

    y : array, shape = [n]
        y coordinates.

    reorder : boolean, optional
        If True, assume that the curve is ascending in the case of ties, as for
        an ROC curve. If the curve is non-ascending, the result will be wrong.

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
    auc_score Computes the area under the ROC curve

    """
    # XXX: Consider using  ``scipy.integrate`` instead, or moving to
    # ``utils.extmath``
    x, y = check_arrays(x, y)
    if x.shape[0] < 2:
        raise ValueError('At least 2 points are needed to compute'
                         ' area under curve, but x.shape = %s' % x.shape)

    if reorder:
        # reorder the data points according to the x axis and using y to
        # break ties
        x, y = np.array(sorted(points for points in zip(x, y))).T
        h = np.diff(x)
    else:
        h = np.diff(x)
        if np.any(h < 0):
            h *= -1
            assert not np.any(h < 0), ("Reordering is not turned on, and "
                                       "The x array is not increasing: %s" % x)

    area = np.sum(h * (y[1:] + y[:-1])) / 2.0
    return area


def unique_labels(*lists_of_labels):
    """Extract an ordered array of unique labels"""
    labels = set().union(*(l.ravel() if hasattr(l, "ravel") else l
                           for l in lists_of_labels))
    return np.asarray(sorted(labels))


###############################################################################
# Binary classification loss
###############################################################################
def hinge_loss(y_true, pred_decision, pos_label=1, neg_label=-1):
    """Average hinge loss (non-regularized)

    Assuming labels in y_true are encoded with +1 and -1, when a prediction
    mistake is made, ``margin = y_true * pred_decision`` is always negative
    (since the signs disagree), implying ``1 - margin`` is always greater than
    1.  The cumulated hinge loss is therefore an upper bound of the number of
    mistakes made by the classifier.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True target (integers).

    pred_decision : array, shape = [n_samples] or [n_samples, n_classes]
        Predicted decisions, as output by decision_function (floats).

    Returns
    -------
    loss : float

    References
    ----------
    http://en.wikipedia.org/wiki/Hinge_loss

    Examples
    --------
    >>> from sklearn import svm
    >>> from sklearn.metrics import hinge_loss
    >>> X = [[0], [1]]
    >>> y = [-1, 1]
    >>> est = svm.LinearSVC(random_state=0)
    >>> est.fit(X, y)
    LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
         intercept_scaling=1, loss='l2', multi_class='ovr', penalty='l2',
         random_state=0, tol=0.0001, verbose=0)
    >>> pred_decision = est.decision_function([[-2], [3], [0.5]])
    >>> pred_decision  # doctest: +ELLIPSIS
    array([-2.18...,  2.36...,  0.09...])
    >>> hinge_loss([-1, 1, 1], pred_decision)  # doctest: +ELLIPSIS
    0.30...

    """
    # TODO: multi-class hinge-loss

    if pos_label != 1 or neg_label != -1:
        # the rest of the code assumes that positive and negative labels
        # are encoded as +1 and -1 respectively
        y_true = y_true.copy()
        y_true[y_true == pos_label] = 1
        y_true[y_true == neg_label] = -1

    margin = y_true * pred_decision
    losses = 1 - margin
    # The hinge doesn't penalize good enough predictions.
    losses[losses <= 0] = 0
    return np.mean(losses)


###############################################################################
# Binary classification scores
###############################################################################
def average_precision_score(y_true, y_score):
    """Compute average precision (AP) from prediction scores

    This score corresponds to the area under the precision-recall curve.

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True binary labels.

    y_score : array, shape = [n_samples]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or binary decisions.

    Returns
    -------
    average_precision : float

    References
    ----------
    http://en.wikipedia.org/wiki/Information_retrieval#Average_precision

    See also
    --------
    auc_score: Area under the ROC curve

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import average_precision_score
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> average_precision_score(y_true, y_scores)  # doctest: +ELLIPSIS
    0.79...

    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    return auc(recall, precision)


def auc_score(y_true, y_score):
    """Compute Area Under the Curve (AUC) from prediction scores

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        True binary labels.

    y_score : array, shape = [n_samples]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or binary decisions.

    Returns
    -------
    auc : float

    References
    ----------
    http://en.wikipedia.org/wiki/Receiver_operating_characteristic

    See also
    --------
    average_precision_score: Area under the precision-recall curve

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import auc_score
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> auc_score(y_true, y_scores)
    0.75

    """

    fpr, tpr, tresholds = roc_curve(y_true, y_score)
    return auc(fpr, tpr, reorder=True)


def matthews_corrcoef(y_true, y_pred):
    """Compute the Matthews correlation coefficient (MCC) for binary classes

    The Matthews correlation coefficient is used in machine learning as a
    measure of the quality of binary (two-class) classifications. It takes into
    account true and false positives and negatives and is generally regarded as
    a balanced measure which can be used even if the classes are of very
    different sizes. The MCC is in essence a correlation coefficient value
    between -1 and +1. A coefficient of +1 represents a perfect prediction, 0
    an average random prediction and -1 an inverse prediction.  The statistic
    is also known as the phi coefficient. [source: Wikipedia]

    Only in the binary case does this relate to information about true and
    false positives and negatives. See references below.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    Returns
    -------
    mcc : float
        The Matthews correlation coefficient (+1 represents a perfect
        prediction, 0 an average random prediction and -1 and inverse
        prediction).

    References
    ----------

    .. [1] `Baldi, Brunak, Chauvin, Andersen and Nielsen, (2000). Assessing the
       accuracy of prediction algorithms for classification: an overview
       <http://dx.doi.org/10.1093/bioinformatics/16.5.412>`_

    .. [2] `Wikipedia entry for the Matthews Correlation Coefficient
       <http://en.wikipedia.org/wiki/Matthews_correlation_coefficient>`_

    Examples
    --------
    >>> from sklearn.metrics import matthews_corrcoef
    >>> y_true = [+1, +1, +1, -1]
    >>> y_pred = [+1, -1, +1, +1]
    >>> matthews_corrcoef(y_true, y_pred)  # doctest: +ELLIPSIS
    -0.33...

    """
    mcc = np.corrcoef(y_true, y_pred)[0, 1]
    if np.isnan(mcc):
        return 0.
    else:
        return mcc


def precision_recall_curve(y_true, probas_pred):
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
    x axis.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets of binary classification in range {-1, 1} or {0, 1}.

    probas_pred : array, shape = [n_samples]
        Estimated probabilities or decision function.

    Returns
    -------
    precision : array, shape = [n + 1]
        Precision values.

    recall : array, shape = [n + 1]
        Recall values.

    thresholds : array, shape = [n]
        Thresholds on y_score used to compute precision and recall.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import precision_recall_curve
    >>> y_true = np.array([0, 0, 1, 1])
    >>> y_scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> precision, recall, threshold = precision_recall_curve(y_true, y_scores)
    >>> precision  # doctest: +ELLIPSIS
    array([ 0.66...,  0.5       ,  1.        ,  1.        ])
    >>> recall
    array([ 1. ,  0.5,  0.5,  0. ])
    >>> threshold
    array([ 0.35,  0.4 ,  0.8 ])

    """
    y_true = np.ravel(y_true)
    probas_pred = np.ravel(probas_pred)

    # Make sure input is boolean
    labels = np.unique(y_true)
    if np.all(labels == np.array([-1, 1])):
        # convert {-1, 1} to boolean {0, 1} repr
        y_true = y_true.copy()
        y_true[y_true == -1] = 0
    elif not np.all(labels == np.array([0, 1])):
        raise ValueError("y_true contains non binary labels: %r" % labels)

    # Sort pred_probas (and corresponding true labels) by pred_proba value
    decreasing_probas_indices = np.argsort(probas_pred, kind="mergesort")[::-1]
    probas_pred = probas_pred[decreasing_probas_indices]
    y_true = y_true[decreasing_probas_indices]

    # probas_pred typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate values for the beginning and end of the curve.
    distinct_value_indices = np.where(np.diff(probas_pred))[0] + 1
    threshold_idxs = np.hstack([0,
                                distinct_value_indices,
                                len(probas_pred)])

    # Initialize true and false positive counts, precision and recall
    total_positive = float(y_true.sum())
    tp_count, fp_count = 0., 0.  # Must remain floats to prevent int division
    precision = [1.]
    recall = [0.]
    thresholds = []

    # Iterate over indices which indicate distinct values (thresholds) of
    # probas_pred. Each of these threshold values will be represented in the
    # curve with a coordinate in precision-recall space. To calculate the
    # precision and recall associated with each point, we use these indices to
    # select all labels associated with the predictions. By incrementally
    # keeping track of the number of positive and negative labels seen so far,
    # we can calculate precision and recall.
    for l_idx, r_idx in izip(threshold_idxs[:-1], threshold_idxs[1:]):
        threshold_labels = y_true[l_idx:r_idx]
        n_at_threshold = r_idx - l_idx
        n_pos_at_threshold = threshold_labels.sum()
        n_neg_at_threshold = n_at_threshold - n_pos_at_threshold
        tp_count += n_pos_at_threshold
        fp_count += n_neg_at_threshold
        fn_count = total_positive - tp_count
        precision.append(tp_count / (tp_count + fp_count))
        recall.append(tp_count / (tp_count + fn_count))
        thresholds.append(probas_pred[l_idx])
        if tp_count == total_positive:
            break

    # sklearn expects these in reverse order
    thresholds = np.array(thresholds)[::-1]
    precision = np.array(precision)[::-1]
    recall = np.array(recall)[::-1]
    return precision, recall, thresholds


def roc_curve(y_true, y_score, pos_label=None):
    """Compute Receiver operating characteristic (ROC)

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        True binary labels in range {0, 1} or {-1, 1}.  If labels are not
        binary, pos_label should be explictly given.

    y_score : array, shape = [n_samples]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or binary decisions.

    pos_label : int
        Label considered as positive and others are considered negative.

    Returns
    -------
    fpr : array, shape = [>2]
        False Positive Rates.

    tpr : array, shape = [>2]
        True Positive Rates.

    thresholds : array, shape = [>2]
        Thresholds on ``y_score`` used to compute ``fpr`` and ``fpr``.

    Notes
    -----
    Since the thresholds are sorted from low to high values, they
    are reversed upon returning them to ensure they correspond to both ``fpr``
    and ``tpr``, which are sorted in reversed order during their calculation.

    References
    ----------
    http://en.wikipedia.org/wiki/Receiver_operating_characteristic

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([1, 1, 2, 2])
    >>> scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=2)
    >>> fpr
    array([ 0. ,  0.5,  0.5,  1. ])

    """
    y_true = np.ravel(y_true)
    y_score = np.ravel(y_score)
    classes = np.unique(y_true)

    # ROC only for binary classification if pos_label not given
    if (pos_label is None and
        not (np.all(classes == [0, 1]) or
             np.all(classes == [-1, 1]) or
             np.all(classes == [0]) or
             np.all(classes == [-1]) or
             np.all(classes == [1]))):
        raise ValueError("ROC is defined for binary classification only or "
                         "pos_label should be explicitly given")
    elif pos_label is None:
        pos_label = 1.

    # y_true will be transformed into a boolean vector
    y_true = (y_true == pos_label)
    n_pos = float(y_true.sum())
    n_neg = y_true.shape[0] - n_pos

    if n_pos == 0:
        warnings.warn("No positive samples in y_true, "
                      "true positive value should be meaningless")
        n_pos = np.nan
    if n_neg == 0:
        warnings.warn("No negative samples in y_true, "
                      "false positive value should be meaningless")
        n_neg = np.nan

    thresholds = np.unique(y_score)
    neg_value, pos_value = False, True

    tpr = np.empty(thresholds.size, dtype=np.float)  # True positive rate
    fpr = np.empty(thresholds.size, dtype=np.float)  # False positive rate

    # Build tpr/fpr vector
    current_pos_count = current_neg_count = sum_pos = sum_neg = idx = 0

    signal = np.c_[y_score, y_true]
    sorted_signal = signal[signal[:, 0].argsort(), :][::-1]
    last_score = sorted_signal[0][0]
    for score, value in sorted_signal:
        if score == last_score:
            if value == pos_value:
                current_pos_count += 1
            else:
                current_neg_count += 1
        else:
            tpr[idx] = (sum_pos + current_pos_count) / n_pos
            fpr[idx] = (sum_neg + current_neg_count) / n_neg
            sum_pos += current_pos_count
            sum_neg += current_neg_count
            current_pos_count = 1 if value == pos_value else 0
            current_neg_count = 1 if value == neg_value else 0
            idx += 1
            last_score = score
    else:
        tpr[-1] = (sum_pos + current_pos_count) / n_pos
        fpr[-1] = (sum_neg + current_neg_count) / n_neg

    thresholds = thresholds[::-1]

    if not (n_pos is np.nan or n_neg is np.nan):
        # add (0,0) and (1, 1)
        if not (fpr[0] == 0 and fpr[-1] == 1):
            fpr = np.r_[0., fpr, 1.]
            tpr = np.r_[0., tpr, 1.]
            thresholds = np.r_[thresholds[0] + 1, thresholds,
                               thresholds[-1] - 1]
        elif not fpr[0] == 0:
            fpr = np.r_[0., fpr]
            tpr = np.r_[0., tpr]
            thresholds = np.r_[thresholds[0] + 1, thresholds]
        elif not fpr[-1] == 1:
            fpr = np.r_[fpr, 1.]
            tpr = np.r_[tpr, 1.]
            thresholds = np.r_[thresholds, thresholds[-1] - 1]
    elif fpr.shape[0] == 2:
        # trivial decisions, add (0,0)
        fpr = np.array([0.0, fpr[0], fpr[1]])
        tpr = np.array([0.0, tpr[0], tpr[1]])
        # trivial decisions, add (0,0) and (1,1)
    elif fpr.shape[0] == 1:
        fpr = np.array([0.0, fpr[0], 1.0])
        tpr = np.array([0.0, tpr[0], 1.0])

    if n_pos is np.nan:
        tpr[0] = np.nan

    if n_neg is np.nan:
        fpr[0] = np.nan

    return fpr, tpr, thresholds


###############################################################################
# Multiclass general function
###############################################################################
def confusion_matrix(y_true, y_pred, labels=None):
    """Compute confusion matrix to evaluate the accuracy of a classification

    By definition a confusion matrix :math:`C` is such that :math:`C_{i, j}`
    is equal to the number of observations known to be in group :math:`i` but
    predicted to be in group :math:`j`.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    labels : array, shape = [n_classes]
        List of all labels occuring in the dataset.
        If none is given, those that appear at least once
        in ``y_true`` or ``y_pred`` are used.

    Returns
    -------
    C : array, shape = [n_classes, n_classes]
        Confusion matrix

    References
    ----------
    http://en.wikipedia.org/wiki/Confusion_matrix

    Examples
    --------
    >>> from sklearn.metrics import confusion_matrix
    >>> y_true = [2, 0, 2, 2, 0, 1]
    >>> y_pred = [0, 0, 2, 2, 0, 2]
    >>> confusion_matrix(y_true, y_pred)
    array([[2, 0, 0],
           [0, 0, 1],
           [1, 0, 2]])

    """
    if labels is None:
        labels = unique_labels(y_true, y_pred)
    else:
        labels = np.asarray(labels, dtype=np.int)

    n_labels = labels.size
    label_to_ind = dict((y, x) for x, y in enumerate(labels))
    # convert yt, yp into index
    y_pred = np.array([label_to_ind.get(x, n_labels + 1) for x in y_pred])
    y_true = np.array([label_to_ind.get(x, n_labels + 1) for x in y_true])

    # intersect y_pred, y_true with labels, eliminate items not in labels
    ind = np.logical_and(y_pred < n_labels, y_true < n_labels)
    y_pred = y_pred[ind]
    y_true = y_true[ind]

    CM = np.asarray(coo_matrix((np.ones(y_true.shape[0]), (y_true, y_pred)),
                               shape=(n_labels, n_labels),
                               dtype=np.int).todense())
    return CM


###############################################################################
# Multiclass loss function
###############################################################################
def zero_one_loss(y_true, y_pred, normalize=True):
    """Zero-One classification loss

    If normalize is ``True``, return the fraction of misclassifications
    (float), else it returns the number of misclassifications (int). The best
    performance is 0.

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

    normalize : bool, optional
        If ``False`` (default), return the number of misclassifications.
        Otherwise, return the fraction of misclassifications.

    Returns
    -------
    loss : float or int,
        If ``normalize == True``, return the fraction of misclassifications
        (float), else it returns the number of misclassifications (int).

    Examples
    --------
    >>> from sklearn.metrics import zero_one_loss
    >>> y_pred = [1, 2, 3, 4]
    >>> y_true = [2, 2, 3, 4]
    >>> zero_one_loss(y_true, y_pred)
    0.25
    >>> zero_one_loss(y_true, y_pred, normalize=False)
    1

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    if not normalize:
        return np.sum(y_pred != y_true)
    else:
        return np.mean(y_pred != y_true)


@deprecated("Function 'zero_one' has been renamed to "
            "'zero_one_loss' and will be removed in release 0.15."
            "Default behavior is changed from 'normalize=False' to "
            "'normalize=True'")
def zero_one(y_true, y_pred, normalize=False):
    """Zero-One classification loss

    If normalize is ``True``, return the fraction of misclassifications
    (float), else it returns the number of misclassifications (int). The best
    performance is 0.

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

    normalize : bool, optional
        If ``False`` (default), return the number of misclassifications.
        Otherwise, return the fraction of misclassifications.

    Returns
    -------
    loss : float
        If normalize is True, return the fraction of misclassifications
        (float), else it returns the number of misclassifications (int).


    Examples
    --------
    >>> from sklearn.metrics import zero_one
    >>> y_pred = [1, 2, 3, 4]
    >>> y_true = [2, 2, 3, 4]
    >>> zero_one(y_true, y_pred)
    1
    >>> zero_one(y_true, y_pred, normalize=True)
    0.25

    """
    return zero_one_loss(y_true, y_pred, normalize)


###############################################################################
# Multiclass score functions
###############################################################################
def accuracy_score(y_true, y_pred):
    """Accuracy classification score

    Parameters
    ----------
    y_true : array-like, shape = n_samples
        Ground truth (correct) labels.

    y_pred : array-like, shape = n_samples
        Predicted labels, as returned by a classifier.

    Returns
    -------
    score : float
        The fraction of correct predictions in ``y_pred``.
        The best performance is 1.

    See also
    --------
    zero_one_loss Zero-One classification loss

    Examples
    --------
    >>> from sklearn.metrics import accuracy_score
    >>> y_pred = [0, 2, 1, 3]
    >>> y_true = [0, 1, 2, 3]
    >>> accuracy_score(y_true, y_pred)
    0.5

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.mean(y_pred == y_true)


def f1_score(y_true, y_pred, labels=None, pos_label=1, average='weighted'):
    """Compute the F1 score, also known as balanced F-score or F-measure

    The F1 score can be interpreted as a weighted average of the precision and
    recall, where an F1 score reaches its best value at 1 and worst score at 0.
    The relative contribution of precision and recall to the F1 score are
    equal. The formula for the F1 score is::

        F1 = 2 * (precision * recall) / (precision + recall)

    In the multi-class case, this is the weighted average of the F1 score of
    each class.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    labels : array
        Integer array of labels.

    pos_label : int
        In the binary classification case, give the label of the positive class
        (default is 1). Everything else but ``pos_label`` is considered to
        belong to the negative class.  Set to ``None`` in the case of
        multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted' (default)]
        In the multiclass classification case, this determines the type of
        averaging performed on the data.

        None:
            Do not perform any averaging, return the score for each class.
        'macro':
            Average over classes (does not take imbalance into account).
        'micro':
            Average over instances (takes imbalance into account).  This
            implies that ``precision == recall == F1``.
        'weighted':
            Average weighted by support (takes imbalance into account).  Can
            result in F-score that is not between precision and recall.

    Returns
    -------
    f1_score : float or array of float, shape = [n_unique_labels]
        F1 score of the positive class in binary classification or weighted
        average of the F1 scores of each class for the multiclass task.

    References
    ----------
    http://en.wikipedia.org/wiki/F1_score

    Examples
    --------
    In the binary case:

    >>> from sklearn.metrics import f1_score
    >>> y_pred = [0, 1, 0, 0]
    >>> y_true = [0, 1, 0, 1]
    >>> f1_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.666...

    In the multiclass case:

    >>> from sklearn.metrics import f1_score
    >>> y_true = [0, 1, 2, 0, 1, 2]
    >>> y_pred = [0, 2, 1, 0, 0, 1]
    >>> f1_score(y_true, y_pred, average='macro')  # doctest: +ELLIPSIS
    0.26...
    >>> f1_score(y_true, y_pred, average='micro')  # doctest: +ELLIPSIS
    0.33...
    >>> f1_score(y_true, y_pred, average='weighted')  # doctest: +ELLIPSIS
    0.26...
    >>> f1_score(y_true, y_pred, average=None)
    array([ 0.8,  0. ,  0. ])

    """
    return fbeta_score(y_true, y_pred, 1, labels=labels,
                       pos_label=pos_label, average=average)


def fbeta_score(y_true, y_pred, beta, labels=None, pos_label=1,
                average='weighted'):
    """Compute the F-beta score

    The F-beta score is the weighted harmonic mean of precision and recall,
    reaching its optimal value at 1 and its worst value at 0.

    The `beta` parameter determines the weight of precision in the combined
    score. ``beta < 1`` lends more weight to precision, while ``beta > 1``
    favors precision (``beta == 0`` considers only precision, ``beta == inf``
    only recall).

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    beta: float
        Weight of precision in harmonic mean.

    labels : array
        Integer array of labels.

    pos_label : int
        In the binary classification case, give the label of the positive class
        (default is 1). Everything else but ``pos_label`` is considered to
        belong to the negative class.  Set to ``None`` in the case of
        multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted' (default)]
        In the multiclass classification case, this determines the type of
        averaging performed on the data.

        ``None``:
            Do not perform any averaging, return the scores for each class.
        ``'macro'``:
            Average over classes (does not take imbalance into account).
        ``'micro'``:
            Average over instances (takes imbalance into account).  This
            implies that ``precision == recall == F1``.
        ``'weighted'``:
            Average weighted by support (takes imbalance into account).  Can
            result in F-score that is not between precision and recall.
            Do not perform any averaging, return the score for each class.

    Returns
    -------
    fbeta_score : float (if average is not None) or array of float, shape =\
        [n_unique_labels]
        F-beta score of the positive class in binary classification or weighted
        average of the F-beta score of each class for the multiclass task.

    References
    ----------
    R. Baeza-Yates and B. Ribeiro-Neto (2011). Modern Information Retrieval.
    Addison Wesley, pp. 327-328.

    http://en.wikipedia.org/wiki/F1_score

    Examples
    --------
    In the binary case:

    >>> from sklearn.metrics import fbeta_score
    >>> y_pred = [0, 1, 0, 0]
    >>> y_true = [0, 1, 0, 1]
    >>> fbeta_score(y_true, y_pred, beta=0.5)  # doctest: +ELLIPSIS
    0.83...
    >>> fbeta_score(y_true, y_pred, beta=1)  # doctest: +ELLIPSIS
    0.66...
    >>> fbeta_score(y_true, y_pred, beta=2)  # doctest: +ELLIPSIS
    0.55...

    In the multiclass case:

    >>> from sklearn.metrics import fbeta_score
    >>> y_true = [0, 1, 2, 0, 1, 2]
    >>> y_pred = [0, 2, 1, 0, 0, 1]
    >>> fbeta_score(y_true, y_pred, average='macro', beta=0.5)\
        # doctest: +ELLIPSIS
    0.23...
    >>> fbeta_score(y_true, y_pred, average='micro', beta=0.5)\
        # doctest: +ELLIPSIS
    0.33...
    >>> fbeta_score(y_true, y_pred, average='weighted', beta=0.5)\
        # doctest: +ELLIPSIS
    0.23...
    >>> fbeta_score(y_true, y_pred, average=None, beta=0.5)\
        # doctest: +ELLIPSIS
    array([ 0.71...,  0.        ,  0.        ])

    """
    _, _, f, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 beta=beta,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return f


def precision_recall_fscore_support(y_true, y_pred, beta=1.0, labels=None,
                                    pos_label=1, average=None):
    """Compute precision, recall, F-measure and support for each class

    The precision is the ratio ``tp / (tp + fp)`` where ``tp`` is the number of
    true positives and ``fp`` the number of false positives. The precision is
    intuitively the ability of the classifier not to label as positive a sample
    that is negative.

    The recall is the ratio ``tp / (tp + fn)`` where ``tp`` is the number of
    true positives and ``fn`` the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The F-beta score can be interpreted as a weighted harmonic mean of
    the precision and recall, where an F-beta score reaches its best
    value at 1 and worst score at 0.

    The F-beta score weights recall more than precision by a factor of
    ``beta``. ``beta == 1.0`` means recall and precsion are equally important.

    The support is the number of occurrences of each class in ``y_true``.

    If ``pos_label is None``, this function returns the average precision,
    recall and F-measure if ``average`` is one of ``'micro'``, ``'macro'``,
    ``'weighted'``.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    beta : float, 1.0 by default
        The strength of recall versus precision in the F-score.

    labels : array
        Integer array of labels.

    pos_label : int
        In the binary classification case, give the label of the positive class
        (default is 1). Everything else but ``pos_label`` is considered to
        belong to the negative class.  Set to ``None`` in the case of
        multiclass classification.

    average : string, [None (default), 'micro', 'macro', 'weighted']
        In the multiclass classification case, this determines the type of
        averaging performed on the data.

        ``None``:
            Do not perform any averaging, return the scores for each class.
        ``'macro'``:
            Average over classes (does not take imbalance into account).
        ``'micro'``:
            Average over instances (takes imbalance into account).  This
            implies that ``precision == recall == F1``.
        ``'weighted'``:
            Average weighted by support (takes imbalance into account).  Can
            result in F-score that is not between precision and recall.

    Returns
    -------
    precision: float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    recall: float (if average is not None) or array of float, , shape =\
        [n_unique_labels]

    f1_score: float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    support: int (if average is not None) or array of int, shape =\
        [n_unique_labels]

    References
    ----------
    http://en.wikipedia.org/wiki/Precision_and_recall

    Examples
    --------
    In the binary case:

    >>> from sklearn.metrics import precision_recall_fscore_support
    >>> y_pred = [0, 1, 0, 0]
    >>> y_true = [0, 1, 0, 1]
    >>> p, r, f, s = precision_recall_fscore_support(y_true, y_pred, beta=0.5)
    >>> p  # doctest: +ELLIPSIS
    array([ 0.66...,  1.        ])
    >>> r
    array([ 1. ,  0.5])
    >>> f  # doctest: +ELLIPSIS
    array([ 0.71...,  0.83...])
    >>> s  # doctest: +ELLIPSIS
    array([2, 2]...)

    In the multiclass case:

    >>> from sklearn.metrics import precision_recall_fscore_support
    >>> y_true = np.array([0, 1, 2, 0, 1, 2])
    >>> y_pred = np.array([0, 2, 1, 0, 0, 1])
    >>> precision_recall_fscore_support(y_true, y_pred, average='macro')\
        # doctest: +ELLIPSIS
    (0.22..., 0.33..., 0.26..., None)
    >>> precision_recall_fscore_support(y_true, y_pred, average='micro')\
        # doctest: +ELLIPSIS
    (0.33..., 0.33..., 0.33..., None)
    >>> precision_recall_fscore_support(y_true, y_pred, average='weighted')\
        # doctest: +ELLIPSIS
    (0.22..., 0.33..., 0.26..., None)

    """
    if beta <= 0:
        raise ValueError("beta should be >0 in the F-beta score")

    y_true, y_pred = check_arrays(y_true, y_pred)
    if labels is None:
        labels = unique_labels(y_true, y_pred)
    else:
        labels = np.asarray(labels, dtype=np.int)

    n_labels = labels.size
    true_pos = np.zeros(n_labels, dtype=np.double)
    false_pos = np.zeros(n_labels, dtype=np.double)
    false_neg = np.zeros(n_labels, dtype=np.double)
    support = np.zeros(n_labels, dtype=np.long)

    for i, label_i in enumerate(labels):
        true_pos[i] = np.sum(y_pred[y_true == label_i] == label_i)
        false_pos[i] = np.sum(y_pred[y_true != label_i] == label_i)
        false_neg[i] = np.sum(y_pred[y_true == label_i] != label_i)
        support[i] = np.sum(y_true == label_i)

    try:
        # oddly, we may get an "invalid" rather than a "divide" error here
        old_err_settings = np.seterr(divide='ignore', invalid='ignore')

        # precision and recall
        precision = true_pos / (true_pos + false_pos)
        recall = true_pos / (true_pos + false_neg)

        # handle division by 0.0 in precision and recall
        precision[(true_pos + false_pos) == 0.0] = 0.0
        recall[(true_pos + false_neg) == 0.0] = 0.0

        # fbeta score
        beta2 = beta ** 2
        fscore = (1 + beta2) * (precision * recall) / (
            beta2 * precision + recall)

        # handle division by 0.0 in fscore
        fscore[(precision + recall) == 0.0] = 0.0
    finally:
        np.seterr(**old_err_settings)

    if not average:
        return precision, recall, fscore, support

    elif n_labels == 2 and pos_label is not None:
        if pos_label not in labels:
            raise ValueError("pos_label=%d is not a valid label: %r" %
                             (pos_label, labels))
        pos_label_idx = list(labels).index(pos_label)
        return (precision[pos_label_idx], recall[pos_label_idx],
                fscore[pos_label_idx], support[pos_label_idx])
    else:
        average_options = (None, 'micro', 'macro', 'weighted')
        if average == 'micro':
            avg_precision = true_pos.sum() / (true_pos.sum() +
                                              false_pos.sum())
            avg_recall = true_pos.sum() / (true_pos.sum() + false_neg.sum())
            avg_fscore = (1 + beta2) * (avg_precision * avg_recall) / \
                         (beta2 * avg_precision + avg_recall)
        elif average == 'macro':
            avg_precision = np.mean(precision)
            avg_recall = np.mean(recall)
            avg_fscore = np.mean(fscore)
        elif average == 'weighted':
            avg_precision = np.average(precision, weights=support)
            avg_recall = np.average(recall, weights=support)
            avg_fscore = np.average(fscore, weights=support)
        else:
            raise ValueError('average has to be one of ' +
                             str(average_options))

        return avg_precision, avg_recall, avg_fscore, None


def precision_score(y_true, y_pred, labels=None, pos_label=1,
                    average='weighted'):
    """Compute the precision

    The precision is the ratio ``tp / (tp + fp)`` where ``tp`` is the number of
    true positives and ``fp`` the number of false positives. The precision is
    intuitively the ability of the classifier not to label as positive a sample
    that is negative.

    The best value is 1 and the worst value is 0.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    labels : array
        Integer array of labels.

    pos_label : int
        In the binary classification case, give the label of the positive class
        (default is 1). Everything else but ``pos_label`` is considered to
        belong to the negative class.  Set to ``None`` in the case of
        multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted' (default)]
        In the multiclass classification case, this determines the type of
        averaging performed on the data.

        ``None``:
            Do not perform any averaging, return the scores for each class.
        ``'macro'``:
            Average over classes (does not take imbalance into account).
        ``'micro'``:
            Average over instances (takes imbalance into account).  This
            implies that ``precision == recall == F1``.
        ``'weighted'``:
            Average weighted by support (takes imbalance into account).  Can
            result in F-score that is not between precision and recall.

    Returns
    -------
    precision : float (if average is not None) or array of float, shape =\
        [n_unique_labels]
        Precision of the positive class in binary classification or weighted
        average of the precision of each class for the multiclass task.

    Examples
    --------
    In the binary case:

    >>> from sklearn.metrics import precision_score
    >>> y_pred = [0, 1, 0, 0]
    >>> y_true = [0, 1, 0, 1]
    >>> precision_score(y_true, y_pred)
    1.0

    In the multiclass case:

    >>> from sklearn.metrics import precision_score
    >>> y_true = [0, 1, 2, 0, 1, 2]
    >>> y_pred = [0, 2, 1, 0, 0, 1]
    >>> precision_score(y_true, y_pred, average='macro')  # doctest: +ELLIPSIS
    0.22...
    >>> precision_score(y_true, y_pred, average='micro')  # doctest: +ELLIPSIS
    0.33...
    >>> precision_score(y_true, y_pred, average='weighted')\
        # doctest: +ELLIPSIS
    0.22...
    >>> precision_score(y_true, y_pred, average=None)  # doctest: +ELLIPSIS
    array([ 0.66...,  0.        ,  0.        ])

    """
    p, _, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return p


def recall_score(y_true, y_pred, labels=None, pos_label=1, average='weighted'):
    """Compute the recall

    The recall is the ratio ``tp / (tp + fn)`` where ``tp`` is the number of
    true positives and ``fn`` the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The best value is 1 and the worst value is 0.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    labels : array
        Integer array of labels.

    pos_label : int
        In the binary classification case, give the label of the positive class
        (default is 1). Everything else but ``pos_label`` is considered to
        belong to the negative class.  Set to ``None`` in the case of
        multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted' (default)]
        In the multiclass classification case, this determines the type of
        averaging performed on the data.

        ``None``:
            Do not perform any averaging, return the scores for each class.
        ``'macro'``:
            Average over classes (does not take imbalance into account).
        ``'micro'``:
            Average over instances (takes imbalance into account).  This
            implies that ``precision == recall == F1``.
        ``'weighted'``:
            Average weighted by support (takes imbalance into account).  Can
            result in F-score that is not between precision and recall.

    Returns
    -------
    recall : float (if average is not None) or array of float, shape =\
        [n_unique_labels]
        Recall of the positive class in binary classification or weighted
        average of the recall of each class for the multiclass task.

    Examples
    --------
    In the binary case:

    >>> from sklearn.metrics import recall_score
    >>> y_pred = [0, 1, 0, 0]
    >>> y_true = [0, 1, 0, 1]
    >>> recall_score(y_true, y_pred)
    0.5

    In the multiclass case:

    >>> from sklearn.metrics import recall_score
    >>> y_true = [0, 1, 2, 0, 1, 2]
    >>> y_pred = [0, 2, 1, 0, 0, 1]
    >>> recall_score(y_true, y_pred, average='macro')  # doctest: +ELLIPSIS
    0.33...
    >>> recall_score(y_true, y_pred, average='micro')  # doctest: +ELLIPSIS
    0.33...
    >>> recall_score(y_true, y_pred, average='weighted')  # doctest: +ELLIPSIS
    0.33...
    >>> recall_score(y_true, y_pred, average=None)
    array([ 1.,  0.,  0.])

    """
    _, r, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return r


@deprecated("Function zero_one_score has been renamed to "
            'accuracy_score'" and will be removed in release 0.15.")
def zero_one_score(y_true, y_pred):
    """Zero-one classification score (accuracy)

    Parameters
    ----------
    y_true : array-like, shape = n_samples
        Ground truth (correct) labels.

    y_pred : array-like, shape = n_samples
        Predicted labels, as returned by a classifier.

    Returns
    -------
    score : float
        Fraction of correct predictions in ``y_pred``. The best performance is
        1.

    """
    return accuracy_score(y_true, y_pred)


###############################################################################
# Multiclass utility function
###############################################################################
def classification_report(y_true, y_pred, labels=None, target_names=None):
    """Build a text report showing the main classification metrics

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        Ground truth (correct) target values.

    y_pred : array, shape = [n_samples]
        Estimated targets as returned by a classifier.

    labels : array, shape = [n_labels]
        Optional list of label indices to include in the report.

    target_names : list of strings
        Optional display names matching the labels (same order).

    Returns
    -------
    report : string
        Text summary of the precision, recall, F1 score for each class.

    Examples
    --------
    >>> from sklearn.metrics import classification_report
    >>> y_true = [0, 1, 2, 2, 0]
    >>> y_pred = [0, 0, 2, 2, 0]
    >>> target_names = ['class 0', 'class 1', 'class 2']
    >>> print(classification_report(y_true, y_pred, target_names=target_names))
                 precision    recall  f1-score   support
    <BLANKLINE>
        class 0       0.67      1.00      0.80         2
        class 1       0.00      0.00      0.00         1
        class 2       1.00      1.00      1.00         2
    <BLANKLINE>
    avg / total       0.67      0.80      0.72         5
    <BLANKLINE>

    """

    if labels is None:
        labels = unique_labels(y_true, y_pred)
    else:
        labels = np.asarray(labels, dtype=np.int)

    last_line_heading = 'avg / total'

    if target_names is None:
        width = len(last_line_heading)
        target_names = ['%d' % l for l in labels]
    else:
        width = max(len(cn) for cn in target_names)
        width = max(width, len(last_line_heading))

    headers = ["precision", "recall", "f1-score", "support"]
    fmt = '%% %ds' % width  # first column: class name
    fmt += '  '
    fmt += ' '.join(['% 9s' for _ in headers])
    fmt += '\n'

    headers = [""] + headers
    report = fmt % tuple(headers)
    report += '\n'

    p, r, f1, s = precision_recall_fscore_support(y_true, y_pred,
                                                  labels=labels,
                                                  average=None)

    for i, label in enumerate(labels):
        values = [target_names[i]]
        for v in (p[i], r[i], f1[i]):
            values += ["%0.2f" % float(v)]
        values += ["%d" % int(s[i])]
        report += fmt % tuple(values)

    report += '\n'

    # compute averages
    values = [last_line_heading]
    for v in (np.average(p, weights=s),
              np.average(r, weights=s),
              np.average(f1, weights=s)):
        values += ["%0.2f" % float(v)]
    values += ['%d' % np.sum(s)]
    report += fmt % tuple(values)
    return report


###############################################################################
# Regression loss functions
###############################################################################
def mean_absolute_error(y_true, y_pred):
    """Mean absolute error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    Returns
    -------
    loss : float
        A positive floating point value (the best value is 0.0).

    Examples
    --------
    >>> from sklearn.metrics import mean_absolute_error
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> mean_absolute_error(y_true, y_pred)
    0.5
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> mean_absolute_error(y_true, y_pred)
    0.75

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.mean(np.abs(y_pred - y_true))


def mean_squared_error(y_true, y_pred):
    """Mean squared error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    Returns
    -------
    loss : float
        A positive floating point value (the best value is 0.0).

    Examples
    --------
    >>> from sklearn.metrics import mean_squared_error
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> mean_squared_error(y_true, y_pred)
    0.375
    >>> y_true = [[0.5, 1],[-1, 1],[7, -6]]
    >>> y_pred = [[0, 2],[-1, 2],[8, -5]]
    >>> mean_squared_error(y_true, y_pred)  # doctest: +ELLIPSIS
    0.708...

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.mean((y_pred - y_true) ** 2)


###############################################################################
# Regression score functions
###############################################################################
def explained_variance_score(y_true, y_pred):
    """Explained variance regression score function

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like
        Ground truth (correct) target values.

    y_pred : array-like
        Estimated target values.

    Returns
    -------
    score : float
        The explained variance.

    Notes
    -----
    This is not a symmetric function.

    Examples
    --------
    >>> from sklearn.metrics import explained_variance_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> explained_variance_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.957...

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    numerator = np.var(y_true - y_pred)
    denominator = np.var(y_true)
    if denominator == 0.0:
        if numerator == 0.0:
            return 1.0
        else:
            # arbitary set to zero to avoid -inf scores, having a constant
            # y_true is not interesting for scoring a regression anyway
            return 0.0
    return 1 - numerator / denominator


def r2_score(y_true, y_pred):
    """R^2 (coefficient of determination) regression score function

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    Returns
    -------
    z : float
        The R^2 score

    Notes
    -----
    This is not a symmetric function.

    References
    ----------
    http://en.wikipedia.org/wiki/Coefficient_of_determination

    Examples
    --------
    >>> from sklearn.metrics import r2_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.948...
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.938...

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    if len(y_true) == 1:
        raise ValueError("r2_score can only be computed given more than one"
                         " sample.")
    numerator = ((y_true - y_pred) ** 2).sum()
    denominator = ((y_true - y_true.mean(axis=0)) ** 2).sum()

    if denominator == 0.0:
        if numerator == 0.0:
            return 1.0
        else:
            # arbitary set to zero to avoid -inf scores, having a constant
            # y_true is not interesting for scoring a regression anyway
            return 0.0

    return 1 - numerator / denominator
