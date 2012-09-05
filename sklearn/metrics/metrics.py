"""Utilities to evaluate the predictive performance of models

Functions named as *_score return a scalar value to maximize: the higher the
better

Function named as *_loss return a scalar value to minimize: the lower the
better
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

import numpy as np

from ..utils import check_arrays
from ..utils import deprecated


def unique_labels(*lists_of_labels):
    """Extract an ordered array of unique labels"""
    labels = set()
    for l in lists_of_labels:
        if hasattr(l, 'ravel'):
            l = l.ravel()
        labels |= set(l)
    return np.unique(sorted(labels))


def confusion_matrix(y_true, y_pred, labels=None):
    """Compute confusion matrix to evaluate the accuracy of a classification

    By definition a confusion matrix cm is such that cm[i, j] is equal
    to the number of observations known to be in group i but predicted
    to be in group j.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        estimated targets

    labels : array, shape = [n_classes]
        lists all labels occuring in the dataset.
        If none is given, those that appear at least once
        in y_true or y_pred are used.

    Returns
    -------
    CM : array, shape = [n_classes, n_classes]
        confusion matrix

    References
    ----------
    http://en.wikipedia.org/wiki/Confusion_matrix
    """
    if labels is None:
        labels = unique_labels(y_true, y_pred)
    else:
        labels = np.asarray(labels, dtype=np.int)

    n_labels = labels.size
    label_to_ind = dict((y, x) for x, y in enumerate(labels))

    if n_labels >= 15:
        CM = np.zeros((n_labels, n_labels), dtype=np.long)
        for yt, yp in zip(y_true, y_pred):
            CM[label_to_ind[yt], label_to_ind[yp]] += 1
    else:
        CM = np.empty((n_labels, n_labels), dtype=np.long)
        for i, label_i in enumerate(labels):
            for j, label_j in enumerate(labels):
                CM[i, j] = np.sum(
                    np.logical_and(y_true == label_i, y_pred == label_j))

    return CM


def roc_curve(y_true, y_score):
    """compute Receiver operating characteristic (ROC)

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        true binary labels

    y_score : array, shape = [n_samples]
        target scores, can either be probability estimates of
        the positive class, confidence values, or binary decisions.

    Returns
    -------
    fpr : array, shape = [>2]
        False Positive Rates

    tpr : array, shape = [>2]
        True Positive Rates

    thresholds : array, shape = [>2]
        Thresholds on y_score used to compute fpr and tpr.

        *Note*: Since the thresholds are sorted from low to high values,
        they are reversed upon returning them to ensure they
        correspond to both fpr and tpr, which are sorted in reversed order
        during their calculation.


    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([1, 1, 2, 2])
    >>> scores = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, scores)
    >>> fpr
    array([ 0. ,  0.5,  0.5,  1. ])

    References
    ----------
    http://en.wikipedia.org/wiki/Receiver_operating_characteristic

    """
    y_true = np.ravel(y_true)
    classes = np.unique(y_true)

    # ROC only for binary classification
    if classes.shape[0] != 2:
        raise ValueError("ROC is defined for binary classification only")

    y_score = np.ravel(y_score)

    n_pos = float(np.sum(y_true == classes[1]))  # nb of true positive
    n_neg = float(np.sum(y_true == classes[0]))  # nb of true negative

    thresholds = np.unique(y_score)
    neg_value, pos_value = classes[0], classes[1]

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

    # hard decisions, add (0,0)
    if fpr.shape[0] == 2:
        fpr = np.array([0.0, fpr[0], fpr[1]])
        tpr = np.array([0.0, tpr[0], tpr[1]])
    # trivial decisions, add (0,0) and (1,1)
    elif fpr.shape[0] == 1:
        fpr = np.array([0.0, fpr[0], 1.0])
        tpr = np.array([0.0, tpr[0], 1.0])

    return fpr, tpr, thresholds[::-1]


def average_precision_score(y_true, y_score):
    """Compute average precision (AP) from prediction scores.

    This score corresponds to the area under the precision-recall curve.

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        true binary labels

    y_score : array, shape = [n_samples]
        target scores, can either be probability estimates of
        the positive class, confidence values, or binary decisions.

    Returns
    -------
    average_precision : float

    References
    ----------
    http://en.wikipedia.org/wiki/Information_retrieval#Average_precision


    See also
    --------
    auc_score: Area under the ROC curve
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    return auc(recall, precision)


def auc_score(y_true, y_score):
    """Compute Area Under the Curve (AUC) from prediction scores.

    Note: this implementation is restricted to the binary classification task.

    Parameters
    ----------

    y_true : array, shape = [n_samples]
        true binary labels

    y_score : array, shape = [n_samples]
        target scores, can either be probability estimates of
        the positive class, confidence values, or binary decisions.

    Returns
    -------
    auc : float

    References
    ----------
    http://en.wikipedia.org/wiki/Receiver_operating_characteristic

    See also
    --------
    average_precision_score: Area under the precision-recall curve
    """

    fpr, tpr, tresholds = roc_curve(y_true, y_score)
    return auc(fpr, tpr)


def auc(x, y):
    """Compute Area Under the Curve (AUC) using the trapezoidal rule

    This is a general fuction, given points on a curve.
    For computing the area under the ROC-curve, see auc_score.

    Parameters
    ----------
    x : array, shape = [n]
        x coordinates

    y : array, shape = [n]
        y coordinates

    Returns
    -------
    auc : float

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import metrics
    >>> y = np.array([1, 1, 2, 2])
    >>> pred = np.array([0.1, 0.4, 0.35, 0.8])
    >>> fpr, tpr, thresholds = metrics.roc_curve(y, pred)
    >>> metrics.auc(fpr, tpr)
    0.75

    See also
    --------
    auc_score Computes the area under the ROC curve

    """
    x, y = check_arrays(x, y)
    if x.shape[0] != y.shape[0]:
        raise ValueError('x and y should have the same shape'
                         ' to compute area under curve,'
                         ' but x.shape = %s and y.shape = %s.'
                         % (x.shape, y.shape))
    if x.shape[0] < 2:
        raise ValueError('At least 2 points are needed to compute'
                         ' area under curve, but x.shape = %s' % x.shape)

    # reorder the data points according to the x axis and using y to break ties
    x, y = np.array(sorted(points for points in zip(x, y))).T

    h = np.diff(x)
    area = np.sum(h * (y[1:] + y[:-1])) / 2.0
    return area


def precision_score(y_true, y_pred, labels=None, pos_label=1,
                    average='weighted'):
    """Compute the precision

    The precision is the ratio :math:`tp / (tp + fp)` where tp is the
    number of true positives and fp the number of false positives. The
    precision is intuitively the ability of the classifier not to
    label as positive a sample that is negative.

    The best value is 1 and the worst value is 0.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Predicted targets

    labels : array
        Integer array of labels

    pos_label : int
        In the binary classification case, give the label of the positive
        class (default is 1). Everything else but 'pos_label'
        is considered to belong to the negative class.
        Set to None in the case of multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted'(default)]
        In the multiclass classification case, this determines the
        type of averaging performed on the data.

        macro:
            Average over classes (does not take imbalance into account).
        micro:
            Average over instances (takes imbalance into account).
            This implies that ``precision == recall == f1``
        weighted:
            Average weighted by support (takes imbalance into account).
            Can result in f1 score that is not between precision and recall.

    Returns
    -------
    precision : float
        Precision of the positive class in binary classification or
        weighted average of the precision of each class for the
        multiclass task

    """
    p, _, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return p


def recall_score(y_true, y_pred, labels=None, pos_label=1, average='weighted'):
    """Compute the recall

    The recall is the ratio :math:`tp / (tp + fn)` where tp is the number of
    true positives and fn the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The best value is 1 and the worst value is 0.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Predicted targets

    labels : array
        Integer array of labels

    pos_label : int
        In the binary classification case, give the label of the positive
        class (default is 1). Everything else but 'pos_label'
        is considered to belong to the negative class.
        Set to None in the case of multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted'(default)]
        In the multiclass classification case, this determines the
        type of averaging performed on the data.

        macro:
            Average over classes (does not take imbalance into account).
        micro:
            Average over instances (takes imbalance into account).
            This implies that ``precision == recall == f1``
        weighted:
            Average weighted by support (takes imbalance into account).
            Can result in f1 score that is not between precision and recall.

    Returns
    -------
    recall : float
        Recall of the positive class in binary classification or weighted
        average of the recall of each class for the multiclass task.

    """
    _, r, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return r


def fbeta_score(y_true, y_pred, beta, labels=None, pos_label=1,
                average='weighted'):
    """Compute fbeta score

    The F_beta score is the weighted harmonic mean of precision and recall,
    reaching its optimal value at 1 and its worst value at 0.

    The beta parameter determines the weight of precision in the combined
    score. ``beta < 1`` lends more weight to precision, while ``beta > 1``
    favors precision (``beta == 0`` considers only precision, ``beta == inf``
    only recall).

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Predicted targets

    beta: float
        Weight of precision in harmonic mean.

    labels : array
        Integer array of labels

    pos_label : int
        In the binary classification case, give the label of the positive
        class (default is 1). Everything else but 'pos_label'
        is considered to belong to the negative class.
        Set to None in the case of multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted'(default)]
        In the multiclass classification case, this determines the
        type of averaging performed on the data.

        macro:
            Average over classes (does not take imbalance into account).
        micro:
            Average over instances (takes imbalance into account).
            This implies that ``precision == recall == f1``
        weighted:
            Average weighted by support (takes imbalance into account).
            Can result in f1 score that is not between precision and recall.

    Returns
    -------
    fbeta_score : float
        fbeta_score of the positive class in binary classification or weighted
        average of the fbeta_score of each class for the multiclass task.

    References
    ----------
    R. Baeza-Yates and B. Ribeiro-Neto (2011). Modern Information Retrieval.
    Addison Wesley, pp. 327-328.

    http://en.wikipedia.org/wiki/F1_score

    """
    _, _, f, _ = precision_recall_fscore_support(y_true, y_pred,
                                                 beta=beta,
                                                 labels=labels,
                                                 pos_label=pos_label,
                                                 average=average)
    return f


def f1_score(y_true, y_pred, labels=None, pos_label=1, average='weighted'):
    """Compute f1 score

    The F1 score can be interpreted as a weighted average of the precision
    and recall, where an F1 score reaches its best value at 1 and worst
    score at 0. The relative contribution of precision and recall to the f1
    score are equal. The formular for the F_1 score is::

        F_1 = 2 * (precision * recall) / (precision + recall)

    See: http://en.wikipedia.org/wiki/F1_score

    In the multi-class case, this is the weighted average of the f1-score of
    each class.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Predicted targets

    labels : array
        Integer array of labels

    pos_label : int
        In the binary classification case, give the label of the positive
        class (default is 1). Everything else but 'pos_label'
        is considered to belong to the negative class.
        Set to None in the case of multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted'(default)]
        In the multiclass classification case, this determines the
        type of averaging performed on the data.

        macro:
            Average over classes (does not take imbalance into account).
        micro:
            Average over instances (takes imbalance into account).
            This implies that ``precision == recall == f1``
        weighted:
            Average weighted by support (takes imbalance into account).
            Can result in f1 score that is not between precision and recall.

    Returns
    -------
    f1_score : float
        f1_score of the positive class in binary classification or weighted
        average of the f1_scores of each class for the multiclass task

    References
    ----------
    http://en.wikipedia.org/wiki/F1_score

    """
    return fbeta_score(y_true, y_pred, 1, labels=labels,
                       pos_label=pos_label, average=average)


def precision_recall_fscore_support(y_true, y_pred, beta=1.0, labels=None,
                                    pos_label=1, average=None):
    """Compute precisions, recalls, f-measures and support for each class

    The precision is the ratio :math:`tp / (tp + fp)` where tp is the number of
    true positives and fp the number of false positives. The precision is
    intuitively the ability of the classifier not to label as positive a sample
    that is negative.

    The recall is the ratio :math:`tp / (tp + fn)` where tp is the number of
    true positives and fn the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The F_beta score can be interpreted as a weighted harmonic mean of
    the precision and recall, where an F_beta score reaches its best
    value at 1 and worst score at 0.

    The F_beta score weights recall beta as much as precision. beta = 1.0 means
    recall and precsion are equally important.

    The support is the number of occurrences of each class in y_true.

    If pos_label is None, this function returns the average precision, recall
    and f-measure if `average` is one of 'micro', 'macro', 'weighted'.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Predicted targets

    beta : float, 1.0 by default
        The strength of recall versus precision in the f-score.

    labels : array
        Integer array of labels

    pos_label : int
        In the binary classification case, give the label of the positive
        class (default is 1). Everything else but 'pos_label'
        is considered to belong to the negative class.
        Set to None in the case of multiclass classification.

    average : string, [None, 'micro', 'macro', 'weighted'(default)]
        In the multiclass classification case, this determines the
        type of averaging performed on the data.

        macro:
            Average over classes (does not take imbalance into account).
        micro:
            Average over instances (takes imbalance into account).
            This implies that ``precision == recall == f1``
        weighted:
            Average weighted by support (takes imbalance into account).
            Can result in f1 score that is not between precision and recall.

    Returns
    -------
    precision: array, shape = [n_unique_labels], dtype = np.double
    recall: array, shape = [n_unique_labels], dtype = np.double
    f1_score: array, shape = [n_unique_labels], dtype = np.double
    support: array, shape = [n_unique_labels], dtype = np.long

    References
    ----------
    http://en.wikipedia.org/wiki/Precision_and_recall

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

    elif n_labels == 2:
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


def matthews_corrcoef(y_true, y_pred):
    """Returns matthew's correlation coefficient for binary classes

    The Matthews correlation coefficient is used in machine learning as a
    measure of the quality of binary (two-class) classifications. It takes
    into account true and false positives and negatives and is generally
    regarded as a balanced measure which can be used even if the classes are
    of very different sizes. The MCC is in essence a correlation coefficient
    value between -1 and +1. A coefficient of +1 represents a perfect
    prediction, 0 an average random prediction and -1 an inverse prediction.
    The statistic is also known as the phi coefficient. [source: Wikipedia]

    Only in the binary case does this relate to information about true and
    false positives and negatives. See references below.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        true targets

    y_pred : array, shape = [n_samples]
        estimated targets

    Returns
    -------
    mcc : float
        matthew's correlation coefficient (+1 represents a perfect prediction,
        0 an average random prediction and -1 and inverse prediction).

    References
    ----------
    http://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    http://dx.doi.org/10.1093/bioinformatics/16.5.412

    """
    mcc = np.corrcoef(y_true, y_pred)[0, 1]
    if np.isnan(mcc):
        return 0.
    else:
        return mcc


def classification_report(y_true, y_pred, labels=None, target_names=None):
    """Build a text report showing the main classification metrics

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets

    y_pred : array, shape = [n_samples]
        Estimated targets

    labels : array, shape = [n_labels]
        Optional list of label indices to include in the report

    target_names : list of strings
        Optional display names matching the labels (same order)

    Returns
    -------
    report : string
        Text summary of the precision, recall, f1-score for each class

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


def precision_recall_curve(y_true, probas_pred):
    """Compute precision-recall pairs for different probability thresholds

    Note: this implementation is restricted to the binary classification task.

    The precision is the ratio :math:`tp / (tp + fp)` where tp is the number of
    true positives and fp the number of false positives. The precision is
    intuitively the ability of the classifier not to label as positive a sample
    that is negative.

    The recall is the ratio :math:`tp / (tp + fn)` where tp is the number of
    true positives and fn the number of false negatives. The recall is
    intuitively the ability of the classifier to find all the positive samples.

    The last precision and recall values are 1. and 0. respectively and do not
    have a corresponding threshold.  This ensures that the graph starts on the
    x axis.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True targets of binary classification in range {-1, 1} or {0, 1}

    probas_pred : array, shape = [n_samples]
        Estimated probabilities

    Returns
    -------
    precision : array, shape = [n + 1]
        Precision values

    recall : array, shape = [n + 1]
        Recall values

    thresholds : array, shape = [n]
        Thresholds on y_score used to compute precision and recall

    """
    y_true = np.ravel(y_true)
    labels = np.unique(y_true)
    if np.all(labels == np.array([-1, 1])):
        # convert {-1, 1} to boolean {0, 1} repr
        y_true = y_true.copy()
        y_true[y_true == -1] = 0
    elif not np.all(labels == np.array([0, 1])):
        raise ValueError("y_true contains non binary labels: %r" % labels)

    probas_pred = np.ravel(probas_pred)
    thresholds = np.sort(np.unique(probas_pred))
    n_thresholds = thresholds.size + 1
    precision = np.empty(n_thresholds)
    recall = np.empty(n_thresholds)
    for i, t in enumerate(thresholds):
        y_pred = (probas_pred >= t).astype(np.int)
        p, r, _, _ = precision_recall_fscore_support(y_true, y_pred)
        precision[i] = p[1]
        recall[i] = r[1]
    precision[-1] = 1.0
    recall[-1] = 0.0
    return precision, recall, thresholds


def explained_variance_score(y_true, y_pred):
    """Explained variance regression score function

    Best possible score is 1.0, lower values are worse.

    Note: the explained variance is not a symmetric function.

    return the explained variance

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

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
    y_true : array-like

    y_pred : array-like

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
    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    if len(y_true) == 1:
        raise ValueError("r2_score can only be computed given more than one"
                " sample.")
    numerator = ((y_true - y_pred) ** 2).sum()
    denominator = ((y_true - y_true.mean()) ** 2).sum()
    if denominator == 0.0:
        if numerator == 0.0:
            return 1.0
        else:
            # arbitary set to zero to avoid -inf scores, having a constant
            # y_true is not interesting for scoring a regression anyway
            return 0.0
    return 1 - numerator / denominator


def zero_one_score(y_true, y_pred):
    """Zero-one classification score (accuracy)

    Positive integer (number of good classifications).
    The best performance is 1.

    Return the fraction of correct predictions in y_pred.

    Parameters
    ----------
    y_true : array-like, shape = n_samples
        Gold standard labels.

    y_pred : array-like, shape = n_samples
        Predicted labels, as returned by a classifier.

    Returns
    -------
    score : float

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.mean(y_pred == y_true)


###############################################################################
# Loss functions

def zero_one(y_true, y_pred):
    """Zero-One classification loss

    Positive integer (number of misclassifications). The best performance
    is 0.

    Return the number of errors

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

    Returns
    -------
    loss : float

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.sum(y_pred != y_true)


def mean_squared_error(y_true, y_pred):
    """Mean squared error regression loss

    Return a a positive floating point value (the best value is 0.0).

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

    Returns
    -------
    loss : float
    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.mean((y_pred - y_true) ** 2)


@deprecated("""Incorrectly returns the cumulated error: use mean_squared_error
            instead; to be removed in v0.13""")
def mean_square_error(y_true, y_pred):
    """Cumulated square error regression loss

    Positive floating point value: the best value is 0.0.

    return the mean square error

    Parameters
    ----------
    y_true : array-like

    y_pred : array-like

    Returns
    -------
    loss : float

    """
    y_true, y_pred = check_arrays(y_true, y_pred)
    return np.linalg.norm(y_pred - y_true) ** 2


def hinge_loss(y_true, pred_decision, pos_label=1, neg_label=-1):
    """
    Cumulated hinge loss (non-regularized).

    Assuming labels in y_true are encoded with +1 and -1,
    when a prediction mistake is made, margin = y_true * pred_decision
    is always negative (since the signs disagree), therefore 1 - margin
    is always greater than 1. The cumulated hinge loss therefore
    upperbounds the number of mistakes made by the classifier.

    Parameters
    ----------
    y_true : array, shape = [n_samples]
        True target (integers)

    pred_decision : array, shape = [n_samples] or [n_samples, n_classes]
        Predicted decisions, as output by decision_function (floats)

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
