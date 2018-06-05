from __future__ import division

import warnings
import numpy as np

from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

from sklearn.preprocessing import LabelBinarizer, label_binarize
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import assert_all_finite
from sklearn.utils import check_array
from sklearn.utils import check_consistent_length
from sklearn.utils import column_or_1d
from sklearn.utils.multiclass import unique_labels
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import _num_samples
from sklearn.utils.sparsefuncs import count_nonzero
from sklearn.exceptions import UndefinedMetricWarning

from sklearn.metrics.classification import _check_targets, _weighted_sum, confusion_matrix, \
        multilabel_confusion_matrix, _prf_divide



@profile
def precision_recall_fscore_support(y_true, y_pred, beta=1.0, labels=None,
                                    pos_label=1, average=None,
                                    warn_for=('precision', 'recall',
                                              'f-score'),
                                    sample_weight=None):
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
    ``beta``. ``beta == 1.0`` means recall and precision are equally important.

    The support is the number of occurrences of each class in ``y_true``.

    If ``pos_label is None`` and in binary classification, this function
    returns the average precision, recall and F-measure if ``average``
    is one of ``'micro'``, ``'macro'``, ``'weighted'`` or ``'samples'``.

    Read more in the :ref:`User Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    y_true : 1d array-like, or label indicator array / sparse matrix
        Ground truth (correct) target values.

    y_pred : 1d array-like, or label indicator array / sparse matrix
        Estimated targets as returned by a classifier.

    beta : float, 1.0 by default
        The strength of recall versus precision in the F-score.

    labels : list, optional
        The set of labels to include when ``average != 'binary'``, and their
        order if ``average is None``. Labels present in the data can be
        excluded, for example to calculate a multiclass average ignoring a
        majority negative class, while labels not present in the data will
        result in 0 components in a macro average. For multilabel targets,
        labels are column indices. By default, all labels in ``y_true`` and
        ``y_pred`` are used in sorted order.

    pos_label : str or int, 1 by default
        The class to report if ``average='binary'`` and the data is binary.
        If the data are multiclass or multilabel, this will be ignored;
        setting ``labels=[pos_label]`` and ``average != 'binary'`` will report
        scores for that label only.

    average : string, [None (default), 'binary', 'micro', 'macro', 'samples', \
                       'weighted']
        If ``None``, the scores for each class are returned. Otherwise, this
        determines the type of averaging performed on the data:

        ``'binary'``:
            Only report results for the class specified by ``pos_label``.
            This is applicable only if targets (``y_{true,pred}``) are binary.
        ``'micro'``:
            Calculate metrics globally by counting the total true positives,
            false negatives and false positives.
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean.  This does not take label imbalance into account.
        ``'weighted'``:
            Calculate metrics for each label, and find their average, weighted
            by support (the number of true instances for each label). This
            alters 'macro' to account for label imbalance; it can result in an
            F-score that is not between precision and recall.
        ``'samples'``:
            Calculate metrics for each instance, and find their average (only
            meaningful for multilabel classification where this differs from
            :func:`accuracy_score`).

    warn_for : tuple or set, for internal use
        This determines which warnings will be made in the case that this
        function is being used to return only one of its metrics.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    precision : float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    recall : float (if average is not None) or array of float, , shape =\
        [n_unique_labels]

    fbeta_score : float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    support : int (if average is not None) or array of int, shape =\
        [n_unique_labels]
        The number of occurrences of each label in ``y_true``.

    References
    ----------
    .. [1] `Wikipedia entry for the Precision and recall
           <https://en.wikipedia.org/wiki/Precision_and_recall>`_

    .. [2] `Wikipedia entry for the F1-score
           <https://en.wikipedia.org/wiki/F1_score>`_

    .. [3] `Discriminative Methods for Multi-labeled Classification Advances
           in Knowledge Discovery and Data Mining (2004), pp. 22-30 by Shantanu
           Godbole, Sunita Sarawagi
           <http://www.godbole.net/shantanu/pubs/multilabelsvm-pakdd04.pdf>`_

    Examples
    --------
    >>> from sklearn.metrics import precision_recall_fscore_support
    >>> y_true = np.array(['cat', 'dog', 'pig', 'cat', 'dog', 'pig'])
    >>> y_pred = np.array(['cat', 'pig', 'dog', 'cat', 'cat', 'dog'])
    >>> precision_recall_fscore_support(y_true, y_pred, average='macro')
    ... # doctest: +ELLIPSIS
    (0.22..., 0.33..., 0.26..., None)
    >>> precision_recall_fscore_support(y_true, y_pred, average='micro')
    ... # doctest: +ELLIPSIS
    (0.33..., 0.33..., 0.33..., None)
    >>> precision_recall_fscore_support(y_true, y_pred, average='weighted')
    ... # doctest: +ELLIPSIS
    (0.22..., 0.33..., 0.26..., None)

    It is possible to compute per-label precisions, recalls, F1-scores and
    supports instead of averaging:

    >>> precision_recall_fscore_support(y_true, y_pred, average=None,
    ... labels=['pig', 'dog', 'cat'])
    ... # doctest: +ELLIPSIS,+NORMALIZE_WHITESPACE
    (array([0.        , 0.        , 0.66...]),
     array([0., 0., 1.]), array([0. , 0. , 0.8]),
     array([2, 2, 2]))

    """
    average_options = (None, 'micro', 'macro', 'weighted', 'samples')
    if average not in average_options and average != 'binary':
        raise ValueError('average has to be one of ' +
                         str(average_options))
    if beta <= 0:
        raise ValueError("beta should be >0 in the F-beta score")

    y_type, y_true, y_pred = _check_targets(y_true, y_pred)
    check_consistent_length(y_true, y_pred, sample_weight)
    present_labels = unique_labels(y_true, y_pred)

    if average == 'binary':
        if y_type == 'binary':
            if pos_label not in present_labels:
                if len(present_labels) < 2:
                    # Only negative labels
                    return (0., 0., 0., 0)
                else:
                    raise ValueError("pos_label=%r is not a valid label: %r" %
                                     (pos_label, present_labels))
            labels = [pos_label]
        else:
            raise ValueError("Target is %s but average='binary'. Please "
                             "choose another average setting." % y_type)
    elif pos_label not in (None, 1):
        warnings.warn("Note that pos_label (set to %r) is ignored when "
                      "average != 'binary' (got %r). You may use "
                      "labels=[pos_label] to specify a single positive class."
                      % (pos_label, average), UserWarning)

    if labels is None:
        labels = present_labels
        n_labels = None
    else:
        n_labels = len(labels)
        labels = np.hstack([labels, np.setdiff1d(present_labels, labels,
                                                 assume_unique=True)])

    # Calculate tp_sum, pred_sum, true_sum ###

    if y_type.startswith('multilabel'):
        sum_axis = 1 if average == 'samples' else 0

        # All labels are index integers for multilabel.
        # Select labels:
        if not np.all(labels == present_labels):
            if np.max(labels) > np.max(present_labels):
                raise ValueError('All labels must be in [0, n labels). '
                                 'Got %d > %d' %
                                 (np.max(labels), np.max(present_labels)))
            if np.min(labels) < 0:
                raise ValueError('All labels must be in [0, n labels). '
                                 'Got %d < 0' % np.min(labels))

        if n_labels is not None:
            y_true = y_true[:, labels[:n_labels]]
            y_pred = y_pred[:, labels[:n_labels]]

        # calculate weighted counts
        true_and_pred = y_true.multiply(y_pred)
        tp_sum = count_nonzero(true_and_pred, axis=sum_axis,
                               sample_weight=sample_weight)
        pred_sum = count_nonzero(y_pred, axis=sum_axis,
                                 sample_weight=sample_weight)
        true_sum = count_nonzero(y_true, axis=sum_axis,
                                 sample_weight=sample_weight)

    elif average == 'samples':
        raise ValueError("Sample-based precision, recall, fscore is "
                         "not meaningful outside multilabel "
                         "classification. See the accuracy_score instead.")
    else:
        le = LabelEncoder()
        le.fit(labels)
        y_true = le.transform(y_true)
        y_pred = le.transform(y_pred)
        sorted_labels = le.classes_

        # labels are now from 0 to len(labels) - 1 -> use bincount
        tp = y_true == y_pred
        tp_bins = y_true[tp]
        if sample_weight is not None:
            tp_bins_weights = np.asarray(sample_weight)[tp]
        else:
            tp_bins_weights = None

        if len(tp_bins):
            tp_sum = np.bincount(tp_bins, weights=tp_bins_weights,
                              minlength=len(labels))
        else:
            # Pathological case
            true_sum = pred_sum = tp_sum = np.zeros(len(labels))
        if len(y_pred):
            pred_sum = np.bincount(y_pred, weights=sample_weight,
                                minlength=len(labels))
        if len(y_true):
            true_sum = np.bincount(y_true, weights=sample_weight,
                                minlength=len(labels))

        # Retain only selected labels
        indices = np.searchsorted(sorted_labels, labels[:n_labels])
        tp_sum = tp_sum[indices]
        true_sum = true_sum[indices]
        pred_sum = pred_sum[indices]

    if average == 'micro':
        tp_sum = np.array([tp_sum.sum()])
        pred_sum = np.array([pred_sum.sum()])
        true_sum = np.array([true_sum.sum()])

    # Finally, we have all our sufficient statistics. Divide! #

    beta2 = beta ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        # Divide, and on zero-division, set scores to 0 and warn:

        # Oddly, we may get an "invalid" rather than a "divide" error
        # here.
        precision = _prf_divide(tp_sum, pred_sum,
                                'precision', 'predicted', average, warn_for)
        recall = _prf_divide(tp_sum, true_sum,
                             'recall', 'true', average, warn_for)
        # Don't need to warn for F: either P or R warned, or tp == 0 where pos
        # and true are nonzero, in which case, F is well-defined and zero
        f_score = ((1 + beta2) * precision * recall /
                   (beta2 * precision + recall))
        f_score[tp_sum == 0] = 0.0

    # Average the results

    if average == 'weighted':
        weights = true_sum
        if weights.sum() == 0:
            return 0, 0, 0, None
    elif average == 'samples':
        weights = sample_weight
    else:
        weights = None

    if average is not None:
        assert average != 'binary' or len(precision) == 1
        precision = np.average(precision, weights=weights)
        recall = np.average(recall, weights=weights)
        f_score = np.average(f_score, weights=weights)
        true_sum = None  # return no support

    return precision, recall, f_score, true_sum


@profile
def precision_recall_fscore_support_with_multilabel_confusion_matrix(
                                    y_true, y_pred,
                                    beta=1.0, labels=None,
                                    pos_label=1, average=None,
                                    warn_for=('precision', 'recall',
                                              'f-score'),
                                    sample_weight=None):
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
    ``beta``. ``beta == 1.0`` means recall and precision are equally important.

    The support is the number of occurrences of each class in ``y_true``.

    If ``pos_label is None`` and in binary classification, this function
    returns the average precision, recall and F-measure if ``average``
    is one of ``'micro'``, ``'macro'``, ``'weighted'`` or ``'samples'``.

    Read more in the :ref:`User Guide <precision_recall_f_measure_metrics>`.

    Parameters
    ----------
    y_true : 1d array-like, or label indicator array / sparse matrix
        Ground truth (correct) target values.

    y_pred : 1d array-like, or label indicator array / sparse matrix
        Estimated targets as returned by a classifier.

    beta : float, 1.0 by default
        The strength of recall versus precision in the F-score.

    labels : list, optional
        The set of labels to include when ``average != 'binary'``, and their
        order if ``average is None``. Labels present in the data can be
        excluded, for example to calculate a multiclass average ignoring a
        majority negative class, while labels not present in the data will
        result in 0 components in a macro average. For multilabel targets,
        labels are column indices. By default, all labels in ``y_true`` and
        ``y_pred`` are used in sorted order.

    pos_label : str or int, 1 by default
        The class to report if ``average='binary'`` and the data is binary.
        If the data are multiclass or multilabel, this will be ignored;
        setting ``labels=[pos_label]`` and ``average != 'binary'`` will report
        scores for that label only.

    average : string, [None (default), 'binary', 'micro', 'macro', 'samples', \
                       'weighted']
        If ``None``, the scores for each class are returned. Otherwise, this
        determines the type of averaging performed on the data:

        ``'binary'``:
            Only report results for the class specified by ``pos_label``.
            This is applicable only if targets (``y_{true,pred}``) are binary.
        ``'micro'``:
            Calculate metrics globally by counting the total true positives,
            false negatives and false positives.
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean.  This does not take label imbalance into account.
        ``'weighted'``:
            Calculate metrics for each label, and find their average, weighted
            by support (the number of true instances for each label). This
            alters 'macro' to account for label imbalance; it can result in an
            F-score that is not between precision and recall.
        ``'samples'``:
            Calculate metrics for each instance, and find their average (only
            meaningful for multilabel classification where this differs from
            :func:`accuracy_score`).

    warn_for : tuple or set, for internal use
        This determines which warnings will be made in the case that this
        function is being used to return only one of its metrics.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    precision : float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    recall : float (if average is not None) or array of float, , shape =\
        [n_unique_labels]

    fbeta_score : float (if average is not None) or array of float, shape =\
        [n_unique_labels]

    support : int (if average is not None) or array of int, shape =\
        [n_unique_labels]
        The number of occurrences of each label in ``y_true``.

    See also
    --------
    precision_score, recall_score, f1_score, fbeta_score,
    jaccard_similarity_score, multilabel_confusion_matrix

    References
    ----------
    .. [1] `Wikipedia entry for the Precision and recall
           <https://en.wikipedia.org/wiki/Precision_and_recall>`_

    .. [2] `Wikipedia entry for the F1-score
           <https://en.wikipedia.org/wiki/F1_score>`_

    .. [3] `Discriminative Methods for Multi-labeled Classification Advances
           in Knowledge Discovery and Data Mining (2004), pp. 22-30 by Shantanu
           Godbole, Sunita Sarawagi
           <http://www.godbole.net/shantanu/pubs/multilabelsvm-pakdd04.pdf>`_

    """
    average_options = (None, 'micro', 'macro', 'weighted', 'samples')
    if average not in average_options and average != 'binary':
        raise ValueError('average has to be one of ' +
                         str(average_options))
    if beta <= 0:
        raise ValueError("beta should be >0 in the F-beta score")

    y_type, y_true, y_pred = _check_targets(y_true, y_pred)
    present_labels = unique_labels(y_true, y_pred)

    if average == 'binary':
        if y_type == 'binary':
            if pos_label not in present_labels:
                if len(present_labels) < 2:
                    # Only negative labels
                    return (0., 0., 0., 0)
                else:
                    raise ValueError("pos_label=%r is not a valid label: %r" %
                                     (pos_label, present_labels))
            labels = [pos_label]
        else:
            raise ValueError("Target is %s but average='binary'. Please "
                             "choose another average setting." % y_type)
    elif pos_label not in (None, 1):
        warnings.warn("Note that pos_label (set to %r) is ignored when "
                      "average != 'binary' (got %r). You may use "
                      "labels=[pos_label] to specify a single positive class."
                      % (pos_label, average), UserWarning)

    # Calculate sufficient statistics: tp_sum, pred_sum, true_sum ###
    C = multilabel_confusion_matrix(y_true, y_pred,
                                    sample_weight=sample_weight,
                                    labels=labels,
                                    samplewise=average == 'samples')

    if average == 'micro':
        C = C.sum(axis=0).reshape(1, 2, 2)

    tp_sum = C[:, 1, 1]
    pred_sum = C[:, 1, 1] + C[:, 0, 1]
    true_sum = C[:, 1, 1] + C[:, 1, 0]

    # Divide!

    beta2 = beta ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        # Divide, and on zero-division, set scores to 0 and warn:

        # Oddly, we may get an "invalid" rather than a "divide" error
        # here.
        precision = _prf_divide(tp_sum, pred_sum,
                                'precision', 'predicted', average, warn_for)
        recall = _prf_divide(tp_sum, true_sum,
                             'recall', 'true', average, warn_for)
        # Don't need to warn for F: either P or R warned, or tp == 0 where pos
        # and true are nonzero, in which case, F is well-defined and zero
        f_score = ((1 + beta2) * precision * recall /
                   (beta2 * precision + recall))
        f_score[tp_sum == 0] = 0.0

    # Average the results

    if average == 'weighted':
        weights = true_sum
        if weights.sum() == 0:
            return 0, 0, 0, None
    elif average == 'samples':
        weights = sample_weight
    else:
        weights = None

    if average is not None:
        assert average != 'binary' or len(precision) == 1
        precision = np.average(precision, weights=weights)
        recall = np.average(recall, weights=weights)
        f_score = np.average(f_score, weights=weights)
        true_sum = None  # return no support

    return precision, recall, f_score, true_sum


if __name__ == '__main__':
    from scipy.stats import bernoulli

    n_samples = 30000
    n_labels = 2000

    y_true = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2,
                           size=(n_samples, n_labels))
    y_pred = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2,
                           size=(n_samples, n_labels))

    #precision_recall_fscore_support_with_multilabel_confusion_matrix(y_true,
    #                                                                 y_pred)
    #precision_recall_fscore_support(y_true,
    #                                                                 y_pred)

    multilabel_confusion_matrix(y_true, y_pred)


