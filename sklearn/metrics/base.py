"""
Common code for all metrics

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


import numpy as np

from ..utils import check_array, check_consistent_length
from ..utils.multiclass import type_of_target


def _average_binary_score(binary_metric, y_true, y_score, average,
                          sample_weight=None):
    """Average a binary metric for multilabel classification

    Parameters
    ----------
    y_true : array, shape = [n_samples] or [n_samples, n_classes]
        True binary labels in binary label indicators.

    y_score : array, shape = [n_samples] or [n_samples, n_classes]
        Target scores, can either be probability estimates of the positive
        class, confidence values, or binary decisions.

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

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    binary_metric : callable, returns shape [n_classes]
        The binary metric function to use.

    Returns
    -------
    score : float or array of shape [n_classes]
        If not ``None``, average the score, else return the score for each
        classes.

    """
    average_options = (None, 'micro', 'macro', 'weighted', 'samples')
    if average not in average_options:
        raise ValueError('average has to be one of {0}'
                         ''.format(average_options))

    y_type = type_of_target(y_true)
    if y_type not in ("binary", "multilabel-indicator"):
        raise ValueError("{0} format is not supported".format(y_type))

    if y_type == "binary":
        return binary_metric(y_true, y_score, sample_weight=sample_weight)

    check_consistent_length(y_true, y_score, sample_weight)
    y_true = check_array(y_true)
    y_score = check_array(y_score)

    not_average_axis = 1
    score_weight = sample_weight
    average_weight = None

    if average == "micro":
        if score_weight is not None:
            score_weight = np.repeat(score_weight, y_true.shape[1])
        y_true = y_true.ravel()
        y_score = y_score.ravel()

    elif average == 'weighted':
        if score_weight is not None:
            average_weight = np.sum(np.multiply(
                y_true, np.reshape(score_weight, (-1, 1))), axis=0)
        else:
            average_weight = np.sum(y_true, axis=0)
        if np.isclose(average_weight.sum(), 0.0):
            return 0

    elif average == 'samples':
        # swap average_weight <-> score_weight
        average_weight = score_weight
        score_weight = None
        not_average_axis = 0

    if y_true.ndim == 1:
        y_true = y_true.reshape((-1, 1))

    if y_score.ndim == 1:
        y_score = y_score.reshape((-1, 1))

    n_classes = y_score.shape[not_average_axis]
    score = np.zeros((n_classes,))
    for c in range(n_classes):
        y_true_c = y_true.take([c], axis=not_average_axis).ravel()
        y_score_c = y_score.take([c], axis=not_average_axis).ravel()
        score[c] = binary_metric(y_true_c, y_score_c,
                                 sample_weight=score_weight)

    # Average the results
    if average is not None:
        return np.average(score, weights=average_weight)
    else:
        return score
