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

from itertools import combinations

import numpy as np

from ..base import is_classifier, is_regressor
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

    sample_weight : array-like of shape (n_samples,), default=None
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
        if average_weight is not None:
            # Scores with 0 weights are forced to be 0, preventing the average
            # score from being affected by 0-weighted NaN elements.
            average_weight = np.asarray(average_weight)
            score[average_weight == 0] = 0
        return np.average(score, weights=average_weight)
    else:
        return score


def _average_multiclass_ovo_score(binary_metric, y_true, y_score,
                                  average='macro'):
    """Average one-versus-one scores for multiclass classification.

    Uses the binary metric for one-vs-one multiclass classification,
    where the score is computed according to the Hand & Till (2001) algorithm.

    Parameters
    ----------
    binary_metric : callable
        The binary metric function to use that accepts the following as input
            y_true_target : array, shape = [n_samples_target]
                Some sub-array of y_true for a pair of classes designated
                positive and negative in the one-vs-one scheme.
            y_score_target : array, shape = [n_samples_target]
                Scores corresponding to the probability estimates
                of a sample belonging to the designated positive class label

    y_true : array-like of shape (n_samples,)
        True multiclass labels.

    y_score : array-like of shape (n_samples, n_classes)
        Target scores corresponding to probability estimates of a sample
        belonging to a particular class

    average : {'macro', 'weighted'}, default='macro'
        Determines the type of averaging performed on the pairwise binary
        metric scores
        ``'macro'``:
            Calculate metrics for each label, and find their unweighted
            mean. This does not take label imbalance into account. Classes
            are assumed to be uniformly distributed.
        ``'weighted'``:
            Calculate metrics for each label, taking into account the
            prevalence of the classes.

    Returns
    -------
    score : float
        Average of the pairwise binary metric scores
    """
    check_consistent_length(y_true, y_score)

    y_true_unique = np.unique(y_true)
    n_classes = y_true_unique.shape[0]
    n_pairs = n_classes * (n_classes - 1) // 2
    pair_scores = np.empty(n_pairs)

    is_weighted = average == "weighted"
    prevalence = np.empty(n_pairs) if is_weighted else None

    # Compute scores treating a as positive class and b as negative class,
    # then b as positive class and a as negative class
    for ix, (a, b) in enumerate(combinations(y_true_unique, 2)):
        a_mask = y_true == a
        b_mask = y_true == b
        ab_mask = np.logical_or(a_mask, b_mask)

        if is_weighted:
            prevalence[ix] = np.average(ab_mask)

        a_true = a_mask[ab_mask]
        b_true = b_mask[ab_mask]

        a_true_score = binary_metric(a_true, y_score[ab_mask, a])
        b_true_score = binary_metric(b_true, y_score[ab_mask, b])
        pair_scores[ix] = (a_true_score + b_true_score) / 2

    return np.average(pair_scores, weights=prevalence)


def _check_classifier_response_method(estimator, response_method):
    """Return prediction method from the response_method

    Parameters
    ----------
    estimator: object
        Classifier to check

    response_method: {'auto', 'predict_proba', 'decision_function', 'predict'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next and :term:`predict` last.

    Returns
    -------
    prediction_method: callable
        prediction method of estimator
    """

    possible_response_methods = (
        "predict", "predict_proba", "decision_function", "auto"
    )
    if response_method not in possible_response_methods:
        raise ValueError(
            f"response_method must be one of "
            f"{','.join(possible_response_methods)}."
        )

    error_msg = "response method {} is not defined in {}"
    if response_method != "auto":
        prediction_method = getattr(estimator, response_method, None)
        if prediction_method is None:
            raise ValueError(
                error_msg.format(response_method, estimator.__class__.__name__)
            )
    else:
        predict_proba = getattr(estimator, 'predict_proba', None)
        decision_function = getattr(estimator, 'decision_function', None)
        predict = getattr(estimator, 'predict', None)
        prediction_method = predict_proba or decision_function or predict
        if prediction_method is None:
            raise ValueError(
                error_msg.format(
                    "decision_function, predict_proba or predict",
                    estimator.__class__.__name__
                )
            )

    return prediction_method


def _get_response(
    estimator,
    X,
    y_true,
    response_method,
    pos_label=None,
    support_multi_class=False,
):
    """Return response and positive label.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    y_true : array-like of shape (n_samples,)
        The true label.

    response_method: {'auto', 'predict_proba', 'decision_function', 'predict'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next and :term:`predict` last.

    pos_label : str or int, default=None
        The class considered as the positive class when computing
        the metrics. By default, `estimators.classes_[1]` is
        considered as the positive class.

    support_multi_class : bool, default=False
        ...

    Returns
    -------
    y_pred : ndarray of shape (n_samples,)
        Target scores calculated from the provided response_method
        and pos_label.

    pos_label : str or int
        The class considered as the positive class when computing
        the metrics.
    """
    if is_regressor(estimator):
        if response_method not in ("predict", "auto"):
            raise ValueError(
                f"{estimator.__class__.__name__} should be a classifier"
            )
        return estimator.predict(X), None
    else:
        y_type = type_of_target(y_true)
        classes = estimator.classes_
        prediction_method = _check_classifier_response_method(
            estimator, response_method
        )
        y_pred = prediction_method(X)

        if pos_label is not None and pos_label not in classes:
            raise ValueError(
                f"The class provided by 'pos_label' is unknown. Got "
                f"{pos_label} instead of one of {classes}."
            )

        if prediction_method.__name__ == "predict_proba":
            if y_type == "binary":
                pos_label = pos_label if pos_label is not None else classes[-1]
                if y_pred.shape[1] == 2:
                    col_idx = np.flatnonzero(classes == pos_label)[0]
                    y_pred = y_pred[:, col_idx]
                else:
                    err_msg = (
                        f"Got predict_proba of shape {y_pred.shape}, but need "
                        f"classifier with two classes"
                    )
                    if support_multi_class and y_pred.shape[1] == 1:
                        raise ValueError(err_msg)
                    elif not support_multi_class:
                        raise ValueError(err_msg)
        elif prediction_method.__name__ == "decision_function":
            if y_type == "binary":
                pos_label = pos_label if pos_label is not None else classes[-1]
                if pos_label == classes[0]:
                    y_pred *= -1

    return y_pred, pos_label
