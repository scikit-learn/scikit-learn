"""
The :mod:`sklearn.metrics.scorer` submodule implements a flexible
interface for model selection and evaluation using
arbitrary score functions.

A scorer object is a callable that can be passed to
:class:`sklearn.grid_search.GridSearchCV` or
:func:`sklearn.cross_validation.cross_val_score` as the ``scoring`` parameter,
to specify how a model should be evaluated.

The signature of the call is ``(estimator, X, y)`` where ``estimator``
is the model to be evaluated, ``X`` is the test data and ``y`` is the
ground truth labeling (or ``None`` in the case of unsupervised models).
"""

# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
# License: Simplified BSD

from abc import ABCMeta, abstractmethod
from warnings import warn

import numpy as np

from . import (r2_score, mean_absolute_error, mean_squared_error,
               accuracy_score, f1_score, roc_auc_score, average_precision_score,
               precision_score, recall_score, log_loss)
from .cluster import adjusted_rand_score
from ..utils.multiclass import type_of_target
from ..externals import six
from  ..base import is_classifier


class _Scorer(object):

    def __init__(self, score_func, greater_is_better=True, needs_proba=False,
                 needs_threshold=False, kwargs={}):
        self.score_func = score_func
        self.greater_is_better = greater_is_better
        self.needs_proba = needs_proba
        self.needs_threshold = needs_threshold
        self.kwargs = kwargs

    def __call__(self, estimator, X, y):
        return _evaluate_scorers(estimator, X, y, [self])[0]


def _evaluate_scorers(estimator, X, y, scorers):
    has_pb = hasattr(estimator, "predict_proba")
    has_df = hasattr(estimator, "decision_function")
    _is_classifier = is_classifier(estimator)
    _type_of_y = type_of_target(y)

    # Make a first pass through scorers to determine if we need
    # predict_proba or decision_function.
    needs_proba = False
    needs_df = False
    for scorer in scorers:
        if scorer.needs_proba:
            if not has_pb:
                raise ValueError("%s needs probabilities but predict_proba is"
                                 "not available in %s." % (scorer, estimator))
            needs_proba = True

        elif scorer.needs_threshold:
            if has_pb:
                # We choose predict_proba first because its interface
                # is more consistent across the project.
                needs_proba = True
                continue

            if _is_classifier and not has_df:
                raise ValueError("%s needs continuous outputs but neither"
                                 "predict_proba nor decision_function "
                                 "are available in %s." % (scorer, estimator))

            if _is_classifier:
                needs_df = True

    # Compute predict_proba if needed.
    y_proba = None
    y_pred = None
    if needs_proba:
        try:
            y_proba = estimator.predict_proba(X)

            y_pred = estimator.classes_[y_proba.argmax(axis=1)]

            if _type_of_y == "binary":
                y_proba = y_proba[:, 1]

        except (NotImplementedError, AttributeError):
            # SVC has predict_proba but it may raise NotImplementedError
            # if probabilities are not enabled.
            needs_proba = False
            needs_df = True

    # Compute decision_function.
    df = None
    if needs_df:
        df = estimator.decision_function(X)

        if len(df.shape) == 2 and df.shape[1] >= 2:
            y_pred = estimator.classes_[df.argmax(axis=1)]
        else:
            y_pred = estimator.classes_[(df >= 0).astype(int)]

    # Compute y_pred if needed.
    if y_pred is None:
        y_pred = estimator.predict(X)

    # Compute scores.
    scores = []
    for scorer in scorers:
        if scorer.needs_proba:
            score = scorer.score_func(y, y_proba, **scorer.kwargs)

        elif scorer.needs_threshold:
            if y_proba is not None:
                score = scorer.score_func(y, y_proba, **scorer.kwargs)
            elif df is not None:
                score = scorer.score_func(y, df, **scorer.kwargs)
            else:
                score = scorer.score_func(y, y_pred, **scorer.kwargs)

        else:
            score = scorer.score_func(y, y_pred, **scorer.kwargs)

        sign = 1 if scorer.greater_is_better else -1
        scores.append(sign * score)

    return np.array(scores)


def get_scorer(scoring):
    if isinstance(scoring, six.string_types):
        try:
            scorer = SCORERS[scoring]
        except KeyError:
            raise ValueError('%r is not a valid scoring value. '
                             'Valid options are %s' % (scoring,
                             sorted(SCORERS.keys())))
    else:
        scorer = scoring
    return scorer


def _passthrough_scorer(estimator, *args, **kwargs):
    """Function that wraps estimator.score"""
    return estimator.score(*args, **kwargs)


def check_scoring(estimator, scoring=None, allow_none=False, loss_func=None,
                  score_func=None, score_overrides_loss=False):
    """Determine scorer from user options.

    A TypeError will be thrown if the estimator cannot be scored.

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    allow_none : boolean, optional, default: False
        If no scoring is specified and the estimator has no score function, we
        can either return None or raise an exception.

    Returns
    -------
    scoring : callable
        A scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
    """
    has_scoring = not (scoring is None and loss_func is None and
                       score_func is None)
    if not hasattr(estimator, 'fit'):
        raise TypeError("estimator should a be an estimator implementing "
                        "'fit' method, %r was passed" % estimator)
    elif hasattr(estimator, 'predict') and has_scoring:
        scorer = None
        if loss_func is not None or score_func is not None:
            if loss_func is not None:
                warn("Passing a loss function is "
                    "deprecated and will be removed in 0.15. "
                    "Either use strings or score objects. "
                    "The relevant new parameter is called ''scoring''. ",
                    category=DeprecationWarning, stacklevel=2)
                scorer = make_scorer(loss_func, greater_is_better=False)
            if score_func is not None:
                warn("Passing function as ``score_func`` is "
                    "deprecated and will be removed in 0.15. "
                    "Either use strings or score objects. "
                    "The relevant new parameter is called ''scoring''.",
                    category=DeprecationWarning, stacklevel=2)
                if loss_func is None or score_overrides_loss:
                    scorer = make_scorer(score_func)
        else:
            scorer = get_scorer(scoring)
        return scorer
    elif hasattr(estimator, 'score'):
        return _passthrough_scorer
    elif not has_scoring:
        if allow_none:
            return None
        raise TypeError(
            "If no scoring is specified, the estimator passed should "
            "have a 'score' method. The estimator %r does not." % estimator)
    else:
        raise TypeError(
            "The estimator passed should have a 'score' or a 'predict' "
            "method. The estimator %r does not." % estimator)


def make_scorer(score_func, greater_is_better=True, needs_proba=False,
                needs_threshold=False, **kwargs):
    """Make a scorer from a performance metric or loss function.

    This factory function wraps scoring functions for use in GridSearchCV
    and cross_val_score. It takes a score function, such as ``accuracy_score``,
    ``mean_squared_error``, ``adjusted_rand_index`` or ``average_precision``
    and returns a callable that scores an estimator's output.

    Parameters
    ----------
    score_func : callable,
        Score function (or loss function) with signature
        ``score_func(y, y_pred, **kwargs)``.

    greater_is_better : boolean, default=True
        Whether score_func is a score function (default), meaning high is good,
        or a loss function, meaning low is good. In the latter case, the
        scorer object will sign-flip the outcome of the score_func.

    needs_proba : boolean, default=False
        Whether score_func requires predict_proba to get probability estimates
        out of a classifier.

    needs_threshold : boolean, default=False
        Whether score_func takes a continuous decision certainty.
        This only works for binary classification using estimators that
        have either a decision_function or predict_proba method.

        For example ``average_precision`` or the area under the roc curve
        can not be computed using discrete predictions alone.

    **kwargs : additional arguments
        Additional parameters to be passed to score_func.

    Returns
    -------
    scorer : callable
        Callable object that returns a scalar score; greater is better.

    Examples
    --------
    >>> from sklearn.metrics import fbeta_score, make_scorer
    >>> ftwo_scorer = make_scorer(fbeta_score, beta=2)
    >>> ftwo_scorer
    make_scorer(fbeta_score, beta=2)
    >>> from sklearn.grid_search import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]},
    ...                     scoring=ftwo_scorer)
    """
    if needs_proba and needs_threshold:
        raise ValueError("Set either needs_proba or needs_threshold to True,"
                         " but not both.")

    return _Scorer(score_func,
                   greater_is_better=greater_is_better,
                   needs_proba=needs_proba,
                   needs_threshold=needs_threshold,
                   kwargs=kwargs)


# Standard regression scores
r2_scorer = make_scorer(r2_score)
mean_squared_error_scorer = make_scorer(mean_squared_error,
                                        greater_is_better=False)
mean_absolute_error_scorer = make_scorer(mean_absolute_error,
                                         greater_is_better=False)

# Standard Classification Scores
accuracy_scorer = make_scorer(accuracy_score)
f1_scorer = make_scorer(f1_score)

# Score functions that need decision values
roc_auc_scorer = make_scorer(roc_auc_score, greater_is_better=True,
                             needs_threshold=True)
average_precision_scorer = make_scorer(average_precision_score,
                                       needs_threshold=True)
precision_scorer = make_scorer(precision_score)
recall_scorer = make_scorer(recall_score)

# Score function for probabilistic classification
log_loss_scorer = make_scorer(log_loss, greater_is_better=False,
                              needs_proba=True)

# Clustering scores
adjusted_rand_scorer = make_scorer(adjusted_rand_score)

SCORERS = dict(r2=r2_scorer,
               mean_absolute_error=mean_absolute_error_scorer,
               mean_squared_error=mean_squared_error_scorer,
               accuracy=accuracy_scorer, f1=f1_scorer, roc_auc=roc_auc_scorer,
               average_precision=average_precision_scorer,
               precision=precision_scorer, recall=recall_scorer,
               log_loss=log_loss_scorer,
               adjusted_rand_score=adjusted_rand_scorer)
