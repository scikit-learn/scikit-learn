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
# License: Simplified BSD

from abc import ABCMeta, abstractmethod
from warnings import warn

import numpy as np

from . import (r2_score, mean_squared_error, accuracy_score, f1_score,
               roc_auc_score, average_precision_score, precision_score,
               recall_score, log_loss)

from .cluster import adjusted_rand_score
from ..externals import six


class _BaseScorer(six.with_metaclass(ABCMeta, object)):
    def __init__(self, score_func, sign, kwargs):
        self._kwargs = kwargs
        self._score_func = score_func
        self._sign = sign

    @abstractmethod
    def __call__(self, estimator, X, y):
        pass

    def __repr__(self):
        kwargs_string = "".join([", %s=%s" % (str(k), str(v))
                                 for k, v in self._kwargs.items()])
        return ("make_scorer(%s%s%s%s)"
                % (self._score_func.__name__,
                   "" if self._sign > 0 else ", greater_is_better=False",
                   self._factory_args(), kwargs_string))

    def _factory_args(self):
        """Return non-default make_scorer arguments for repr."""
        return ""


class _PredictScorer(_BaseScorer):
    def __call__(self, estimator, X, y_true):
        """Evaluate predicted target values for X relative to y_true.

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring. Must have a predict_proba
            method; the output of that is used to compute the score.

        X : array-like or sparse matrix
            Test data that will be fed to estimator.predict.

        y_true : array-like
            Gold standard target values for X.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        y_pred = estimator.predict(X)
        return self._sign * self._score_func(y_true, y_pred, **self._kwargs)


class _ProbaScorer(_BaseScorer):
    def __call__(self, clf, X, y):
        """Evaluate predicted probabilities for X relative to y_true.

        Parameters
        ----------
        clf : object
            Trained classifier to use for scoring. Must have a predict_proba
            method; the output of that is used to compute the score.

        X : array-like or sparse matrix
            Test data that will be fed to clf.predict_proba.

        y : array-like
            Gold standard target values for X. These must be class labels,
            not probabilities.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        y_pred = clf.predict_proba(X)
        return self._sign * self._score_func(y, y_pred, **self._kwargs)

    def _factory_args(self):
        return ", needs_proba=True"


class _ThresholdScorer(_BaseScorer):
    def __call__(self, clf, X, y):
        """Evaluate decision function output for X relative to y_true.

        Parameters
        ----------
        clf : object
            Trained classifier to use for scoring. Must have either a
            decision_function method or a predict_proba method; the output of
            that is used to compute the score.

        X : array-like or sparse matrix
            Test data that will be fed to clf.decision_function or
            clf.predict_proba.

        y : array-like
            Gold standard target values for X. These must be class labels,
            not decision function values.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        if len(np.unique(y)) > 2:
            raise ValueError("This classification score only "
                             "supports binary classification.")
        try:
            y_pred = clf.decision_function(X).ravel()
        except (NotImplementedError, AttributeError):
            y_pred = clf.predict_proba(X)[:, 1]
        return self._sign * self._score_func(y, y_pred, **self._kwargs)

    def _factory_args(self):
        return ", needs_threshold=True"


def _deprecate_loss_and_score_funcs(
        loss_func=None, score_func=None, scoring=None,
        score_overrides_loss=False):

    scorer = None
    if loss_func is not None or score_func is not None:

        if loss_func is not None:
            warn("Passing a loss function is "
                 "deprecated and will be removed in 0.15. "
                 "Either use strings or score objects."
                 "The relevant new parameter is called ''scoring''. ",
                 category=DeprecationWarning, stacklevel=2)
            scorer = make_scorer(loss_func, greater_is_better=False)
        if score_func is not None:
            warn("Passing function as ``score_func`` is "
                 "deprecated and will be removed in 0.15. "
                 "Either use strings or score objects."
                 "The relevant new parameter is called ''scoring''.",
                 category=DeprecationWarning, stacklevel=2)
            if loss_func is None or score_overrides_loss:
                scorer = make_scorer(score_func)

    elif isinstance(scoring, six.string_types):
        try:
            scorer = SCORERS[scoring]
        except KeyError:
            raise ValueError('%r is not a valid scoring value. '
                             'Valid options are %s' % (scoring,
                             sorted(SCORERS.keys())))
    else:
        scorer = scoring
    return scorer


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
    sign = 1 if greater_is_better else -1
    if needs_proba and needs_threshold:
        raise ValueError("Set either needs_proba or needs_threshold to True,"
                         " but not both.")
    if needs_proba:
        cls = _ProbaScorer
    elif needs_threshold:
        cls = _ThresholdScorer
    else:
        cls = _PredictScorer
    return cls(score_func, sign, kwargs)


# Standard regression scores
r2_scorer = make_scorer(r2_score)
mean_squared_error_scorer = make_scorer(mean_squared_error,
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
               mean_squared_error=mean_squared_error_scorer,
               accuracy=accuracy_scorer, f1=f1_scorer, roc_auc=roc_auc_scorer,
               average_precision=average_precision_scorer,
               precision=precision_scorer, recall=recall_scorer,
               log_loss=log_loss_scorer,
               adjusted_rand_score=adjusted_rand_scorer)
