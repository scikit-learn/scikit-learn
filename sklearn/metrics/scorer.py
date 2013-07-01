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
from collections import namedtuple

import numpy as np

from . import (r2_score, mean_squared_error, accuracy_score,
               auc_score, average_precision_score, precision_score,
               recall_score, precision_recall_fscore_support)

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
    def _score(self, clf, X, y):
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
mse_scorer = make_scorer(mean_squared_error, greater_is_better=False)

# Standard classification scores
accuracy_scorer = make_scorer(accuracy_score)


class _FPR(namedtuple("fpr", ["f_score", "precision", "recall"])):
    __slots__ = ()

    def __str__(self):
        return ("F = {0:.4f}, precision = {1:.4f}, recall = {2:.4f}"
                .format(self.f_score, self.precision, self.recall))


def f_scorer(clf, X, y_true, beta=1.):
    """Evaluate a classifier's predictions for X according to F1/F-beta score.

    Parameters
    ----------
    clf : object
        Trained classifier to evaluate.

    X : array-like or sparse matrix
        Test data that will be fed to clf.predict.

    y_true : array-like
        Gold standard target values for X.

    beta : float, optional
        The strength of recall versus precision in the F-score.

    Returns
    -------
    (fscore, precision, recall) : tuple of floats
        F-score and the scores it is based on for estimator's predictions on X
        relative to y_true. When this scorer is used inside GridSearchCV or
        similar, the first value is used for optimization.

    See also
    --------
    sklearn.metrics.precision_recall_fscore_support
    """
    # TODO support the various weightings for precision_recall_fscore_support
    # offers
    p, r, f, _ = precision_recall_fscore_support(y_true, clf.predict(X),
                                                 beta, average="weighted")
    return _FPR(f, p, r)


# Score functions that need decision values
auc_scorer = make_scorer(auc_score, greater_is_better=True,
                         needs_threshold=True)
average_precision_scorer = make_scorer(average_precision_score,
                                       needs_threshold=True)
precision_scorer = make_scorer(precision_score)
recall_scorer = make_scorer(recall_score)

# Clustering scores
ari_scorer = make_scorer(adjusted_rand_score)

SCORERS = dict(r2=r2_scorer, mse=mse_scorer, accuracy=accuracy_scorer,
               f1=f_scorer, roc_auc=auc_scorer,
               average_precision=average_precision_scorer,
               precision=precision_scorer, recall=recall_scorer,
               ari=ari_scorer)
