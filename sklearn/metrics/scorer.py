"""
The :mod:`sklearn.metrics.scorer` submodule implements a flexible
interface for model selection and evaluation using
arbitrary score functions.

A Scorer object is a callable that can be passed to
:class:`sklearn.grid_search.GridSearchCV` or
:func:`sklearn.cross_validation.cross_val_score` as the ``scoring`` parameter,
to specify how a model should be evaluated.

The signature of the call is ``(estimator, X, y)`` where ``estimator``
is the model to be evaluated, ``X`` is the test data and ``y`` is the
ground truth labeling (or ``None`` in the case of unsupervised models).
"""

# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
# Liscence: Simplified BSD

from abc import ABCMeta, abstractmethod

import numpy as np

from . import (r2_score, mean_squared_error, accuracy_score,
               auc_score, average_precision_score, precision_score,
               recall_score, precision_recall_fscore_support)

from .cluster import adjusted_rand_score

class BaseScorer(object):
    __metaclass__ = ABCMeta

    def __init__(self, greater_is_better=True):
        self.greater_is_better = greater_is_better

    def calc_scores(self, estimator, X, y=None):
        """Calculate one or more scores for X against y using the provided
        estimator. While __call__ calculates a single score, this may return
        multiple.

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring.
            If ``needs_threshold`` is True, estimator needs
            to provide ``decision_function`` or ``predict_proba``.
            Otherwise, estimator needs to provide ``predict``.

        X : array-like or sparse matrix
            Test data that will be scored by the estimator.

        y : array-like
            True prediction for X.

        Returns
        -------
        scores : iterable of (name, score) pairs
            Scores of the estimator's predictions of X with respect to y.
            Names must be distinct, and exactly one name must be 'score', whose
            score corresponds to the result of `__call__`.
        """
        yield ('score', self(estimator, X, y))

    @abstractmethod
    def __call__(self, estimator, X, y=None):
        """Score X and y using the provided estimator.

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring.
            If ``needs_threshold`` is True, estimator needs
            to provide ``decision_function`` or ``predict_proba``.
            Otherwise, estimator needs to provide ``predict``.

        X : array-like or sparse matrix
            Test data that will be scored by the estimator.

        y : array-like
            True prediction for X.

        Returns
        -------
        score : float
            The score of estimator's prediction of X.
        """
        pass


class Scorer(BaseScorer):
    """Flexible scores for any estimator.

    This class wraps estimator scoring functions for the use in GridSearchCV
    and cross_val_score. It takes a score function, such as ``accuracy_score``,
    ``mean_squared_error``, ``adjusted_rand_index`` or ``average_precision``
    and provides a call method.

    Parameters
    ----------
    score_func : callable,
        Score function (or loss function) with signature
        ``score_func(y, y_pred, **kwargs)``.

    greater_is_better : boolean, default=True
        Whether score_func is a score function (default), meaning high is good,
        or a loss function, meaning low is good.

    needs_threshold : bool, default=False
        Whether score_func takes a continuous decision certainty.
        For example ``average_precision`` or the area under the roc curve
        can not be computed using predictions alone, but need the output of
        ``decision_function`` or ``predict_proba``.

    **kwargs : additional arguments
        Additional parameters to be passed to score_func.

    Examples
    --------
    >>> from sklearn.metrics import fbeta_score, Scorer
    >>> ftwo_scorer = Scorer(fbeta_score, beta=2)
    >>> from sklearn.grid_search import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]},
    ...                     scoring=ftwo_scorer)
    """
    def __init__(self, score_func, greater_is_better=True,
                 needs_threshold=False, **kwargs):
        super(Scorer, self).__init__(greater_is_better)
        self.score_func = score_func
        self.needs_threshold = needs_threshold
        self.kwargs = kwargs

    def __repr__(self):
        kwargs_string = "".join([", %s=%s" % (str(k), str(v))
                                 for k, v in self.kwargs.items()])
        return ("Scorer(score_func=%s, greater_is_better=%s, needs_thresholds="
                "%s%s)" % (self.score_func.__name__, self.greater_is_better,
                           self.needs_threshold, kwargs_string))

    def __call__(self, estimator, X, y=None):
        """Score X and y using the provided estimator.

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring.
            If ``needs_threshold`` is True, estimator needs
            to provide ``decision_function`` or ``predict_proba``.
            Otherwise, estimator needs to provide ``predict``.

        X : array-like or sparse matrix
            Test data that will be scored by the estimator.

        y : array-like
            True prediction for X.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        if self.needs_threshold:
            if len(np.unique(y)) > 2:
                raise ValueError("This classification score only "
                                 "supports binary classification.")
            try:
                y_pred = estimator.decision_function(X).ravel()
            except (NotImplementedError, AttributeError):
                y_pred = estimator.predict_proba(X)[:, 1]
        else:
            y_pred = estimator.predict(X)
        if y is not None:
            return self.score_func(y, y_pred, **self.kwargs)
        else:
            return self.score_func(y_pred, **self.kwargs)


class PRFScorer(Scorer):
    """Scorer to optimise F score while also providing precision and recall.

    Parameters
    ----------
    **kwargs : additional arguments
        Additional parameters to be passed to
        `metrics.precision_recall_fscore_support`.
    """

    def __init__(self, **kwargs):
        if 'average' not in kwargs:
            kwargs['average'] = 'weighted'
        super(PRFScorer, self).__init__(precision_recall_fscore_support, **kwargs)

    def __repr__(self):
        kwargs_string = "".join([", %s=%s" % (str(k), str(v))
                                 for k, v in self.kwargs.items()])
        return 'PRFScorer(%s)' % kwargs_string

    def calc_scores(self, estimator, X, y):
        """
        Calculates F score, precision and recall

        Parameters
        ----------
        estimator : object
            Trained estimator to use for scoring.
            If ``needs_threshold`` is True, estimator needs
            to provide ``decision_function`` or ``predict_proba``.
            Otherwise, estimator needs to provide ``predict``.

        X : array-like or sparse matrix
            Test data that will be scored by the estimator.

        y : array-like
            True prediction for X.

        Returns
        -------
        scores : list of (name, score) pairs
            providing names 'score', 'precision' and 'recall'
        """
        p, r, f, support = super(PRFScorer, self).__call__(estimator, X, y)
        return [
                ('score', f),
                ('precision', p),
                ('recall', r),
        ]

    def __call__(self, estimator, X, y):
        p, r, f, support = super(PRFScorer, self).__call__(estimator, X, y)
        return f


class WrapScorer(BaseScorer):
    """Scores by passing the estimator and data to a given function
    
    Parameters
    ----------
    score_fn : function with signature of `Scorer.__call__`
        A function which returns a score given an estimator, instances and
        ground truth if available.

    greater_is_better : boolean, default=True
        Whether score_func is a score function (default), meaning high is good,
        or a loss function, meaning low is good.
    """

    def __init__(self, score_fn, greater_is_better=True):
        super(WrapScorer, self).__init__(greater_is_better)
        self.score_fn = score_fn

    def __call__(self, estimator, X, y=None):
        if y is None:
            return self.score_fn(estimator, X)
        return self.score_fn(estimator, X, y)


class EstimatorScorer(BaseScorer):
    """Scores by calling the estimator's score method."""

    def __call__(self, estimator, X, y=None):
        if y is None:
            return estimator.score(X)
        return estimator.score(X, y)


# Standard regression scores
r2_scorer = Scorer(r2_score)
mse_scorer = Scorer(mean_squared_error, greater_is_better=False)

# Standard Classification Scores
accuracy_scorer = Scorer(accuracy_score)
f1_scorer = PRFScorer()

# Score functions that need decision values
auc_scorer = Scorer(auc_score, greater_is_better=True, needs_threshold=True)
average_precision_scorer = Scorer(average_precision_score,
                                  needs_threshold=True)
precision_scorer = Scorer(precision_score)
recall_scorer = Scorer(recall_score)

# Clustering scores
ari_scorer = Scorer(adjusted_rand_score)

SCORERS = dict(r2=r2_scorer, mse=mse_scorer, accuracy=accuracy_scorer,
               f1=PRFScorer(), roc_auc=auc_scorer,
               average_precision=average_precision_scorer,
               precision=precision_scorer, recall=recall_scorer,
               ari=ari_scorer)
