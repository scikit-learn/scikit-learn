"""
The :mod:`sklearn.metrics.scorer` submodule implements a flexible
interface for model selection and evaluation using
arbitrary score functions.

A scorer object is a callable that can be passed to
:class:`sklearn.model_selection.GridSearchCV` or
:func:`sklearn.model_selection.cross_val_score` as the ``scoring``
parameter, to specify how a model should be evaluated.

The signature of the call is ``(estimator, X, y)`` where ``estimator``
is the model to be evaluated, ``X`` is the test data and ``y`` is the
ground truth labeling (or ``None`` in the case of unsupervised models).
"""

# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Lars Buitinck
#          Arnaud Joly <arnaud.v.joly@gmail.com>
# License: Simplified BSD

from abc import ABCMeta, abstractmethod
import warnings

import numpy as np

from . import (r2_score, median_absolute_error, mean_absolute_error,
               mean_squared_error, mean_squared_log_error, accuracy_score,
               f1_score, roc_auc_score, average_precision_score,
               precision_score, recall_score, log_loss)

from .cluster import adjusted_rand_score
from .cluster import homogeneity_score
from .cluster import completeness_score
from .cluster import v_measure_score
from .cluster import mutual_info_score
from .cluster import adjusted_mutual_info_score
from .cluster import normalized_mutual_info_score
from .cluster import fowlkes_mallows_score

from ..utils.multiclass import type_of_target
from ..externals import six
from ..base import is_regressor


class _BaseScorer(six.with_metaclass(ABCMeta, object)):
    def __init__(self, score_func, sign, kwargs):
        self._kwargs = kwargs
        self._score_func = score_func
        self._sign = sign
        # XXX After removing the deprecated scorers (v0.20) remove the
        # XXX deprecation_msg property again and remove __call__'s body again
        self._deprecation_msg = None

    @abstractmethod
    def __call__(self, estimator, X, y, sample_weight=None):
        if self._deprecation_msg is not None:
            warnings.warn(self._deprecation_msg,
                          category=DeprecationWarning,
                          stacklevel=2)

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
    def __call__(self, estimator, X, y_true, sample_weight=None):
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

        sample_weight : array-like, optional (default=None)
            Sample weights.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        super(_PredictScorer, self).__call__(estimator, X, y_true,
                                             sample_weight=sample_weight)
        y_pred = estimator.predict(X)
        if sample_weight is not None:
            return self._sign * self._score_func(y_true, y_pred,
                                                 sample_weight=sample_weight,
                                                 **self._kwargs)
        else:
            return self._sign * self._score_func(y_true, y_pred,
                                                 **self._kwargs)


class _ProbaScorer(_BaseScorer):
    def __call__(self, clf, X, y, sample_weight=None):
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

        sample_weight : array-like, optional (default=None)
            Sample weights.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        super(_ProbaScorer, self).__call__(clf, X, y,
                                           sample_weight=sample_weight)
        y_pred = clf.predict_proba(X)
        if sample_weight is not None:
            return self._sign * self._score_func(y, y_pred,
                                                 sample_weight=sample_weight,
                                                 **self._kwargs)
        else:
            return self._sign * self._score_func(y, y_pred, **self._kwargs)

    def _factory_args(self):
        return ", needs_proba=True"


class _ThresholdScorer(_BaseScorer):
    def __call__(self, clf, X, y, sample_weight=None):
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

        sample_weight : array-like, optional (default=None)
            Sample weights.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        super(_ThresholdScorer, self).__call__(clf, X, y,
                                               sample_weight=sample_weight)
        y_type = type_of_target(y)
        if y_type not in ("binary", "multilabel-indicator"):
            raise ValueError("{0} format is not supported".format(y_type))

        if is_regressor(clf):
            y_pred = clf.predict(X)
        else:
            try:
                y_pred = clf.decision_function(X)

                # For multi-output multi-class estimator
                if isinstance(y_pred, list):
                    y_pred = np.vstack(p for p in y_pred).T

            except (NotImplementedError, AttributeError):
                y_pred = clf.predict_proba(X)

                if y_type == "binary":
                    y_pred = y_pred[:, 1]
                elif isinstance(y_pred, list):
                    y_pred = np.vstack([p[:, -1] for p in y_pred]).T

        if sample_weight is not None:
            return self._sign * self._score_func(y, y_pred,
                                                 sample_weight=sample_weight,
                                                 **self._kwargs)
        else:
            return self._sign * self._score_func(y, y_pred, **self._kwargs)

    def _factory_args(self):
        return ", needs_threshold=True"


def get_scorer(scoring):
    if isinstance(scoring, six.string_types):
        try:
            scorer = SCORERS[scoring]
        except KeyError:
            scorers = [scorer for scorer in SCORERS
                       if SCORERS[scorer]._deprecation_msg is None]
            raise ValueError('%r is not a valid scoring value. '
                             'Valid options are %s'
                             % (scoring, sorted(scorers)))
    else:
        scorer = scoring
    return scorer


def _passthrough_scorer(estimator, *args, **kwargs):
    """Function that wraps estimator.score"""
    return estimator.score(*args, **kwargs)


def check_scoring(estimator, scoring=None, allow_none=False):
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
    has_scoring = scoring is not None
    if not hasattr(estimator, 'fit'):
        raise TypeError("estimator should be an estimator implementing "
                        "'fit' method, %r was passed" % estimator)
    if isinstance(scoring, six.string_types):
        return get_scorer(scoring)
    elif has_scoring:
        # Heuristic to ensure user has not passed a metric
        module = getattr(scoring, '__module__', None)
        if hasattr(module, 'startswith') and \
           module.startswith('sklearn.metrics.') and \
           not module.startswith('sklearn.metrics.scorer') and \
           not module.startswith('sklearn.metrics.tests.'):
            raise ValueError('scoring value %r looks like it is a metric '
                             'function rather than a scorer. A scorer should '
                             'require an estimator as its first parameter. '
                             'Please use `make_scorer` to convert a metric '
                             'to a scorer.' % scoring)
        return get_scorer(scoring)
    elif hasattr(estimator, 'score'):
        return _passthrough_scorer
    elif allow_none:
        return None
    else:
        raise TypeError(
            "If no scoring is specified, the estimator passed should "
            "have a 'score' method. The estimator %r does not." % estimator)


def make_scorer(score_func, greater_is_better=True, needs_proba=False,
                needs_threshold=False, **kwargs):
    """Make a scorer from a performance metric or loss function.

    This factory function wraps scoring functions for use in GridSearchCV
    and cross_val_score. It takes a score function, such as ``accuracy_score``,
    ``mean_squared_error``, ``adjusted_rand_index`` or ``average_precision``
    and returns a callable that scores an estimator's output.

    Read more in the :ref:`User Guide <scoring>`.

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
    >>> from sklearn.model_selection import GridSearchCV
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
neg_mean_squared_error_scorer = make_scorer(mean_squared_error,
                                            greater_is_better=False)
deprecation_msg = ('Scoring method mean_squared_error was renamed to '
                   'neg_mean_squared_error in version 0.18 and will '
                   'be removed in 0.20.')
mean_squared_error_scorer = make_scorer(mean_squared_error,
                                        greater_is_better=False)
mean_squared_error_scorer._deprecation_msg = deprecation_msg
neg_mean_squared_log_error_scorer = make_scorer(mean_squared_log_error,
                                                greater_is_better=False)
neg_mean_absolute_error_scorer = make_scorer(mean_absolute_error,
                                             greater_is_better=False)
deprecation_msg = ('Scoring method mean_absolute_error was renamed to '
                   'neg_mean_absolute_error in version 0.18 and will '
                   'be removed in 0.20.')
mean_absolute_error_scorer = make_scorer(mean_absolute_error,
                                         greater_is_better=False)
mean_absolute_error_scorer._deprecation_msg = deprecation_msg
neg_median_absolute_error_scorer = make_scorer(median_absolute_error,
                                               greater_is_better=False)
deprecation_msg = ('Scoring method median_absolute_error was renamed to '
                   'neg_median_absolute_error in version 0.18 and will '
                   'be removed in 0.20.')
median_absolute_error_scorer = make_scorer(median_absolute_error,
                                           greater_is_better=False)
median_absolute_error_scorer._deprecation_msg = deprecation_msg


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
neg_log_loss_scorer = make_scorer(log_loss, greater_is_better=False,
                                  needs_proba=True)
deprecation_msg = ('Scoring method log_loss was renamed to '
                   'neg_log_loss in version 0.18 and will be removed in 0.20.')
log_loss_scorer = make_scorer(log_loss, greater_is_better=False,
                              needs_proba=True)
log_loss_scorer._deprecation_msg = deprecation_msg


# Clustering scores
adjusted_rand_scorer = make_scorer(adjusted_rand_score)
homogeneity_scorer = make_scorer(homogeneity_score)
completeness_scorer = make_scorer(completeness_score)
v_measure_scorer = make_scorer(v_measure_score)
mutual_info_scorer = make_scorer(mutual_info_score)
adjusted_mutual_info_scorer = make_scorer(adjusted_mutual_info_score)
normalized_mutual_info_scorer = make_scorer(normalized_mutual_info_score)
fowlkes_mallows_scorer = make_scorer(fowlkes_mallows_score)


SCORERS = dict(r2=r2_scorer,
               neg_median_absolute_error=neg_median_absolute_error_scorer,
               neg_mean_absolute_error=neg_mean_absolute_error_scorer,
               neg_mean_squared_error=neg_mean_squared_error_scorer,
               neg_mean_squared_log_error=neg_mean_squared_log_error_scorer,
               median_absolute_error=median_absolute_error_scorer,
               mean_absolute_error=mean_absolute_error_scorer,
               mean_squared_error=mean_squared_error_scorer,
               accuracy=accuracy_scorer, roc_auc=roc_auc_scorer,
               average_precision=average_precision_scorer,
               log_loss=log_loss_scorer,
               neg_log_loss=neg_log_loss_scorer,
               # Cluster metrics that use supervised evaluation
               adjusted_rand_score=adjusted_rand_scorer,
               homogeneity_score=homogeneity_scorer,
               completeness_score=completeness_scorer,
               v_measure_score=v_measure_scorer,
               mutual_info_score=mutual_info_scorer,
               adjusted_mutual_info_score=adjusted_mutual_info_scorer,
               normalized_mutual_info_score=normalized_mutual_info_scorer,
               fowlkes_mallows_score=fowlkes_mallows_scorer)


for name, metric in [('precision', precision_score),
                     ('recall', recall_score), ('f1', f1_score)]:
    SCORERS[name] = make_scorer(metric)
    for average in ['macro', 'micro', 'samples', 'weighted']:
        qualified_name = '{0}_{1}'.format(name, average)
        SCORERS[qualified_name] = make_scorer(metric, pos_label=None,
                                              average=average)
