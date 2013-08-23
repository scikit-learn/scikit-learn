"""Bagging meta-estimator."""

# Author: Gilles Louppe <g.louppe@gmail.com>
# License: BSD 3 clause

from __future__ import division

import itertools
import numbers
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod
from inspect import getargspec

from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..externals import six
from ..metrics import r2_score
from ..utils import check_random_state, check_arrays
from ..utils.fixes import bincount, unique

from .base import BaseEnsemble

__all__ = ["BaggingClassifier",
           "BaggingRegressor"]

MAX_INT = np.iinfo(np.int32).max


def _parallel_build_estimators(n_estimators, ensemble, X, y, sample_weight,
                               seeds, verbose):
    """Private function used to build a batch of estimators within a job."""
    # Retrieve settings
    n_samples, n_features = X.shape
    max_samples = ensemble.max_samples
    max_features = ensemble.max_features

    if (not isinstance(max_samples, (numbers.Integral, np.integer)) and
        (0.0 < max_samples <= 1.0)):
        max_samples = int(max_samples * n_samples)

    if (not isinstance(max_features, (numbers.Integral, np.integer)) and
        (0.0 < max_features <= 1.0)):
        max_features = int(max_features * n_features)

    bootstrap = ensemble.bootstrap
    bootstrap_features = ensemble.bootstrap_features
    support_sample_weight = ("sample_weight" in
                             getargspec(ensemble.base_estimator.fit)[0])

    # Build estimators
    estimators = []
    estimators_features = []

    for i in range(n_estimators):
        if verbose > 1:
            print("building estimator %d of %d" % (i + 1, n_estimators))

        random_state = check_random_state(seeds[i])
        seed = check_random_state(random_state.randint(MAX_INT))
        estimator = ensemble._make_estimator(append=False)

        try:  # Not all estimator accept a random_state
            estimator.set_params(random_state=seed)
        except ValueError:
            pass

        # Draw features
        if bootstrap_features:
            features = random_state.randint(0, n_features, max_features)
        else:
            features = random_state.permutation(n_features)[:max_features]

        # Draw samples, using sample weights, and then fit
        if support_sample_weight:
            if sample_weight is None:
                curr_sample_weight = np.ones((n_samples,))
            else:
                curr_sample_weight = sample_weight.copy()

            if bootstrap:
                samples = random_state.randint(0, n_samples, max_samples)
                sample_counts = bincount(samples, minlength=n_samples)
                curr_sample_weight *= sample_counts

            else:
                samples = random_state.permutation(n_samples)[max_samples:]
                curr_sample_weight[samples] = 0

            estimator.fit(X[:, features], y, sample_weight=curr_sample_weight)
            estimator.indices_ = curr_sample_weight > 0.

        # Draw samples, using a mask, and then fit
        else:
            if bootstrap:
                samples = random_state.randint(0, n_samples, max_samples)
            else:
                samples = random_state.permutation(n_samples)[:max_samples]

            sample_counts = bincount(samples, minlength=n_samples)

            estimator.fit((X[samples])[:, features], y[samples])
            estimator.indices_ = sample_counts > 0.

        estimators.append(estimator)
        estimators_features.append(features)

    return estimators, estimators_features


def _parallel_predict(estimators, estimators_features, X, n_classes):
    """Private function used to compute predictions within a job."""
    n_samples = X.shape[0]
    counts = np.zeros((n_samples, n_classes))

    for estimator, features in zip(estimators, estimators_features):
        predictions = estimator.predict(X[:, features])

        for i in range(n_samples):
            counts[i, predictions[i]] += 1

    return counts


def _parallel_predict_proba(estimators, estimators_features, X, n_classes):
    """Private function used to compute (proba-)predictions within a job."""
    n_samples = X.shape[0]
    proba = np.zeros((n_samples, n_classes))

    for estimator, features in zip(estimators, estimators_features):
        proba_estimator = estimator.predict_proba(X[:, features])

        if n_classes == len(estimator.classes_):
            proba += proba_estimator

        else:
            for j, c in enumerate(estimator.classes_):
                proba[:, c] += proba_estimator[:, j]

    return proba


def _parallel_predict_log_proba(estimators, estimators_features, X, n_classes):
    """Private function used to compute the log proba within a job."""
    if hasattr(estimators[0], "predict_log_proba"):
        log_proba = estimators[0].predict_log_proba(
            X[:, estimators_features[0]])

        for estimator, features in zip(estimators[1:],
                                       estimators_features[1:]):
            log_proba += estimator.predict_log_proba(X[:, features])

    else:
        log_proba = np.log(estimators[0].predict_proba(
            X[:, estimators_features[0]]))

        for estimator, features in zip(estimators[1:],
                                       estimators_features[1:]):
            log_proba += np.log(estimator.predict_proba(X[:, features]))

    return log_proba


def _parallel_decision_function(estimators, estimators_features, X, n_classes):
    """Private function used to compute (proba-)predictions within a job."""
    decisions = estimators[0].decision_function(X[:, estimators_features[0]])

    for estimator, features in zip(estimators[1:],
                                   estimators_features[1:]):
        decisions += estimator.decision_function(X[:, features])

    return decisions


def _parallel_predict_regression(estimators, estimators_features, X):
    """Private function used to compute predictions within a job."""
    return sum(estimator.predict(X[:, features])
               for estimator, features in zip(estimators,
                                              estimators_features))


def _partition_estimators(ensemble):
    """Private function used to partition estimators between jobs."""
    # Compute the number of jobs
    if ensemble.n_jobs == -1:
        n_jobs = min(cpu_count(), ensemble.n_estimators)

    else:
        n_jobs = min(ensemble.n_jobs, ensemble.n_estimators)

    # Partition estimators between jobs
    n_estimators = [ensemble.n_estimators // n_jobs] * n_jobs

    for i in range(ensemble.n_estimators % n_jobs):
        n_estimators[i] += 1

    starts = [0] * (n_jobs + 1)

    for i in range(1, n_jobs + 1):
        starts[i] = starts[i - 1] + n_estimators[i - 1]

    return n_jobs, n_estimators, starts


class BaseBagging(six.with_metaclass(ABCMeta, BaseEnsemble)):
    """Base class for Bagging meta-estimator.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    @abstractmethod
    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 max_samples=1.0,
                 max_features=1.0,
                 bootstrap=True,
                 bootstrap_features=False,
                 oob_score=False,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):
        super(BaseBagging, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators)

        self.max_samples = max_samples
        self.max_features = max_features
        self.bootstrap = bootstrap
        self.bootstrap_features = bootstrap_features
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = random_state

        #self.n_features_ = None
        #self.classes_ = None
        #self.n_classes_ = None

        self.verbose = verbose

    def fit(self, X, y, sample_weight=None):
        """Build a Bagging ensemble of estimators from the training
           set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.

        Returns
        -------
        self : object
            Returns self.
        """
        random_state = check_random_state(self.random_state)

        # Convert data
        X, y = check_arrays(X, y)

        # Remap output
        n_samples, self.n_features_ = X.shape
        y = self._validate_y(y)

        # Check parameters
        if isinstance(self.max_samples, (numbers.Integral, np.integer)):
            max_samples = self.max_samples
        else:  # float
            max_samples = int(self.max_samples * X.shape[0])

        if not (0 < max_samples <= X.shape[0]):
            raise ValueError("max_samples must be in (0, n_samples]")

        if isinstance(self.max_features, (numbers.Integral, np.integer)):
            max_features = self.max_features
        else:  # float
            max_features = int(self.max_features * self.n_features_)

        if not (0 < max_features <= self.n_features_):
            raise ValueError("max_features must be in (0, n_features]")

        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available"
                             " if bootstrap=True")

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)
        seeds = random_state.randint(MAX_INT, size=self.n_estimators)

        all_results = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_build_estimators)(
                n_estimators[i],
                self,
                X,
                y,
                sample_weight,
                seeds[starts[i]:starts[i + 1]],
                verbose=self.verbose)
            for i in range(n_jobs))

        # Reduce
        self.estimators_ = list(itertools.chain(
            *(t[0] for t in all_results)))
        self.estimators_features_ = list(itertools.chain(
            *(t[1] for t in all_results)))

        if self.oob_score:
            self._set_oob_score(X, y)

        return self

    @abstractmethod
    def _set_oob_score(self, X, y):
        """Calculate out of bag predictions and score."""

    def _validate_y(self, y):
        # Default implementation
        return y


class BaggingClassifier(BaseBagging, ClassifierMixin):
    """A Bagging classifier.

    A Bagging classifier is an ensemble meta-estimator that fits base
    classifiers each on random subsets of the original dataset and then
    aggregate their individual predictions (either by voting or by averaging)
    to form a final prediction. Such a meta-estimator can typically be used as
    a way to reduce the variance of a black-box estimator (e.g., a decision
    tree), by introducing randomization into its construction procedure and
    then making an ensemble out of it.

    This algorithm encompasses several works from the literature. When random
    subsets of the dataset are drawn as random subsets of the samples, then
    this algorithm is known as Pasting [1]. If samples are drawn with
    replacement, then the method is known as Bagging [2]. When random subsets
    of the dataset are drawn as random subsets of the features, then the method
    is known as Random Subspaces [3]. Finally, when base estimators are built
    on subsets of both samples and features, then the method is known as
    Random Patches [4].

    Parameters
    ----------
    base_estimator : object
        The base estimator to fit on random subsets of the dataset.

    n_estimators : int, optional (default=10)
        The number of base estimators in the ensemble.

    max_samples : int or float, optional (default=1.0)
        The number of samples to draw from X to train each base estimator.
            - If int, then draw `max_samples` samples.
            - If float, then draw `max_samples * X.shape[0]` samples.

    max_features : int or float, optional (default=1.0)
        The number of features to draw from X to train each base estimator.
            - If int, then draw `max_features` features.
            - If float, then draw `max_features * X.shape[1]` features.

    bootstrap : boolean, optional (default=False)
        Whether samples are drawn with replacement.

    bootstrap_features : boolean, optional (default=False)
        Whether features are drawn with replacement.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : int, optional (default=1)
        The number of jobs to run in parallel for both `fit` and `predict`.
        If -1, then the number of jobs is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controls the verbosity of the building process.

    Attributes
    ----------
    `estimators_`: list of estimators
        The collection of fitted sub-estimators.

    `classes_`: array of shape = [n_classes]
        The classes labels.

    `n_classes_`: int or list
        The number of classes.

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array of shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_decision_function_` might contain NaN.

    References
    ----------

    .. [1] L. Breiman, "Pasting small votes for classification in large
           databases and on-line", Machine Learning, 36(1), 85-103, 1999.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning, 24(2), 123-140,
           1996.

    .. [3] T. Ho, "The random subspace method for constructing decision
           forests", 1998.

    .. [4] G. Louppe and P. Geurts, "Ensembles on Random Patches", 2012.
    """
    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 max_samples=1.0,
                 max_features=1.0,
                 bootstrap=True,
                 bootstrap_features=False,
                 oob_score=False,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):

        super(BaggingClassifier, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            max_samples=max_samples,
            max_features=max_features,
            bootstrap=bootstrap,
            bootstrap_features=bootstrap_features,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def _set_oob_score(self, X, y):
        n_classes_ = self.n_classes_
        classes_ = self.classes_
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples, n_classes_))

        for estimator, features in zip(self.estimators_,
                                       self.estimators_features_):
            mask = np.ones(n_samples, dtype=np.bool)
            mask[estimator.indices_] = False

            try:
                predictions[mask, :] += estimator.predict_proba(
                    (X[mask, :])[:, features])

            except:
                p = estimator.predict((X[mask, :])[:, features])
                j = 0

                for i in range(n_samples):
                    if mask[i]:
                        predictions[i, p[j]] += 1
                        j += 1

        if (predictions.sum(axis=1) == 0).any():
            warn("Some inputs do not have OOB scores. "
                 "This probably means too few estimators were used "
                 "to compute any reliable oob estimates.")

        oob_decision_function = (predictions /
                                 predictions.sum(axis=1)[:, np.newaxis])
        oob_score = np.mean((y == classes_.take(np.argmax(predictions, axis=1),
                                                axis=0)))

        self.oob_decision_function_ = oob_decision_function
        self.oob_score_ = oob_score

    def _validate_y(self, y):
        self.classes_, y = unique(y, return_inverse=True)
        self.n_classes_ = len(self.classes_)

        return y

    def predict(self, X):
        """Predict class for X.

        The predicted class of an input sample is computed as the majority
        prediction of the estimators in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted classes.
        """
        try:
            return self.classes_.take(np.argmax(self.predict_proba(X), axis=1),
                                      axis=0)

        except:
            # Check data
            X, = check_arrays(X)

            # Parallel loop
            n_jobs, n_estimators, starts = _partition_estimators(self)

            all_counts = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
                delayed(_parallel_predict)(
                    self.estimators_[starts[i]:starts[i + 1]],
                    self.estimators_features_[starts[i]:starts[i + 1]],
                    X,
                    self.n_classes_)
                for i in range(n_jobs))

            # Reduce
            counts = all_counts[0]

            for j in range(1, len(all_counts)):
                counts += all_counts[j]

            return self.classes_.take(np.argmax(counts, axis=1), axis=0)

    def predict_proba(self, X):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the mean predicted class probabilities of the estimators in the
        ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes]
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        # Check data
        X, = check_arrays(X)

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)

        all_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_proba)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X,
                self.n_classes_)
            for i in range(n_jobs))

        # Reduce
        proba = all_proba[0]

        for j in range(1, len(all_proba)):
            proba += all_proba[j]

        proba /= self.n_estimators

        return proba

    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the mean predicted class log-probabilities of the estimators in the
        ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        p : array of shape = [n_samples, n_classes]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        # Check data
        X, = check_arrays(X)

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)

        all_log_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_log_proba)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X,
                self.n_classes_)
            for i in range(n_jobs))

        # Reduce
        log_proba = all_log_proba[0]

        for j in range(1, len(all_log_proba)):
            log_proba += all_log_proba[j]

        log_proba /= self.n_estimators

        return log_proba

    def decision_function(self, X):
        """Average of the decision functions of the base classifiers.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : array, shape = [n_samples, k]
            The decision function of the input samples. Classes are
            ordered by arithmetical order. Regression and binary
            classification are special cases with ``k == 1``,
            otherwise ``k==n_classes``.
        """
        # Check data
        X, = check_arrays(X)

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)

        all_decisions = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_decision_function)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X,
                self.n_classes_)
            for i in range(n_jobs))

        # Reduce
        decisions = all_decisions[0]

        for j in range(1, len(all_decisions)):
            decisions += all_decisions[j]

        decisions /= self.n_estimators

        return decisions


class BaggingRegressor(BaseBagging, RegressorMixin):
    """A Bagging regressor.

    A Bagging regressor is an ensemble meta-estimator that fits base
    regressors each on random subsets of the original dataset and then
    aggregate their individual predictions (either by voting or by averaging)
    to form a final prediction. Such a meta-estimator can typically be used as
    a way to reduce the variance of a black-box estimator (e.g., a decision
    tree), by introducing randomization into its construction procedure and
    then making an ensemble out of it.

    This algorithm encompasses several works from the literature. When random
    subsets of the dataset are drawn as random subsets of the samples, then
    this algorithm is known as Pasting [1]. If samples are drawn with
    replacement, then the method is known as Bagging [2]. When random subsets
    of the dataset are drawn as random subsets of the features, then the method
    is known as Random Subspaces [3]. Finally, when base estimators are built
    on subsets of both samples and features, then the method is known as
    Random Patches [4].

    Parameters
    ----------
    base_estimator : object
        The base estimator to fit on random subsets of the dataset.

    n_estimators : int, optional (default=10)
        The number of base estimators in the ensemble.

    max_samples : int or float, optional (default=1.0)
        The number of samples to draw from X to train each base estimator.
            - If int, then draw `max_samples` samples.
            - If float, then draw `max_samples * X.shape[0]` samples.

    max_features : int or float, optional (default=1.0)
        The number of features to draw from X to train each base estimator.
            - If int, then draw `max_features` features.
            - If float, then draw `max_features * X.shape[1]` features.

    bootstrap : boolean, optional (default=False)
        Whether samples are drawn with replacement.

    bootstrap_features : boolean, optional (default=False)
        Whether features are drawn with replacement.

    oob_score : bool
        Whether to use out-of-bag samples to estimate
        the generalization error.

    n_jobs : int, optional (default=1)
        The number of jobs to run in parallel for both `fit` and `predict`.
        If -1, then the number of jobs is set to the number of cores.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    verbose : int, optional (default=0)
        Controls the verbosity of the building process.

    Attributes
    ----------
    `estimators_`: list of estimators
        The collection of fitted sub-estimators.

    `oob_score_` : float
        Score of the training dataset obtained using an out-of-bag estimate.

    `oob_decision_function_` : array of shape = [n_samples, n_classes]
        Decision function computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_decision_function_` might contain NaN.

    References
    ----------

    .. [1] L. Breiman, "Pasting small votes for classification in large
           databases and on-line", Machine Learning, 36(1), 85-103, 1999.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning, 24(2), 123-140,
           1996.

    .. [3] T. Ho, "The random subspace method for constructing decision
           forests", 1998.

    .. [4] G. Louppe and P. Geurts, "Ensembles on Random Patches", 2012.
    """

    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 max_samples=1.0,
                 max_features=1.0,
                 bootstrap=True,
                 bootstrap_features=False,
                 oob_score=False,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):
        super(BaggingRegressor, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            max_samples=max_samples,
            max_features=max_features,
            bootstrap=bootstrap,
            bootstrap_features=bootstrap_features,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def predict(self, X):
        """Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the estimators in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y: array of shape = [n_samples]
            The predicted values.
        """
        # Check data
        X, = check_arrays(X)

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)

        all_y_hat = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]],
                self.estimators_features_[starts[i]:starts[i + 1]],
                X)
            for i in range(n_jobs))

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat

    def _set_oob_score(self, X, y):
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples,))
        n_predictions = np.zeros((n_samples,))

        for estimator, features in zip(self.estimators_,
                                       self.estimators_features_):
            mask = np.ones(n_samples, dtype=np.bool)
            mask[estimator.indices_] = False

            predictions[mask] += estimator.predict((X[mask, :])[:, features])
            n_predictions[mask] += 1

        if (n_predictions == 0).any():
            warn("Some inputs do not have OOB scores. "
                 "This probably means too few estimators were used "
                 "to compute any reliable oob estimates.")
            n_predictions[n_predictions == 0] = 1

        predictions /= n_predictions

        self.oob_prediction_ = predictions
        self.oob_score_ = r2_score(y, predictions)
