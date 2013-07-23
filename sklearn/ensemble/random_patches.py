"""RandomPatches meta-estimator."""

# Author: Gilles Louppe <g.louppe@gmail.com>
# License: BSD 3 clause

from __future__ import division

import inspect
import itertools
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod

from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..externals import six
from ..externals.six.moves import xrange
from ..metrics import r2_score
from ..utils import array2d, check_random_state, check_arrays, safe_asarray
from ..utils.fixes import bincount, unique
from ..utils.validation import check_arrays

from .base import BaseEnsemble

__all__ = ["RandomPatchesClassifier",
           "RandomPatchesRegressor"]

MAX_INT = np.iinfo(np.int32).max


def _parallel_build_estimators(n_estimators, ensemble, X, y,
                               support_sample_weight, sample_weight, seeds,
                               verbose):
    """Private function used to build a batch of estimators within a job."""
    n_samples = X.shape[0]
    estimators = []

    for i in range(n_estimators):
        if verbose > 1:
            print("building estimator %d of %d" % (i + 1, n_estimators))

        random_state = check_random_state(seeds[i])
        seed = random_state.randint(MAX_INT)

        estimator = ensemble._make_estimator(append=False)

        try: # Not all estimator accept a random_state
            estimator.set_params(random_state=check_random_state(seed))
        except ValueError:
            pass

        if support_sample_weight:
            if ensemble.bootstrap:
                indices = random_state.randint(0, n_samples, n_samples)
                sample_counts = bincount(indices, minlength=n_samples)

                if sample_weight is None:
                    curr_sample_weight = np.ones((n_samples,))
                else:
                    curr_sample_weight = sample_weight.copy()

                curr_sample_weight *= sample_counts

                estimator.fit(X, y, sample_weight=curr_sample_weight)
                estimator.indices_ = sample_counts > 0.

            else:
                estimator.fit(X, y, sample_weight=sample_weight)

        else:
            if ensemble.bootstrap:
                indices = random_state.randint(0, n_samples, n_samples)
                sample_counts = bincount(indices, minlength=n_samples)

                estimator.fit(X[indices], y[indices])
                estimator.indices_ = sample_counts > 0.

            else:
                estimator.fit(X, y)

        estimators.append(estimator)

    return estimators


def _parallel_predict(estimators, X, n_classes):
    """Private function used to compute a batch of predictions within a job."""
    n_samples = X.shape[0]
    counts = np.zeros((n_samples, n_classes))

    for estimator in estimators:
        predictions = estimator.predict(X)

        for i in xrange(n_samples):
            counts[i, predictions[i]] += 1

    return counts


def _parallel_predict_proba(estimators, X, n_classes):
    """Private function used to compute a batch of (proba-)predictions within
       a job."""
    n_samples = X.shape[0]
    proba = np.zeros((n_samples, n_classes))

    for estimator in estimators:
        proba_estimator = estimator.predict_proba(X)

        if n_classes == len(estimator.classes_):
            proba += proba_estimator

        else:
            for j, c in enumerate(estimator.classes_):
                proba[:, c] += proba_estimator[:, j]

    return proba


def _parallel_predict_regression(estimators, X):
    """Private function used to compute a batch of predictions within a job."""
    return sum(estimator.predict(X) for estimator in estimators)


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


class BaseRandomPatches(six.with_metaclass(ABCMeta, BaseEnsemble)):
    """Base class for Random Patches meta-estimator.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    @abstractmethod
    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 bootstrap=True,
                 oob_score=True,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):
        super(BaseRandomPatches, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators)

        self.bootstrap = bootstrap
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = random_state

        #self.n_features_ = None
        #self.classes_ = None
        #self.n_classes_ = None

        self.verbose = verbose

    def fit(self, X, y, sample_weight=None):
        """Build a Random Patches ensemble of estimators from the training
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
        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available"
                             " if bootstrap=True")

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, _ = _partition_estimators(self)

        # Precalculate the random states
        seeds = [random_state.randint(MAX_INT, size=i) for i in n_estimators]

        # Parallel loop
        support_sample_weight = "sample_weight" in inspect.getargspec(self.base_estimator.fit)[0]

        all_estimators = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_build_estimators)(
                n_estimators[i],
                self,
                X,
                y,
                support_sample_weight,
                sample_weight,
                seeds[i],
                verbose=self.verbose)
            for i in range(n_jobs))

        # Reduce
        self.estimators_ = list(itertools.chain(*all_estimators))

        if self.oob_score:
            self._set_oob_score(X, y)

        return self

    @abstractmethod
    def _set_oob_score(self, X, y):
        """Calculate out of bag predictions and score."""

    def _validate_y(self, y):
        # Default implementation
        return y


class RandomPatchesClassifier(BaseRandomPatches, ClassifierMixin):
    """RandomPatches Classifier."""

    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 bootstrap=True,
                 oob_score=True,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):

        super(RandomPatchesClassifier, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose)

    def _set_oob_score(self, X, y):
        n_classes_ = self.n_classes_
        classes_ = self.classes_
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples, self.n_classes_))

        for estimator in self.estimators_:
            mask = np.ones(n_samples, dtype=np.bool)
            mask[estimator.indices_] = False

            try:
                predictions[mask, :] += estimator.predict_proba(X[mask, :])
            except:
                p = estimator.predict(X[mask, :])
                j = 0

                for i in xrange(n_samples):
                    if mask[i]:
                        predictions[i, p[j]] += 1
                        j += 1

        if (predictions.sum(axis=1) == 0).any():
            warn("Some inputs do not have OOB scores. "
                 "This probably means too few estimators were used "
                 "to compute any reliable oob estimates.")

        oob_decision_function = (predictions / predictions.sum(axis=1)[:, np.newaxis])
        oob_score = np.mean((y == classes_.take(
                np.argmax(predictions, axis=1),
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
            return self.classes_.take(np.argmax(self.predict_proba(X), axis=1), axis=0)

        except:
            # Check data
            X, = check_arrays(X)

            # Assign chunk of estimators to jobs
            n_jobs, n_estimators, starts = _partition_estimators(self)

            # Parallel loop
            all_counts = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
                delayed(_parallel_predict)(
                    self.estimators_[starts[i]:starts[i + 1]],
                    X,
                    self.n_classes_)
                for i in range(n_jobs))

            # Reduce
            counts = all_counts[0]

            for j in xrange(1, len(all_counts)):
                counts += all_counts[j]

            return self.classes_.take(np.argmax(counts, axis=1), axis=0)

    def predict_proba(self, X):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the mean predicted class probabilities of the estimators in the ensemble.

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

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, starts = _partition_estimators(self)

        # Parallel loop
        all_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_proba)(
                self.estimators_[starts[i]:starts[i + 1]],
                X,
                self.n_classes_)
            for i in range(n_jobs))

        # Reduce
        proba = all_proba[0]

        for j in xrange(1, len(all_proba)):
            proba += all_proba[j]

        proba /= self.n_estimators

        return proba

    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the mean predicted class log-probabilities of the estimators in the ensemble.

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
        return np.log(self.predict_proba(X))


class RandomPatchesRegressor(BaseRandomPatches, RegressorMixin):
    """RandomPatches Regressor."""

    def __init__(self,
                 base_estimator,
                 n_estimators=10,
                 bootstrap=True,
                 oob_score=True,
                 n_jobs=1,
                 random_state=None,
                 verbose=0):
        super(RandomPatchesRegressor, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            bootstrap=bootstrap,
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

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, starts = _partition_estimators(self)

        # Parallel loop
        all_y_hat = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]], X)
            for i in range(n_jobs))

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat

    def _set_oob_score(self, X, y):
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples,))
        n_predictions = np.zeros((n_samples,))

        for estimator in self.estimators_:
            mask = np.ones(n_samples, dtype=np.bool)
            mask[estimator.indices_] = False

            predictions[mask] = estimator.predict(X[mask, :])
            n_predictions[mask, :] += 1

        if (n_predictions == 0).any():
            warn("Some inputs do not have OOB scores. "
                 "This probably means too few estimators were used "
                 "to compute any reliable oob estimates.")
            n_predictions[n_predictions == 0] = 1

        predictions /= n_predictions

        self.oob_prediction_ = predictions
        self.oob_score_ = r2_score(y, predictions)
