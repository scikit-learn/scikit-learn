"""Bagged ensemble of estimators."""

# Authors: Gilles Louppe
# License: BSD 3

import itertools
import numpy as np

from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..utils import check_random_state
from ..metrics import r2_score

from .base import BaseEnsemble

__all__ = ["BaggedClassifier",
           "BaggedRegressor"]

MAX_INT = np.iinfo(np.int32).max


def _parallel_build_estimators(n_estimators, ensemble, X, y, fit_params, seed):
    """Private function used to build a batch of estimators within a job."""
    random_state = check_random_state(seed)
    estimators = []

    for i in xrange(n_estimators):
        seed = random_state.randint(MAX_INT)

        estimator = ensemble._make_estimator(append=False)
        estimator.set_params(random_state=check_random_state(seed))

        if ensemble.bootstrap:
            n_samples = X.shape[0]
            indices = random_state.randint(0, n_samples, n_samples)
            estimator.fit(X[indices], y[indices], **fit_params)
            estimator.indices_ = indices

        else:
            estimator.fit(X, y)

        estimators.append(estimator)

    return estimators


def _parallel_predict_proba(estimators, X, n_classes):
    """Private function used to compute a batch of predictions within a job."""
    p = np.zeros((X.shape[0], n_classes))

    for estimator in estimators:
        if n_classes == estimator.n_classes_:
            p += estimator.predict_proba(X)

        else:
            proba = estimator.predict_proba(X)

            for j, c in enumerate(estimator.classes_):
                p[:, c] += proba[:, j]

    return p


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
    n_estimators = [ensemble.n_estimators / n_jobs] * n_jobs

    for i in xrange(ensemble.n_estimators % n_jobs):
        n_estimators[i] += 1

    starts = [0] * (n_jobs + 1)

    for i in xrange(1, n_jobs + 1):
        starts[i] = starts[i - 1] + n_estimators[i - 1]

    return n_jobs, n_estimators, starts


class BaseBagging(BaseEnsemble):
    """Base class for ensembles of estimators.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=True,
                       oob_score=True,
                       n_jobs=1,
                       random_state=None):
        super(BaseBagging, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params)

        self.bootstrap = bootstrap
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = check_random_state(random_state)

    def fit(self, X, y):
        """Build a ensemble of estimators from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        Returns
        -------
        self : object
            Returns self.
        """
        # Precompute some data
        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available"
                    " if bootstrap=True")

        if isinstance(self.base_estimator, ClassifierMixin):
            self.classes_ = np.unique(y)
            self.n_classes_ = len(self.classes_)
            y = np.searchsorted(self.classes_, y)

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, _ = _partition_estimators(self)

        # Parallel loop
        all_estimators = Parallel(n_jobs=n_jobs)(
            delayed(_parallel_build_estimators)(
                n_estimators[i],
                self,
                X,
                y,
                self._pre_fit(X, y),
                self.random_state.randint(MAX_INT))
            for i in xrange(n_jobs))

        # Reduce
        self.estimators_ = [e for e in itertools.chain(*all_estimators)]

        # Post-fit
        self._post_fit(X, y)

        # Calculate out of bag predictions and score
        if self.oob_score:
            if isinstance(self, ClassifierMixin):
                predictions = np.zeros((X.shape[0], self.n_classes_))

                for estimator in self.estimators_:
                    mask = np.ones(X.shape[0], dtype=np.bool)
                    mask[estimator.indices_] = False
                    predictions[mask, :] += estimator.predict_proba(X[mask, :])

                self.oob_decision_function_ = (predictions
                        / predictions.sum(axis=1)[:, np.newaxis])
                self.oob_score_ = np.mean(y == np.argmax(predictions, axis=1))

            else:
                # Regression:
                predictions = np.zeros(X.shape[0])
                n_predictions = np.zeros(X.shape[0])

                for estimator in self.estimators_:
                    mask = np.ones(X.shape[0], dtype=np.bool)
                    mask[estimator.indices_] = False
                    predictions[mask] += estimator.predict(X[mask, :])
                    n_predictions[mask] += 1

                predictions /= n_predictions

                self.oob_prediction_ = predictions
                self.oob_score_ = r2_score(y, predictions)

        return self

    def _pre_fit(self, X, y):
        return {}

    def _post_fit(self, X, y):
        pass


class BaggedClassifier(BaseBagging, ClassifierMixin):
    """Base class for ensemble of estimators-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=True,
                       oob_score=True,
                       n_jobs=1,
                       random_state=None):
        super(BaggedClassifier, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state)

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
        return self.classes_.take(
            np.argmax(self.predict_proba(X), axis=1),  axis=0)

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
        p : array of shape = [n_samples]
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        # Check data
        X = np.atleast_2d(X)

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, starts = _partition_estimators(self)

        # Parallel loop
        all_p = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_predict_proba)(
                self.estimators_[starts[i]:starts[i + 1]],
                X, self.n_classes_)
            for i in xrange(n_jobs))

        # Reduce
        p = sum(all_p) / self.n_estimators

        return p

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
        p : array of shape = [n_samples]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        return np.log(self.predict_proba(X))


class BaggedRegressor(BaseBagging, RegressorMixin):
    """Base class for ensemble of estimators-based regressors.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator,
                       n_estimators=10,
                       estimator_params=[],
                       bootstrap=True,
                       oob_score=True,
                       n_jobs=1,
                       random_state=None):
        super(BaggedRegressor, self).__init__(
            base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state)

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
        X = np.atleast_2d(X)

        # Assign chunk of estimators to jobs
        n_jobs, n_estimators, starts = _partition_estimators(self)

        # Parallel loop
        all_y_hat = Parallel(n_jobs=self.n_jobs)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i]:starts[i + 1]], X)
            for i in xrange(n_jobs))

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat
