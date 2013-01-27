"""Weight Boosting

This module contains weight boosting estimators for both classification and
regression.

The module structure is the following:

- The ``BaseAdaBoost`` base class implements a common ``fit`` method
  for all the estimators in the module. Regression and classification
  only differ from each other in the loss function that is optimized.

- ``AdaBoostClassifier`` implements adaptive boosting (AdaBoost-SAMME) for
  classification problems.

- ``AdaBoostRegressor`` implements adaptive boosting (AdaBoost.R2) for
  regression problems.
"""

# Authors: Noel Dawe, Gilles Louppe
# License: BSD Style

from abc import ABCMeta, abstractmethod

import numpy as np
from numpy.core.umath_tests import inner1d

from .base import BaseEnsemble
from ..base import ClassifierMixin, RegressorMixin
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor
from ..utils import check_arrays
from ..metrics import accuracy_score, r2_score


__all__ = [
    'AdaBoostClassifier',
    'AdaBoostRegressor',
]


class BaseWeightBoosting(BaseEnsemble):
    """Base class for AdaBoost estimators.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self,
                 base_estimator,
                 n_estimators=50,
                 estimator_params=tuple(),
                 learning_rate=0.5,
                 compute_importances=False):
        super(BaseWeightBoosting, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params)

        self.estimator_weights_ = None
        self.errors_ = None
        self.learning_rate = learning_rate
        self.compute_importances = compute_importances
        self.feature_importances_ = None

    def fit(self, X, y, sample_weight=None):
        """Build a boosted classifier/regressor from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights. If None, the sample weights are initialized to
            1 / n_samples.

        Returns
        -------
        self : object
            Returns self.
        """
        # Check parameters
        if self.learning_rate <= 0:
            raise ValueError("``learning_rate`` must be greater than zero")

        if self.compute_importances:
            self.base_estimator.set_params(compute_importances=True)

        # Check data
        X, y = check_arrays(X, y, sparse_format="dense")

        if sample_weight is None:
            # Initialize weights to 1 / n_samples
            sample_weight = np.empty(X.shape[0], dtype=np.float)
            sample_weight[:] = 1. / X.shape[0]
        else:
            # Normalize existing weights
            sample_weight = np.copy(sample_weight) / sample_weight.sum()

            # Check that the sample weights sum is positive
            if sample_weight.sum() <= 0:
                raise ValueError(
                    "Attempting to fit with a non-positive "
                    "weighted number of samples.")

        # Clear any previous fit results
        self.estimators_ = []
        self.estimator_weights_ = np.zeros(self.n_estimators, dtype=np.float)
        self.errors_ = np.ones(self.n_estimators, dtype=np.float)

        for iboost in xrange(self.n_estimators):
            # Boosting step
            sample_weight, weight, error = self._boost(
                iboost,
                X, y,
                sample_weight)

            # Early termination
            if sample_weight is None:
                break

            self.estimator_weights_[iboost] = weight
            self.errors_[iboost] = error

            # Stop if error is zero
            if error == 0:
                break

            # Stop if the sum of sample weights has become non-positive
            if np.sum(sample_weight) <= 0:
                break

            if iboost < self.n_estimators - 1:
                # Normalize
                sample_weight /= sample_weight.sum()

        # Sum the importances
        try:
            if self.compute_importances:
                norm = self.estimator_weights_.sum()
                self.feature_importances_ = (
                    sum(weight * clf.feature_importances_ for weight, clf
                        in zip(self.estimator_weights_, self.estimators_))
                    / norm)

        except AttributeError:
            raise AttributeError(
                "Unable to compute feature importances "
                "since base_estimator does not have a "
                "``feature_importances_`` attribute")

        return self

    @abstractmethod
    def _boost(self, iboost, X, y, sample_weight):
        """Implement a single boost.

        Warning: This method needs to be overriden by subclasses.

        Parameters
        ----------
        iboost : int
            The index of the current boost iteration.

        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (integers that correspond to classes).

        sample_weight : array-like of shape = [n_samples]
            The current sample weights.

        Returns
        -------
        sample_weight : array-like of shape = [n_samples] or None
            The reweighted sample weights.
            If None then boosting has terminated early.

        weight : float
            The weight for the current boost.
            If None then boosting has terminated early.

        error : float
            The classification error for the current boost.
            If None then boosting has terminated early.
        """
        pass

    def staged_score(self, X, y, n_estimators=-1):
        """Return staged scores for X, y.

        This generator method yields the ensemble score after each iteration of
        boosting and therefore allows monitoring, such as to determine the
        score on a test set after each boost.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training set.

        y : array-like, shape = [n_samples]
            Labels for X.

        Returns
        -------
        z : float

        """
        for y_pred in self.staged_predict(X, n_estimators=n_estimators):
            if isinstance(self, ClassifierMixin):
                yield accuracy_score(y, y_pred)
            else:
                yield r2_score(y, y_pred)


def _samme_proba(estimator, n_classes, X):
    """
    Calculate algorithm 4, step 2, equation c) of Zhu et al [1]

    References
    ----------
    .. [1] J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class AdaBoost", 2009.

    """
    proba = estimator.predict_proba(X)

    # Displace zero probabilities so the log is defined.
    # Also fix negative elements which may occur with
    # negative sample weights.
    proba[proba <= 0] = 1e-5
    log_proba = np.log(proba)

    return (n_classes - 1) * (log_proba - (1. / n_classes)
                           * log_proba.sum(axis=1)[:, np.newaxis])


class AdaBoostClassifier(BaseWeightBoosting, ClassifierMixin):
    """An AdaBoost classifier.

    An AdaBoost classifier is a meta-estimator that begins by fitting a
    classifier on the original dataset and then fits additional copies of the
    classifer on the same dataset but where the weights of incorrectly
    classified instances are adjusted such that subsequent classifiers focus
    more on difficult cases.

    This class implements the algorithm known as AdaBoost-SAMME [2].

    Parameters
    ----------
    base_estimator : object, optional (default=DecisionTreeClassifier)
        The base estimator from which the boosted ensemble is built.
        Support for sample weighting is required, as well as proper `classes_`
        and `n_classes_` attributes.

    n_estimators : integer, optional (default=50)
        The maximum number of estimators at which boosting is terminated.
        In case of perfect fit, the learning procedure is stopped early.

    learning_rate : float, optional (default=0.1)
        Learning rate shrinks the contribution of each classifier by
        ``learning_rate``. There is a trade-off between ``learning_rate`` and
        ``n_estimators``.

    algorithm : string, optional (default="SAMME.R")
        If "SAMME.R" then use the SAMME.R real boosting algorithm.
        ``base_estimator`` must support calculation of class probabilities.
        If "SAMME" then use the SAMME discrete boosting algorithm.

    compute_importances : boolean, optional (default=False)
        Whether feature importances are computed and stored in the
        ``feature_importances_`` attribute when calling fit.

    Attributes
    ----------
    `estimators_` : list of classifiers
        The collection of fitted sub-estimators.

    `classes_` : array of shape = [n_classes]
        The classes labels.

    `n_classes_` : int
        The number of classes.

    `estimator_weights_` : list of floats
        Weights for each estimator in the boosted ensemble.

    `errors_` : list of floats
        Classification error for each estimator in the boosted
        ensemble.

    `feature_importances_` : array of shape = [n_features]
        The feature importances if supported by the ``base_estimator``.
        Only computed if ``compute_importances=True``.

    See also
    --------
    AdaBoostRegressor, GradientBoostingClassifier, DecisionTreeClassifier

    References
    ----------
    .. [1] Y. Freund, R. Schapire, "A Decision-Theoretic Generalization of
           on-Line Learning and an Application to Boosting", 1995.

    .. [2] J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class AdaBoost", 2009.

    """
    def __init__(self,
                 base_estimator=DecisionTreeClassifier(max_depth=1),
                 n_estimators=50,
                 learning_rate=0.5,
                 algorithm="SAMME.R",
                 compute_importances=False):
        super(AdaBoostClassifier, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            compute_importances=compute_importances)

        self.algorithm = algorithm

    def fit(self, X, y, sample_weight=None):
        """Build a boosted classifier from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (integers that correspond to classes).

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights. If None, the sample weights are initialized to
            1 / n_samples.

        Returns
        -------
        self : object
            Returns self.
        """
        # Check that the base estimator is a classifier
        if not isinstance(self.base_estimator, ClassifierMixin):
            raise TypeError("``base_estimator`` must be a "
                            "subclass of ``ClassifierMixin``")

        # Check that algorithm is supported
        if self.algorithm != "SAMME" and self.algorithm != "SAMME.R":
            raise ValueError("``algorithm`` %s is not supported"
                             % self.algorithm)

        #  SAMME-R requires predict_proba-enabled base estimators
        if self.algorithm == "SAMME.R":
            if not hasattr(self.base_estimator, "predict_proba"):
                raise TypeError(
                    "The real AdaBoost algorithm requires that the weak"
                    "learner supports the calculation of class probabilities")

        return super(AdaBoostClassifier, self).fit(X, y, sample_weight)

    def _boost(self, iboost, X, y, sample_weight):
        """Implement a single boost.

        Perform a single boost according to the real multi-class SAMME.R
        algorithm or to the discrete SAMME algorithm and return the updated
        sample weights.

        Parameters
        ----------
        iboost : int
            The index of the current boost iteration.

        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (integers that correspond to classes).

        sample_weight : array-like of shape = [n_samples]
            The current sample weights.

        Returns
        -------
        sample_weight : array-like of shape = [n_samples] or None
            The reweighted sample weights.
            If None then boosting has terminated early.

        weight : float
            The weight for the current boost.
            If None then boosting has terminated early.

        error : float
            The classification error for the current boost.
            If None then boosting has terminated early.
        """
        if self.algorithm == "SAMME.R":
            return self._boost_real(iboost, X, y, sample_weight)
        else:  # elif self.algorithm == "SAMME":
            return self._boost_discrete(iboost, X, y, sample_weight)

    def _boost_real(self, iboost, X, y, sample_weight):
        """Implement a single boost using the SAMME.R real algorithm."""
        estimator = self._make_estimator()

        if hasattr(estimator, 'fit_predict_proba'):
            # Optimization for estimators that are able to save redundant
            # computations when calling fit + predict_proba
            # on the same input X
            y_predict_proba = estimator.fit_predict_proba(
                X, y, sample_weight=sample_weight)
        else:
            y_predict_proba = estimator.fit(
                X, y, sample_weight=sample_weight).predict_proba(X)

        if iboost == 0:
            self.classes_ = getattr(estimator, 'classes_', None)
            self.n_classes_ = getattr(estimator, 'n_classes_',
                                      getattr(estimator, 'n_classes', 1))

        y_predict = self.classes_.take(np.argmax(y_predict_proba, axis=1),
                                       axis=0)

        # Instances incorrectly classified
        incorrect = y_predict != y

        # Error fraction
        error = np.mean(np.average(incorrect, weights=sample_weight, axis=0))

        # Stop if classification is perfect
        if error == 0:
            return sample_weight, 1., 0.

        # Negative sample weights can yield an overall negative error...
        if error < 0:
            # use the absolute value
            # if you have a better idea of how to handle negative
            # sample weights let me know
            error = abs(error)

        # Construct y coding
        n_classes = self.n_classes_
        classes = self.classes_
        y_codes = np.array([-1. / (n_classes - 1), 1.])
        y_coding = y_codes.take(classes == y.reshape(y.shape[0], 1))

        # Displace zero probabilities so the log is defined.
        # Also fix negative elements which may occur with
        # negative sample weights.
        y_predict_proba[y_predict_proba <= 0] = 1e-5

        # Boost weight using multi-class AdaBoost SAMME.R alg
        weight = -1. * self.learning_rate * (
            ((n_classes - 1.) / n_classes) *
            inner1d(y_coding, np.log(y_predict_proba)))

        # Only boost the weights if I will fit again
        if not iboost == self.n_estimators - 1:
            # Only boost positive weights
            sample_weight *= np.exp(weight * ((sample_weight > 0) |
                                              (weight < 0)))

        return sample_weight, 1., error

    def _boost_discrete(self, iboost, X, y, sample_weight):
        """Implement a single boost using the SAMME discrete algorithm."""
        estimator = self._make_estimator()

        if hasattr(estimator, 'fit_predict'):
            # Optimization for estimators that are able to save redundant
            # computations when calling fit + predict
            # on the same input X
            y_predict = estimator.fit_predict(
                X, y, sample_weight=sample_weight)
        else:
            y_predict = estimator.fit(
                X, y, sample_weight=sample_weight).predict(X)

        if iboost == 0:
            self.classes_ = getattr(estimator, 'classes_', None)
            self.n_classes_ = getattr(estimator, 'n_classes_',
                                      getattr(estimator, 'n_classes', 1))

        # Instances incorrectly classified
        incorrect = y_predict != y

        # Error fraction
        error = np.mean(np.average(incorrect, weights=sample_weight, axis=0))

        # Stop if classification is perfect
        if error == 0:
            return sample_weight, 1., 0.

        # Negative sample weights can yield an overall negative error...
        if error < 0:
            # Use the absolute value
            error = abs(error)

        n_classes = self.n_classes_

        # Stop if the error is at least as bad as random guessing
        if error >= 1. - (1. / n_classes):
            self.estimators_.pop(-1)
            return None, None, None

        # Boost weight using multi-class AdaBoost SAMME alg
        weight = self.learning_rate * (
            np.log((1. - error) / error) +
            np.log(n_classes - 1.))

        # Only boost the weights if I will fit again
        if not iboost == self.n_estimators - 1:
            # Only boost positive weights
            sample_weight *= np.exp(weight * incorrect * ((sample_weight > 0) |
                                                          (weight < 0)))

        return sample_weight, weight, error

    def predict(self, X, n_estimators=-1):
        """Predict classes for X.

        The predicted class of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted classes.
        """
        pred = self.decision_function(X, n_estimators)

        if self.n_classes_ == 2:
            return self.classes_.take(pred > 0, axis=0)

        return self.classes_.take(np.argmax(pred, axis=1), axis=0)

    def staged_predict(self, X, n_estimators=-1):
        """Return staged predictions for X.

        The predicted class of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        This generator method yields the ensemble prediction after each
        iteration of boosting and therefore allows monitoring, such as to
        determine the prediction on a test set after each boost.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : generator of array, shape = [n_samples]
            The predicted classes.
        """
        n_classes = self.n_classes_
        classes = self.classes_

        if n_classes == 2:
            for pred in self.staged_decision_function(X, n_estimators):
                yield np.array(classes.take(pred > 0, axis=0))

        else:
            for pred in self.staged_decision_function(X, n_estimators):
                yield np.array(classes.take(
                    np.argmax(pred, axis=1), axis=0))

    def decision_function(self, X, n_estimators=-1):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        score : array, shape = [n_samples, k]
            The decision function of the input samples. Classes are
            ordered by arithmetical order. Binary classification is a
            special cases with ``k == 1``, otherwise ``k==n_classes``.
            For binary classification, values closer to -1 or 1 mean more
            like the first or second class in ``classes_``, respectively.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                ("{0} is not initialized. "
                 "Perform a fit first").format(self.__class__.__name__))

        n_classes = self.n_classes_
        classes = self.classes_[:, np.newaxis]
        pred = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):

            if i == n_estimators:
                break

            norm += weight

            if self.algorithm == "SAMME.R":
                current_pred = _samme_proba(estimator, n_classes, X)
            else:  # elif self.algorithm == "SAMME":
                current_pred = estimator.predict(X)
                current_pred = (current_pred == classes).T * weight

            if pred is None:
                pred = current_pred
            else:
                pred += current_pred

        pred /= norm
        if n_classes == 2:
            pred[:, 0] *= -1
            return pred.sum(axis=1)
        return pred

    def staged_decision_function(self, X, n_estimators=-1):
        """Compute decision function of ``X`` for each boosting iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each boosting iteration.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        score : generator of array, shape = [n_samples, k]
            The decision function of the input samples. Classes are
            ordered by arithmetical order. Binary classification is a
            special cases with ``k == 1``, otherwise ``k==n_classes``.
            For binary classification, values closer to -1 or 1 mean more
            like the first or second class in ``classes_``, respectively.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                ("{0} is not initialized. "
                 "Perform a fit first").format(self.__class__.__name__))

        n_classes = self.n_classes_
        classes = self.classes_[:, np.newaxis]
        pred = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):

            if i == n_estimators:
                break

            norm += weight

            if self.algorithm == "SAMME.R":
                current_pred = _samme_proba(estimator, n_classes, X)
            else:  # elif self.algorithm == "SAMME":
                current_pred = estimator.predict(X)
                current_pred = (current_pred == classes).T * weight

            if pred is None:
                pred = current_pred
            else:
                pred += current_pred

            if n_classes == 2:
                tmp_pred = np.copy(pred)
                tmp_pred[:, 0] *= -1
                yield (tmp_pred / norm).sum(axis=1)
            else:
                yield pred / norm

    def predict_proba(self, X, n_estimators=-1):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the weighted mean predicted class probabilities
        of the classifiers in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        n_classes = self.n_classes_
        proba = None

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):

            if i == n_estimators:
                break

            current_proba = _samme_proba(estimator, n_classes, X)

            if proba is None:
                proba = current_proba
            else:
                proba += current_proba

        proba = np.exp((1. / (n_classes - 1)) * proba)
        normalizer = proba.sum(axis=1)[:, np.newaxis]
        normalizer[normalizer == 0.0] = 1.0
        proba /= normalizer

        return proba

    def staged_predict_proba(self, X, n_estimators=-1):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the weighted mean predicted class probabilities
        of the classifiers in the ensemble.

        This generator method yields the ensemble predicted class probabilities
        after each iteration of boosting and therefore allows monitoring, such
        as to determine the predicted class probabilities on a test set after
        each boost.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        p : generator of array, shape = [n_samples]
            The class probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        n_classes = self.n_classes_
        proba = None

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):

            if i == n_estimators:
                break

            current_proba = _samme_proba(estimator, n_classes, X)

            if proba is None:
                proba = current_proba
            else:
                proba += current_proba

            real_proba = np.exp((1. / (n_classes - 1)) * proba)
            normalizer = real_proba.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            real_proba /= normalizer

            yield real_proba

    def predict_log_proba(self, X, n_estimators=-1):
        """Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the weighted mean predicted class log-probabilities
        of the classifiers in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        p : array of shape = [n_samples]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        return np.log(self.predict_proba(X, n_estimators=n_estimators))


class AdaBoostRegressor(BaseWeightBoosting, RegressorMixin):
    """An AdaBoost regressor.

    An AdaBoost regressor is a meta-estimator that begins by fitting a
    regressor on the original dataset and then fits additional copies of the
    regressor on the same dataset but where the weights of instances are
    adjusted according to the error of the current prediction. As such,
    subsequent regressors focus more on difficult cases.

    This class implements the algorithm known as AdaBoost.R2 [2].

    Parameters
    ----------
    base_estimator : object, optional (default=DecisionTreeRegressor)
        The base estimator from which the boosted ensemble is built.
        Support for sample weighting is required.

    n_estimators : integer, optional (default=50)
        The maximum number of estimators at which boosting is terminated.
        In case of perfect fit, the learning procedure is stopped early.

    learning_rate : float, optional (default=0.1)
        Learning rate shrinks the contribution of each regressor by
        ``learning_rate``. There is a trade-off between ``learning_rate`` and
        ``n_estimators``.

    compute_importances : boolean, optional (default=False)
        Whether feature importances are computed and stored in the
        ``feature_importances_`` attribute when calling fit.

    Attributes
    ----------
    `estimators_` : list of classifiers
        The collection of fitted sub-estimators.

    `estimator_weights_` : list of floats
        Weights for each estimator in the boosted ensemble.

    `errors_` : list of floats
        Regression error for each estimator in the boosted ensemble.

    `feature_importances_` : array of shape = [n_features]
        The feature importances if supported by the ``base_estimator``.
        Only computed if ``compute_importances=True``.

    See also
    --------
    AdaBoostClassifier, GradientBoostingRegressor, DecisionTreeRegressor

    References
    ----------
    .. [1] Y. Freund, R. Schapire, "A Decision-Theoretic Generalization of
           on-Line Learning and an Application to Boosting", 1995.

    .. [2] H. Drucker, "Improving Regressor using Boosting Techniques", 1997.

    """
    def __init__(self,
                 base_estimator=DecisionTreeRegressor(max_depth=3),
                 n_estimators=50,
                 learning_rate=0.1,
                 compute_importances=False):
        super(AdaBoostRegressor, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            compute_importances=compute_importances)

    def fit(self, X, y, sample_weight=None):
        """Build a boosted regressor from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (real numbers).

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights. If None, the sample weights are initialized to
            1 / n_samples.

        Returns
        -------
        self : object
            Returns self.
        """
        # Check that the base estimator is a regressor
        if not isinstance(self.base_estimator, RegressorMixin):
            raise TypeError("``base_estimator`` must be a "
                            "subclass of ``RegressorMixin``")

        # Fit
        return super(AdaBoostRegressor, self).fit(X, y, sample_weight)

    def _boost(self, iboost, X, y, sample_weight):
        """Implement a single boost for regression

        Perform a single boost according to the AdaBoost.R2 algorithm and
        return the updated sample weights.

        Parameters
        ----------
        iboost : int
            The index of the current boost iteration.

        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like of shape = [n_samples]
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        sample_weight : array-like of shape = [n_samples]
            The current sample weights.

        Returns
        -------
        sample_weight : array-like of shape = [n_samples] or None
            The reweighted sample weights.
            If None then boosting has terminated early.

        weight : float
            The weight for the current boost.
            If None then boosting has terminated early.

        error : float
            The regression error for the current boost.
            If None then boosting has terminated early.
        """
        estimator = self._make_estimator()

        if hasattr(estimator, 'fit_predict'):
            # Optimization for estimators that are able to save redundant
            # computations when calling fit + predict
            # on the same input X
            y_predict = estimator.fit_predict(
                X, y, sample_weight=sample_weight)
        else:
            y_predict = estimator.fit(
                X, y, sample_weight=sample_weight).predict(X)

        error_vect = np.abs(y_predict - y)
        error_max = error_vect.max()

        if error_max != 0.:
            error_vect /= error_vect.max()

        error = (sample_weight * error_vect).sum()

        # Stop if fit is perfect
        if error == 0:
            return sample_weight, 1., 0.

        # Negative sample weights can yield an overall negative error...
        if error < 0:
            # Use the absolute value
            error = abs(error)

        beta = error / (1. - error)

        # Boost weight using AdaBoost.R2 alg
        weight = self.learning_rate * np.log(1. / beta)

        if not iboost == self.n_estimators - 1:
            sample_weight *= np.power(
                beta,
                (1. - error_vect) * self.learning_rate)

        return sample_weight, weight, error

    def predict(self, X, n_estimators=-1):
        """Predict regression value for X.

        The predicted regression value of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted regression values.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                ("{0} is not initialized. "
                 "Perform a fit first").format(self.__class__.__name__))

        pred = None

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):
            if i == n_estimators:
                break

            current_pred = estimator.predict(X)

            if pred is None:
                pred = current_pred * weight
            else:
                pred += current_pred * weight

        pred /= self.estimator_weights_.sum()

        return pred

    def staged_predict(self, X, n_estimators=-1):
        """Return staged predictions for X.

        The predicted regression value of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        This generator method yields the ensemble prediction after each
        iteration of boosting and therefore allows monitoring, such as to
        determine the prediction on a test set after each boost.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter
            since it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : generator of array, shape = [n_samples]
            The predicted regression values.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                ("{0} is not initialized. "
                 "Perform a fit first").format(self.__class__.__name__))

        pred = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.estimator_weights_, self.estimators_)):
            if i == n_estimators:
                break

            current_pred = estimator.predict(X)

            if pred is None:
                pred = current_pred * weight
            else:
                pred += current_pred * weight

            norm += weight
            normed_pred = pred / norm

            yield normed_pred
