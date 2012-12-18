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

import numpy as np

from .base import BaseEnsemble
from ..base import ClassifierMixin, RegressorMixin
from ..base import WeightedClassifierMixin, WeightedRegressorMixin
from ..metrics import weighted_r2_score
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor
from ..tree._tree import DTYPE
from ..utils import array2d, check_arrays


__all__ = [
    'AdaBoostClassifier',
    'AdaBoostRegressor',
]


class BaseWeightBoosting(BaseEnsemble):
    """Base class for weight boosting.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """
    def __init__(self, base_estimator=None,
                 n_estimators=10,
                 learning_rate=.5,
                 compute_importances=False):
        self.weights_ = []
        self.errors_ = []
        self.learning_rate = learning_rate
        self.compute_importances = compute_importances
        self.feature_importances_ = None

        if base_estimator is None:
            if isinstance(self, ClassifierMixin):
                base_estimator = DecisionTreeClassifier(max_depth=3)
            else:
                base_estimator = DecisionTreeRegressor(max_depth=3)
        elif (isinstance(self, ClassifierMixin)
              and not isinstance(base_estimator, ClassifierMixin)):
            raise TypeError("``base_estimator`` must be a "
                            "subclass of ``ClassifierMixin``")
        elif (isinstance(self, RegressorMixin)
              and not isinstance(base_estimator, RegressorMixin)):
            raise TypeError("``base_estimator`` must be a "
                            "subclass of ``RegressorMixin``")

        super(BaseWeightBoosting, self).__init__(
            base_estimator=base_estimator,
            n_estimators=n_estimators)

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
            Sample weights.

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
        X, y = check_arrays(X, y, sparse_format='dense')

        if sample_weight is None:
            # initialize weights to 1 / n_samples
            sample_weight = np.ones(X.shape[0], dtype=np.float64) / X.shape[0]
        else:
            sample_weight = np.copy(sample_weight)
            # normalize
            sample_weight /= sample_weight.sum()

        # Clear any previous fit results
        self.estimators_ = []
        self.weights_ = []
        self.errors_ = []

        for iboost in xrange(self.n_estimators):
            estimator = self._make_estimator()

            if hasattr(estimator, 'fit_predict'):
                # optim for estimators that are able to save redundant
                # computations when calling fit + predict
                # on the same input X
                p = estimator.fit_predict(X, y, sample_weight=sample_weight)
            else:
                p = estimator.fit(X, y, sample_weight=sample_weight).predict(X)

            if iboost == 0:
                if hasattr(estimator, 'classes_'):
                    self.classes_ = estimator.classes_
                    self.n_classes_ = estimator.n_classes_
                else:
                    self.n_classes_ = 1

            sample_weight, weight, error = self._boost(sample_weight, p, y,
                    iboost == self.n_estimators - 1)

            # early termination
            if sample_weight is None:
                break

            self.weights_.append(weight)
            self.errors_.append(error)

            if iboost < self.n_estimators - 1:
                # normalize
                sample_weight /= sample_weight.sum()

        # Sum the importances
        try:
            if self.compute_importances:
                norm = sum(self.weights_)
                self.feature_importances_ = (
                    sum(weight * clf.feature_importances_ for weight, clf
                      in zip(self.weights_, self.estimators_))
                    / norm)
        except AttributeError:
            raise AttributeError(
                    "Unable to compute feature importances "
                    "since base_estimator does not have a "
                    "``feature_importances_`` attribute")

        return self

    def predict(self, X, n_estimators=-1):
        """Predict class or regression value for X.

        The predicted class or regression value of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter since
            it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted classes.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                    ("{0} is not initialized. "
                     "Perform a fit first").format(self.__class__.__name__))

        X = array2d(X)
        n_samples, n_features = X.shape

        pred = None

        for i, (weight, estimator) in enumerate(
                zip(self.weights_, self.estimators_)):
            if i == n_estimators:
                break

            if isinstance(self, ClassifierMixin):
                current_pred = estimator.predict_proba(X) * weight
            else:
                current_pred = estimator.predict(X) * weight

            if pred is None:
                pred = current_pred
            else:
                pred += current_pred

        pred /= sum(self.weights_)

        if isinstance(self, ClassifierMixin):
            return self.classes_.take(np.argmax(pred, axis=1), axis=0)

        else:
            return pred

    def staged_predict(self, X, n_estimators=-1):
        """Return staged predictions for X.

        The predicted class or regression value of an input sample is computed
        as the weighted mean prediction of the classifiers in the ensemble.

        This generator method yields the ensemble prediction after each boost
        and therefore allows monitoring, such as to determine the prediction on
        a test set after each boost.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter since
            it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted classes.
        """
        if n_estimators == 0:
            raise ValueError("``n_estimators`` must not equal zero")

        if not self.estimators_:
            raise RuntimeError(
                    ("{0} is not initialized. "
                     "Perform a fit first").format(self.__class__.__name__))

        X = array2d(X)
        n_samples, n_features = X.shape

        pred = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.weights_, self.estimators_)):
            if i == n_estimators:
                break

            if isinstance(self, ClassifierMixin):
                current_pred = estimator.predict_proba(X) * weight
            else:
                current_pred = estimator.predict(X) * weight

            if pred is None:
                pred = current_pred
            else:
                pred += current_pred

            norm += weight
            normed_pred = pred / norm

            if isinstance(self, ClassifierMixin):
                yield self.classes_.take(np.argmax(normed_pred, axis=1),
                                         axis=0)
            else:
                yield normed_pred

    def staged_score(self, X, y, sample_weight=None, n_estimators=-1):
        """Return staged scores for X, y.

        This generator method yields the ensemble score after each boost
        and therefore allows monitoring, such as to determine the score on a
        test set after each boost.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training set.

        y : array-like, shape = [n_samples]
            Labels for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        z : float

        """
        for y_pred in self.staged_predict(X, n_estimators=n_estimators):
            if isinstance(self, ClassifierMixin):
                yield np.average((y_pred == y), weights=sample_weight)
            else:
                yield weighted_r2_score(y, y_pred, weights=sample_weight)


class AdaBoostClassifier(BaseWeightBoosting, WeightedClassifierMixin):
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
        and `n_classes_` attributes in case of classification.

    n_estimators : integer, optional (default=10)
        The maximum number of estimators at which boosting is terminated.

    learning_rate : float, optional (default=0.5)
        Learning rate shrinks the contribution of each classifier by
        ``learning_rate``. There is a trade-off between ``learning_rate`` and
        ``n_estimators``.

    compute_importances : boolean, optional (default=False)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    Attributes
    ----------
    `estimators_` : list of classifiers
        The collection of fitted sub-estimators.

    `classes_`: array of shape = [n_classes]
        The classes labels.

    `n_classes_`: int
        The number of classes.

    `weights_` : list of floats
        Weights for each estimator in the boosted ensemble.

    `errors_` : list of floats
        Classification error for each estimator in the boosted
        ensemble.

    `feature_importances_` : array of shape = [n_features]
        The feature importances if supported by the ``base_estimator``.

    See also
    --------
    AdaBoostRegressor, GradientBoostingClassifier, DecisionTreeClassifier

    References
    ----------

    .. [1] Yoav Freund, Robert E. Schapire. "A Decision-Theoretic
           Generalization of on-Line Learning and an Application
           to Boosting", 1995.

    .. [2] Ji Zhu, Hui Zou, Saharon Rosset, Trevor Hastie.
           "Multi-class AdaBoost", 2009.
    """
    def _boost(self, sample_weight, y_predict, y_true, is_last):
        """Implement a single boost

        Perform a single boost according to the multi-class SAMME algorithm and
        return the updated sample weights.

        Parameters
        ----------
        sample_weight : array-like of shape = [n_samples]
            The current sample weights.

        y_predict : array-like of shape = [n_samples]
            The predicted class labels.

        y_true : array-like of shape = [n_samples]
            The true class labels.

        is_last : Boolean
            True if this is the last boost.

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
        # instances incorrectly classified
        incorrect = y_predict != y_true

        # error fraction
        error = np.mean(np.average(incorrect, weights=sample_weight, axis=0))

        # stop if classification is perfect
        if error == 0:
            self.weights_.append(1.)
            self.errors_.append(error)
            return None, None, None

        # negative sample weights can yield an overall negative error...
        if error < 0:
            # use the absolute value
            # if you have a better idea of how to handle negative
            # sample weights let me know
            error = abs(error)

        n_classes = self.n_classes_

        # stop if the error is at least as bad as random guessing
        if error >= 1. - (1. / n_classes):
            self.estimators_.pop(-1)
            return None, None, None

        # boost weight using multi-class AdaBoost SAMME alg
        weight = self.learning_rate * (
                np.log((1. - error) / error) +
                np.log(n_classes - 1.))

        # only boost the weights if I will fit again
        if not is_last:
            sample_weight *= np.exp(weight * incorrect)

        return sample_weight, weight, error

    def predict_proba(self, X, n_estimators=-1):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the weighted mean predicted class probabilities
        of the classifiers in the ensemble.

        This method allows monitoring (i.e. determine error on testing set)
        after each boost. See examples/ensemble/plot_adaboost_error.py

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter since
            it is not necessary to fit separately for all choices of
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

        proba = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.weights_, self.estimators_)):

            if i == n_estimators:
                break

            current_proba = estimator.predict_proba(X) * weight

            if proba is None:
                proba = current_proba
            else:
                proba += current_proba

            norm += weight

        return proba / norm

    def staged_predict_proba(self, X, n_estimators=-1):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the weighted mean predicted class probabilities
        of the classifiers in the ensemble.

        This generator method yields the ensemble predicted class probabilities
        after each boost and therefore allows monitoring, such as to determine
        the predicted class probabilities on a test set after each boost.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        n_estimators : int, optional (default=-1)
            Use only the first ``n_estimators`` classifiers for the prediction.
            This is useful for grid searching the ``n_estimators`` parameter since
            it is not necessary to fit separately for all choices of
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

        proba = None
        norm = 0.

        for i, (weight, estimator) in enumerate(
                zip(self.weights_, self.estimators_)):
            if i == n_estimators:
                break

            current_proba = estimator.predict_proba(X) * weight

            if proba is None:
                proba = current_proba
            else:
                proba += current_proba

            norm += weight

            yield proba / norm

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
            This is useful for grid searching the ``n_estimators`` parameter since
            it is not necessary to fit separately for all choices of
            ``n_estimators``, but only the highest ``n_estimators``. Any
            negative value will result in all estimators being used.

        Returns
        -------
        p : array of shape = [n_samples]
            The class log-probabilities of the input samples. Classes are
            ordered by arithmetical order.
        """
        return np.log(self.predict_proba(X, n_estimators=n_estimators))


class AdaBoostRegressor(BaseWeightBoosting, WeightedRegressorMixin):
    """An AdaBoost regressor.

    An AdaBoost regressor is a meta-estimator that begins by fitting a
    regressor on the original dataset and then fits additional copies of the
    regressor on the same dataset but where the weights of instances are
    adjusted according to the error of the current prediction. As such,
    subsequent regressors focus more on difficult cases.

    This class implements the algorithm known as AdaBoost.R2 [2].

    Parameters
    ----------
    base_estimator : object, optional (default=DecisionTreeClassifier)
        The base estimator from which the boosted ensemble is built.
        Support for sample weighting is required, as well as proper `classes_`
        and `n_classes_` attributes in case of classification.

    n_estimators : integer, optional (default=10)
        The maximum number of estimators at which boosting is terminated.

    learning_rate : float, optional (default=0.5)
        Learning rate shrinks the contribution of each regressor by
        ``learning_rate``. There is a trade-off between ``learning_rate`` and
        ``n_estimators``.

    compute_importances : boolean, optional (default=False)
        Whether feature importances are computed and stored into the
        ``feature_importances_`` attribute when calling fit.

    Attributes
    ----------
    `estimators_` : list of classifiers
        The collection of fitted sub-estimators.

    `weights_` : list of floats
        Weights for each estimator in the boosted ensemble.

    `errors_` : list of floats
        Regression error for each estimator in the boosted ensemble.

    `feature_importances_` : array of shape = [n_features]
        The feature importances if supported by the ``base_estimator``.

    See also
    --------
    AdaBoostClassifier, GradientBoostingRegressor, DecisionTreeRegressor

    References
    ----------

    .. [1] Yoav Freund, Robert E. Schapire. "A Decision-Theoretic
           Generalization of on-Line Learning and an Application
           to Boosting", 1995.

    .. [2] Harris Drucker. "Improving Regressor using Boosting Techniques",
           1997.
    """
    def _boost(self, sample_weight, y_predict, y_true, is_last):
        """Implement a single boost for regression

        Perform a single boost according to the AdaBoost.R2 algorithm and
        return the updated sample weights.

        Parameters
        ----------
        sample_weight : array-like of shape = [n_samples]
            The current sample weights.

        y_predict : array-like of shape = [n_samples]
            The predicted class labels.

        y_true : array-like of shape = [n_samples]
            The true class labels.

        is_last : Boolean
            True if this is the last boost.

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
        error_vect = np.abs(y_predict - y_true)
        error_max = error_vect.max()

        if error_max != 0.:
            error_vect /= error_vect.max()

        error = (sample_weight * error_vect).sum()

        # stop if fit is perfect
        if error == 0:
            self.weights_.append(1.)
            self.errors_.append(error)
            return None, None, None

        # negative sample weights can yield an overall negative error...
        if error < 0:
            # use the absolute value
            # if you have a better idea of how to handle negative
            # sample weights let me know
            error = abs(error)

        beta = error / (1. - error)

        # boost weight using AdaBoost.R2 alg
        weight = self.learning_rate * np.log(1. / beta)
        if not is_last:
            sample_weight *= np.power(beta, (1. - error_vect) * self.learning_rate)

        return sample_weight, weight, error
