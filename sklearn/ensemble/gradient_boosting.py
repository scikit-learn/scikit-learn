"""Gradient Boosted Regression Trees

This module contains methods for fitting gradient boosted regression trees for
both classification and regression.

The module structure is the following:

- The ``BaseGradientBoosting`` base class implements a common ``fit`` method
  for all the estimators in the module. Regression and classification
  only differ in the concrete ``LossFunction`` used.

- ``GradientBoostingClassifier`` implements gradient boosting for
  classification problems.

- ``GradientBoostingRegressor`` implements gradient boosting for
  regression problems.
"""

# Authors: Peter Prettenhofer, Scott White, Gilles Louppe, Emanuele Olivetti,
#          Arnaud Joly, Jacob Schreiber
# License: BSD 3 clause

from abc import ABCMeta
from abc import abstractmethod
import warnings

from .base import BaseEnsemble
from ..base import ClassifierMixin
from ..base import RegressorMixin
from ..base import BaseEstimator
from ..base import is_classifier

from ._gradient_boosting import predict_stages
from ._gradient_boosting import predict_stage
from ._gradient_boosting import _random_sample_mask

import numbers
import numpy as np

from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import issparse
from scipy.special import expit

from time import time
from ..model_selection import train_test_split
from ..tree.tree import DecisionTreeRegressor
from ..tree._tree import DTYPE, DOUBLE
from ..tree._tree import TREE_LEAF
from . import _gb_losses

from ..utils import check_random_state
from ..utils import check_array
from ..utils import column_or_1d
from ..utils import check_consistent_length
from ..utils import deprecated
from ..utils.fixes import logsumexp
from ..utils.stats import _weighted_percentile
from ..utils.validation import check_is_fitted
from ..utils.multiclass import check_classification_targets
from ..exceptions import NotFittedError


# FIXME: 0.23
# All the losses and corresponding init estimators have been moved to the
# _losses module in 0.21. We deprecate them and keep them here for now in case
# someone has imported them. None of these losses can be used as a parameter
# to a GBDT estimator anyway (loss param only accepts strings).

@deprecated("QuantileEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class QuantileEstimator:
    """An estimator predicting the alpha-quantile of the training targets.

    Parameters
    ----------
    alpha : float
        The quantile
    """
    def __init__(self, alpha=0.9):
        if not 0 < alpha < 1.0:
            raise ValueError("`alpha` must be in (0, 1.0) but was %r" % alpha)
        self.alpha = alpha

    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : array, shape (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary

        sample_weight : numpy array of shape (n_samples,)
            Individual weights for each sample
        """
        if sample_weight is None:
            self.quantile = np.percentile(y, self.alpha * 100.0)
        else:
            self.quantile = _weighted_percentile(y, sample_weight,
                                                 self.alpha * 100.0)

    def predict(self, X):
        """Predict labels

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y : array, shape (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self)

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.quantile)
        return y


@deprecated("MeanEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class MeanEstimator:
    """An estimator predicting the mean of the training targets."""
    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : array, shape (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary

        sample_weight : numpy array of shape (n_samples,)
            Individual weights for each sample
        """
        if sample_weight is None:
            self.mean = np.mean(y)
        else:
            self.mean = np.average(y, weights=sample_weight)

    def predict(self, X):
        """Predict labels

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y : array, shape (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self)

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.mean)
        return y


@deprecated("LogOddsEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class LogOddsEstimator:
    """An estimator predicting the log odds ratio."""
    scale = 1.0

    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : array, shape (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary

        sample_weight : numpy array of shape (n_samples,)
            Individual weights for each sample
        """
        # pre-cond: pos, neg are encoded as 1, 0
        if sample_weight is None:
            pos = np.sum(y)
            neg = y.shape[0] - pos
        else:
            pos = np.sum(sample_weight * y)
            neg = np.sum(sample_weight * (1 - y))

        if neg == 0 or pos == 0:
            raise ValueError('y contains non binary labels.')
        self.prior = self.scale * np.log(pos / neg)

    def predict(self, X):
        """Predict labels

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y : array, shape (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self)

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.prior)
        return y


@deprecated("ScaledLogOddsEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class ScaledLogOddsEstimator(LogOddsEstimator):
    """Log odds ratio scaled by 0.5 -- for exponential loss. """
    scale = 0.5


@deprecated("PriorProbablityEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class PriorProbabilityEstimator:
    """An estimator predicting the probability of each
    class in the training data.
    """
    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : array, shape (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary

        sample_weight : array, shape (n_samples,)
            Individual weights for each sample
        """
        if sample_weight is None:
            sample_weight = np.ones_like(y, dtype=np.float64)
        class_counts = np.bincount(y, weights=sample_weight)
        self.priors = class_counts / class_counts.sum()

    def predict(self, X):
        """Predict labels

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y : array, shape (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self)

        y = np.empty((X.shape[0], self.priors.shape[0]), dtype=np.float64)
        y[:] = self.priors
        return y


@deprecated("Using ZeroEstimator is deprecated in version "
            "0.21 and will be removed in version 0.23.")
class ZeroEstimator:
    """An estimator that simply predicts zero.

    .. deprecated:: 0.21
        Using ``ZeroEstimator`` or ``init='zero'`` is deprecated in version
        0.21 and will be removed in version 0.23.

    """

    def fit(self, X, y, sample_weight=None):
        """Fit the estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data

        y : numpy, shape (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary

        sample_weight : array, shape (n_samples,)
            Individual weights for each sample
        """
        if np.issubdtype(y.dtype, np.signedinteger):
            # classification
            self.n_classes = np.unique(y).shape[0]
            if self.n_classes == 2:
                self.n_classes = 1
        else:
            # regression
            self.n_classes = 1

    def predict(self, X):
        """Predict labels

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y : array, shape (n_samples,)
            Returns predicted values.
        """
        check_is_fitted(self)

        y = np.empty((X.shape[0], self.n_classes), dtype=np.float64)
        y.fill(0.0)
        return y

    def predict_proba(self, X):
        return self.predict(X)


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class LossFunction(metaclass=ABCMeta):
    """Abstract base class for various loss functions.

    Parameters
    ----------
    n_classes : int
        Number of classes

    Attributes
    ----------
    K : int
        The number of regression trees to be induced;
        1 for regression and binary classification;
        ``n_classes`` for multi-class classification.
    """

    is_multi_class = False

    def __init__(self, n_classes):
        self.K = n_classes

    def init_estimator(self):
        """Default ``init`` estimator for loss function. """
        raise NotImplementedError()

    @abstractmethod
    def __call__(self, y, pred, sample_weight=None):
        """Compute the loss.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """

    @abstractmethod
    def negative_gradient(self, y, y_pred, **kargs):
        """Compute the negative gradient.

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        y_pred : array, shape (n_samples,)
            The predictions.
        """

    def update_terminal_regions(self, tree, X, y, residual, y_pred,
                                sample_weight, sample_mask,
                                learning_rate=0.1, k=0):
        """Update the terminal regions (=leaves) of the given tree and
        updates the current predictions of the model. Traverses tree
        and invokes template method `_update_terminal_region`.

        Parameters
        ----------
        tree : tree.Tree
            The tree object.
        X : array, shape (n, m)
            The data array.
        y : array, shape (n,)
            The target labels.
        residual : array, shape (n,)
            The residuals (usually the negative gradient).
        y_pred : array, shape (n,)
            The predictions.
        sample_weight : array, shape (n,)
            The weight of each sample.
        sample_mask : array, shape (n,)
            The sample mask to be used.
        learning_rate : float, default=0.1
            learning rate shrinks the contribution of each tree by
             ``learning_rate``.
        k : int, default 0
            The index of the estimator being updated.

        """
        # compute leaf for each sample in ``X``.
        terminal_regions = tree.apply(X)

        # mask all which are not in sample mask.
        masked_terminal_regions = terminal_regions.copy()
        masked_terminal_regions[~sample_mask] = -1

        # update each leaf (= perform line search)
        for leaf in np.where(tree.children_left == TREE_LEAF)[0]:
            self._update_terminal_region(tree, masked_terminal_regions,
                                         leaf, X, y, residual,
                                         y_pred[:, k], sample_weight)

        # update predictions (both in-bag and out-of-bag)
        y_pred[:, k] += (learning_rate
                         * tree.value[:, 0, 0].take(terminal_regions, axis=0))

    @abstractmethod
    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        """Template method for updating terminal regions (=leaves). """


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class RegressionLossFunction(LossFunction, metaclass=ABCMeta):
    """Base class for regression loss functions.

    Parameters
    ----------
    n_classes : int
        Number of classes
    """
    def __init__(self, n_classes):
        if n_classes != 1:
            raise ValueError("``n_classes`` must be 1 for regression but "
                             "was %r" % n_classes)
        super().__init__(n_classes)


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class LeastSquaresError(RegressionLossFunction):
    """Loss function for least squares (LS) estimation.
    Terminal regions need not to be updated for least squares.

    Parameters
    ----------
    n_classes : int
        Number of classes
    """

    def init_estimator(self):
        return MeanEstimator()

    def __call__(self, y, pred, sample_weight=None):
        """Compute the least squares loss.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        if sample_weight is None:
            return np.mean((y - pred.ravel()) ** 2.0)
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * ((y - pred.ravel()) ** 2.0)))

    def negative_gradient(self, y, pred, **kargs):
        """Compute the negative gradient.

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        pred : array, shape (n_samples,)
            The predictions.
        """
        return y - pred.ravel()

    def update_terminal_regions(self, tree, X, y, residual, y_pred,
                                sample_weight, sample_mask,
                                learning_rate=0.1, k=0):
        """Least squares does not need to update terminal regions.

        But it has to update the predictions.

        Parameters
        ----------
        tree : tree.Tree
            The tree object.
        X : array, shape (n, m)
            The data array.
        y : array, shape (n,)
            The target labels.
        residual : array, shape (n,)
            The residuals (usually the negative gradient).
        y_pred : array, shape (n,)
            The predictions.
        sample_weight : array, shape (n,)
            The weight of each sample.
        sample_mask : array, shape (n,)
            The sample mask to be used.
        learning_rate : float, default=0.1
            learning rate shrinks the contribution of each tree by
             ``learning_rate``.
        k : int, default 0
            The index of the estimator being updated.
        """
        # update predictions
        y_pred[:, k] += learning_rate * tree.predict(X).ravel()

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        pass


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class LeastAbsoluteError(RegressionLossFunction):
    """Loss function for least absolute deviation (LAD) regression.

    Parameters
    ----------
    n_classes : int
        Number of classes
    """
    def init_estimator(self):
        return QuantileEstimator(alpha=0.5)

    def __call__(self, y, pred, sample_weight=None):
        """Compute the least absolute error.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        if sample_weight is None:
            return np.abs(y - pred.ravel()).mean()
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * np.abs(y - pred.ravel())))

    def negative_gradient(self, y, pred, **kargs):
        """Compute the negative gradient.

        1.0 if y - pred > 0.0 else -1.0

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        pred : array, shape (n_samples,)
            The predictions.
        """
        pred = pred.ravel()
        return 2.0 * (y - pred > 0.0) - 1.0

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        """LAD updates terminal regions to median estimates. """
        terminal_region = np.where(terminal_regions == leaf)[0]
        sample_weight = sample_weight.take(terminal_region, axis=0)
        diff = y.take(terminal_region, axis=0) - pred.take(terminal_region, axis=0)
        tree.value[leaf, 0, 0] = _weighted_percentile(diff, sample_weight, percentile=50)


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class HuberLossFunction(RegressionLossFunction):
    """Huber loss function for robust regression.

    M-Regression proposed in Friedman 2001.

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.

    Parameters
    ----------
    n_classes : int
        Number of classes

    alpha : float
        Percentile at which to extract score
    """

    def __init__(self, n_classes, alpha=0.9):
        super().__init__(n_classes)
        self.alpha = alpha
        self.gamma = None

    def init_estimator(self):
        return QuantileEstimator(alpha=0.5)

    def __call__(self, y, pred, sample_weight=None):
        """Compute the Huber loss.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        pred = pred.ravel()
        diff = y - pred
        gamma = self.gamma
        if gamma is None:
            if sample_weight is None:
                gamma = np.percentile(np.abs(diff), self.alpha * 100)
            else:
                gamma = _weighted_percentile(np.abs(diff), sample_weight, self.alpha * 100)

        gamma_mask = np.abs(diff) <= gamma
        if sample_weight is None:
            sq_loss = np.sum(0.5 * diff[gamma_mask] ** 2.0)
            lin_loss = np.sum(gamma * (np.abs(diff[~gamma_mask]) - gamma / 2.0))
            loss = (sq_loss + lin_loss) / y.shape[0]
        else:
            sq_loss = np.sum(0.5 * sample_weight[gamma_mask] * diff[gamma_mask] ** 2.0)
            lin_loss = np.sum(gamma * sample_weight[~gamma_mask] *
                              (np.abs(diff[~gamma_mask]) - gamma / 2.0))
            loss = (sq_loss + lin_loss) / sample_weight.sum()
        return loss

    def negative_gradient(self, y, pred, sample_weight=None, **kargs):
        """Compute the negative gradient.

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        pred : array, shape (n_samples,)
            The predictions.

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        pred = pred.ravel()
        diff = y - pred
        if sample_weight is None:
            gamma = np.percentile(np.abs(diff), self.alpha * 100)
        else:
            gamma = _weighted_percentile(np.abs(diff), sample_weight, self.alpha * 100)
        gamma_mask = np.abs(diff) <= gamma
        residual = np.zeros((y.shape[0],), dtype=np.float64)
        residual[gamma_mask] = diff[gamma_mask]
        residual[~gamma_mask] = gamma * np.sign(diff[~gamma_mask])
        self.gamma = gamma
        return residual

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        terminal_region = np.where(terminal_regions == leaf)[0]
        sample_weight = sample_weight.take(terminal_region, axis=0)
        gamma = self.gamma
        diff = (y.take(terminal_region, axis=0)
                - pred.take(terminal_region, axis=0))
        median = _weighted_percentile(diff, sample_weight, percentile=50)
        diff_minus_median = diff - median
        tree.value[leaf, 0] = median + np.mean(
            np.sign(diff_minus_median) *
            np.minimum(np.abs(diff_minus_median), gamma))


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class QuantileLossFunction(RegressionLossFunction):
    """Loss function for quantile regression.

    Quantile regression allows to estimate the percentiles
    of the conditional distribution of the target.

    Parameters
    ----------
    n_classes : int
        Number of classes.

    alpha : float, optional (default = 0.9)
        The percentile
    """
    def __init__(self, n_classes, alpha=0.9):
        super().__init__(n_classes)
        self.alpha = alpha
        self.percentile = alpha * 100.0

    def init_estimator(self):
        return QuantileEstimator(self.alpha)

    def __call__(self, y, pred, sample_weight=None):
        """Compute the Quantile loss.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        pred = pred.ravel()
        diff = y - pred
        alpha = self.alpha

        mask = y > pred
        if sample_weight is None:
            loss = (alpha * diff[mask].sum() -
                    (1.0 - alpha) * diff[~mask].sum()) / y.shape[0]
        else:
            loss = ((alpha * np.sum(sample_weight[mask] * diff[mask]) -
                    (1.0 - alpha) * np.sum(sample_weight[~mask] * diff[~mask])) /
                    sample_weight.sum())
        return loss

    def negative_gradient(self, y, pred, **kargs):
        """Compute the negative gradient.

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        pred : array, shape (n_samples,)
            The predictions.
        """
        alpha = self.alpha
        pred = pred.ravel()
        mask = y > pred
        return (alpha * mask) - ((1.0 - alpha) * ~mask)

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        terminal_region = np.where(terminal_regions == leaf)[0]
        diff = (y.take(terminal_region, axis=0)
                - pred.take(terminal_region, axis=0))
        sample_weight = sample_weight.take(terminal_region, axis=0)

        val = _weighted_percentile(diff, sample_weight, self.percentile)
        tree.value[leaf, 0] = val


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class ClassificationLossFunction(LossFunction, metaclass=ABCMeta):
    """Base class for classification loss functions. """

    def _score_to_proba(self, score):
        """Template method to convert scores to probabilities.

         the does not support probabilities raises AttributeError.
        """
        raise TypeError('%s does not support predict_proba' % type(self).__name__)

    @abstractmethod
    def _score_to_decision(self, score):
        """Template method to convert scores to decisions.

        Returns int arrays.
        """


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class BinomialDeviance(ClassificationLossFunction):
    """Binomial deviance loss function for binary classification.

    Binary classification is a special case; here, we only need to
    fit one tree instead of ``n_classes`` trees.

    Parameters
    ----------
    n_classes : int
        Number of classes.
    """
    def __init__(self, n_classes):
        if n_classes != 2:
            raise ValueError("{0:s} requires 2 classes; got {1:d} class(es)"
                             .format(self.__class__.__name__, n_classes))
        # we only need to fit one tree for binary clf.
        super().__init__(1)

    def init_estimator(self):
        return LogOddsEstimator()

    def __call__(self, y, pred, sample_weight=None):
        """Compute the deviance (= 2 * negative log-likelihood).

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        # logaddexp(0, v) == log(1.0 + exp(v))
        pred = pred.ravel()
        if sample_weight is None:
            return -2.0 * np.mean((y * pred) - np.logaddexp(0.0, pred))
        else:
            return (-2.0 / sample_weight.sum() *
                    np.sum(sample_weight * ((y * pred) - np.logaddexp(0.0, pred))))

    def negative_gradient(self, y, pred, **kargs):
        """Compute the residual (= negative gradient).

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels
        """
        return y - expit(pred.ravel())

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        """Make a single Newton-Raphson step.

        our node estimate is given by:

            sum(w * (y - prob)) / sum(w * prob * (1 - prob))

        we take advantage that: y - prob = residual
        """
        terminal_region = np.where(terminal_regions == leaf)[0]
        residual = residual.take(terminal_region, axis=0)
        y = y.take(terminal_region, axis=0)
        sample_weight = sample_weight.take(terminal_region, axis=0)

        numerator = np.sum(sample_weight * residual)
        denominator = np.sum(sample_weight * (y - residual) * (1 - y + residual))

        # prevents overflow and division by zero
        if abs(denominator) < 1e-150:
            tree.value[leaf, 0, 0] = 0.0
        else:
            tree.value[leaf, 0, 0] = numerator / denominator

    def _score_to_proba(self, score):
        proba = np.ones((score.shape[0], 2), dtype=np.float64)
        proba[:, 1] = expit(score.ravel())
        proba[:, 0] -= proba[:, 1]
        return proba

    def _score_to_decision(self, score):
        proba = self._score_to_proba(score)
        return np.argmax(proba, axis=1)


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class MultinomialDeviance(ClassificationLossFunction):
    """Multinomial deviance loss function for multi-class classification.

    For multi-class classification we need to fit ``n_classes`` trees at
    each stage.

    Parameters
    ----------
    n_classes : int
        Number of classes
    """

    is_multi_class = True

    def __init__(self, n_classes):
        if n_classes < 3:
            raise ValueError("{0:s} requires more than 2 classes.".format(
                self.__class__.__name__))
        super().__init__(n_classes)

    def init_estimator(self):
        return PriorProbabilityEstimator()

    def __call__(self, y, pred, sample_weight=None):
        """Compute the Multinomial deviance.

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        # create one-hot label encoding
        Y = np.zeros((y.shape[0], self.K), dtype=np.float64)
        for k in range(self.K):
            Y[:, k] = y == k

        if sample_weight is None:
            return np.sum(-1 * (Y * pred).sum(axis=1) +
                          logsumexp(pred, axis=1))
        else:
            return np.sum(-1 * sample_weight * (Y * pred).sum(axis=1) +
                          logsumexp(pred, axis=1))

    def negative_gradient(self, y, pred, k=0, **kwargs):
        """Compute negative gradient for the ``k``-th class.

        Parameters
        ----------
        y : array, shape (n_samples,)
            The target labels.

        pred : array, shape (n_samples,)
            The predictions.

        k : int, optional (default=0)
            The index of the class
        """
        return y - np.nan_to_num(np.exp(pred[:, k] -
                                        logsumexp(pred, axis=1)))

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        """Make a single Newton-Raphson step. """
        terminal_region = np.where(terminal_regions == leaf)[0]
        residual = residual.take(terminal_region, axis=0)
        y = y.take(terminal_region, axis=0)
        sample_weight = sample_weight.take(terminal_region, axis=0)

        numerator = np.sum(sample_weight * residual)
        numerator *= (self.K - 1) / self.K

        denominator = np.sum(sample_weight * (y - residual) *
                             (1.0 - y + residual))

        # prevents overflow and division by zero
        if abs(denominator) < 1e-150:
            tree.value[leaf, 0, 0] = 0.0
        else:
            tree.value[leaf, 0, 0] = numerator / denominator

    def _score_to_proba(self, score):
        return np.nan_to_num(
            np.exp(score - (logsumexp(score, axis=1)[:, np.newaxis])))

    def _score_to_decision(self, score):
        proba = self._score_to_proba(score)
        return np.argmax(proba, axis=1)


@deprecated("All Losses in sklearn.ensemble.gradient_boosting are "
            "deprecated in version "
            "0.21 and will be removed in version 0.23.")
class ExponentialLoss(ClassificationLossFunction):
    """Exponential loss function for binary classification.

    Same loss as AdaBoost.

    References
    ----------
    Greg Ridgeway, Generalized Boosted Models: A guide to the gbm package, 2007

    Parameters
    ----------
    n_classes : int
        Number of classes.
    """
    def __init__(self, n_classes):
        if n_classes != 2:
            raise ValueError("{0:s} requires 2 classes; got {1:d} class(es)"
                             .format(self.__class__.__name__, n_classes))
        # we only need to fit one tree for binary clf.
        super().__init__(1)

    def init_estimator(self):
        return ScaledLogOddsEstimator()

    def __call__(self, y, pred, sample_weight=None):
        """Compute the exponential loss

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.
        """
        pred = pred.ravel()
        if sample_weight is None:
            return np.mean(np.exp(-(2. * y - 1.) * pred))
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * np.exp(-(2 * y - 1) * pred)))

    def negative_gradient(self, y, pred, **kargs):
        """Compute the residual (= negative gradient).

        Parameters
        ----------
        y : array, shape (n_samples,)
            True labels

        pred : array, shape (n_samples,)
            Predicted labels
        """
        y_ = -(2. * y - 1.)
        return y_ * np.exp(y_ * pred.ravel())

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        terminal_region = np.where(terminal_regions == leaf)[0]
        pred = pred.take(terminal_region, axis=0)
        y = y.take(terminal_region, axis=0)
        sample_weight = sample_weight.take(terminal_region, axis=0)

        y_ = 2. * y - 1.

        numerator = np.sum(y_ * sample_weight * np.exp(-y_ * pred))
        denominator = np.sum(sample_weight * np.exp(-y_ * pred))

        # prevents overflow and division by zero
        if abs(denominator) < 1e-150:
            tree.value[leaf, 0, 0] = 0.0
        else:
            tree.value[leaf, 0, 0] = numerator / denominator

    def _score_to_proba(self, score):
        proba = np.ones((score.shape[0], 2), dtype=np.float64)
        proba[:, 1] = expit(2.0 * score.ravel())
        proba[:, 0] -= proba[:, 1]
        return proba

    def _score_to_decision(self, score):
        return (score.ravel() >= 0.0).astype(np.int)


class VerboseReporter:
    """Reports verbose output to stdout.

    Parameters
    ----------
    verbose : int
        Verbosity level. If ``verbose==1`` output is printed once in a while
        (when iteration mod verbose_mod is zero).; if larger than 1 then output
        is printed for each update.
    """

    def __init__(self, verbose):
        self.verbose = verbose

    def init(self, est, begin_at_stage=0):
        """Initialize reporter

        Parameters
        ----------
        est : Estimator
            The estimator

        begin_at_stage : int
            stage at which to begin reporting
        """
        # header fields and line format str
        header_fields = ['Iter', 'Train Loss']
        verbose_fmt = ['{iter:>10d}', '{train_score:>16.4f}']
        # do oob?
        if est.subsample < 1:
            header_fields.append('OOB Improve')
            verbose_fmt.append('{oob_impr:>16.4f}')
        header_fields.append('Remaining Time')
        verbose_fmt.append('{remaining_time:>16s}')

        # print the header line
        print(('%10s ' + '%16s ' *
               (len(header_fields) - 1)) % tuple(header_fields))

        self.verbose_fmt = ' '.join(verbose_fmt)
        # plot verbose info each time i % verbose_mod == 0
        self.verbose_mod = 1
        self.start_time = time()
        self.begin_at_stage = begin_at_stage

    def update(self, j, est):
        """Update reporter with new iteration.

        Parameters
        ----------
        j : int
            The new iteration
        est : Estimator
            The estimator
        """
        do_oob = est.subsample < 1
        # we need to take into account if we fit additional estimators.
        i = j - self.begin_at_stage  # iteration relative to the start iter
        if (i + 1) % self.verbose_mod == 0:
            oob_impr = est.oob_improvement_[j] if do_oob else 0
            remaining_time = ((est.n_estimators - (j + 1)) *
                              (time() - self.start_time) / float(i + 1))
            if remaining_time > 60:
                remaining_time = '{0:.2f}m'.format(remaining_time / 60.0)
            else:
                remaining_time = '{0:.2f}s'.format(remaining_time)
            print(self.verbose_fmt.format(iter=j + 1,
                                          train_score=est.train_score_[j],
                                          oob_impr=oob_impr,
                                          remaining_time=remaining_time))
            if self.verbose == 1 and ((i + 1) // (self.verbose_mod * 10) > 0):
                # adjust verbose frequency (powers of 10)
                self.verbose_mod *= 10


class BaseGradientBoosting(BaseEnsemble, metaclass=ABCMeta):
    """Abstract base class for Gradient Boosting. """

    @abstractmethod
    def __init__(self, loss, learning_rate, n_estimators, criterion,
                 min_samples_split, min_samples_leaf, min_weight_fraction_leaf,
                 max_depth, min_impurity_decrease, min_impurity_split,
                 init, subsample, max_features, ccp_alpha,
                 random_state, alpha=0.9, verbose=0, max_leaf_nodes=None,
                 warm_start=False, presort='auto',
                 validation_fraction=0.1, n_iter_no_change=None,
                 tol=1e-4):

        self.n_estimators = n_estimators
        self.learning_rate = learning_rate
        self.loss = loss
        self.criterion = criterion
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.subsample = subsample
        self.max_features = max_features
        self.max_depth = max_depth
        self.min_impurity_decrease = min_impurity_decrease
        self.min_impurity_split = min_impurity_split
        self.ccp_alpha = ccp_alpha
        self.init = init
        self.random_state = random_state
        self.alpha = alpha
        self.verbose = verbose
        self.max_leaf_nodes = max_leaf_nodes
        self.warm_start = warm_start
        self.presort = presort
        self.validation_fraction = validation_fraction
        self.n_iter_no_change = n_iter_no_change
        self.tol = tol

    def _fit_stage(self, i, X, y, raw_predictions, sample_weight, sample_mask,
                   random_state, X_idx_sorted, X_csc=None, X_csr=None):
        """Fit another stage of ``n_classes_`` trees to the boosting model. """

        assert sample_mask.dtype == np.bool
        loss = self.loss_
        original_y = y

        # Need to pass a copy of raw_predictions to negative_gradient()
        # because raw_predictions is partially updated at the end of the loop
        # in update_terminal_regions(), and gradients need to be evaluated at
        # iteration i - 1.
        raw_predictions_copy = raw_predictions.copy()

        for k in range(loss.K):
            if loss.is_multi_class:
                y = np.array(original_y == k, dtype=np.float64)

            residual = loss.negative_gradient(y, raw_predictions_copy, k=k,
                                              sample_weight=sample_weight)

            # induce regression tree on residuals
            tree = DecisionTreeRegressor(
                criterion=self.criterion,
                splitter='best',
                max_depth=self.max_depth,
                min_samples_split=self.min_samples_split,
                min_samples_leaf=self.min_samples_leaf,
                min_weight_fraction_leaf=self.min_weight_fraction_leaf,
                min_impurity_decrease=self.min_impurity_decrease,
                min_impurity_split=self.min_impurity_split,
                max_features=self.max_features,
                max_leaf_nodes=self.max_leaf_nodes,
                random_state=random_state,
                presort=self.presort,
                ccp_alpha=self.ccp_alpha)

            if self.subsample < 1.0:
                # no inplace multiplication!
                sample_weight = sample_weight * sample_mask.astype(np.float64)

            X = X_csr if X_csr is not None else X
            tree.fit(X, residual, sample_weight=sample_weight,
                     check_input=False, X_idx_sorted=X_idx_sorted)

            # update tree leaves
            loss.update_terminal_regions(
                tree.tree_, X, y, residual, raw_predictions, sample_weight,
                sample_mask, learning_rate=self.learning_rate, k=k)

            # add tree to ensemble
            self.estimators_[i, k] = tree

        return raw_predictions

    def _check_params(self):
        """Check validity of parameters and raise ValueError if not valid. """
        if self.n_estimators <= 0:
            raise ValueError("n_estimators must be greater than 0 but "
                             "was %r" % self.n_estimators)

        if self.learning_rate <= 0.0:
            raise ValueError("learning_rate must be greater than 0 but "
                             "was %r" % self.learning_rate)

        if (self.loss not in self._SUPPORTED_LOSS
                or self.loss not in _gb_losses.LOSS_FUNCTIONS):
            raise ValueError("Loss '{0:s}' not supported. ".format(self.loss))

        if self.loss == 'deviance':
            loss_class = (_gb_losses.MultinomialDeviance
                          if len(self.classes_) > 2
                          else _gb_losses.BinomialDeviance)
        else:
            loss_class = _gb_losses.LOSS_FUNCTIONS[self.loss]

        if self.loss in ('huber', 'quantile'):
            self.loss_ = loss_class(self.n_classes_, self.alpha)
        else:
            self.loss_ = loss_class(self.n_classes_)

        if not (0.0 < self.subsample <= 1.0):
            raise ValueError("subsample must be in (0,1] but "
                             "was %r" % self.subsample)

        if self.init is not None:
            # init must be an estimator or 'zero'
            if isinstance(self.init, BaseEstimator):
                self.loss_.check_init_estimator(self.init)
            elif not (isinstance(self.init, str) and self.init == 'zero'):
                raise ValueError(
                    "The init parameter must be an estimator or 'zero'. "
                    "Got init={}".format(self.init)
                )

        if not (0.0 < self.alpha < 1.0):
            raise ValueError("alpha must be in (0.0, 1.0) but "
                             "was %r" % self.alpha)

        if isinstance(self.max_features, str):
            if self.max_features == "auto":
                # if is_classification
                if self.n_classes_ > 1:
                    max_features = max(1, int(np.sqrt(self.n_features_)))
                else:
                    # is regression
                    max_features = self.n_features_
            elif self.max_features == "sqrt":
                max_features = max(1, int(np.sqrt(self.n_features_)))
            elif self.max_features == "log2":
                max_features = max(1, int(np.log2(self.n_features_)))
            else:
                raise ValueError("Invalid value for max_features: %r. "
                                 "Allowed string values are 'auto', 'sqrt' "
                                 "or 'log2'." % self.max_features)
        elif self.max_features is None:
            max_features = self.n_features_
        elif isinstance(self.max_features, numbers.Integral):
            max_features = self.max_features
        else:  # float
            if 0. < self.max_features <= 1.:
                max_features = max(int(self.max_features *
                                       self.n_features_), 1)
            else:
                raise ValueError("max_features must be in (0, n_features]")

        self.max_features_ = max_features

        if not isinstance(self.n_iter_no_change,
                          (numbers.Integral, type(None))):
            raise ValueError("n_iter_no_change should either be None or an "
                             "integer. %r was passed"
                             % self.n_iter_no_change)

        allowed_presort = ('auto', True, False)
        if self.presort not in allowed_presort:
            raise ValueError("'presort' should be in {}. Got {!r} instead."
                             .format(allowed_presort, self.presort))

    def _init_state(self):
        """Initialize model state and allocate model state data structures. """

        self.init_ = self.init
        if self.init_ is None:
            self.init_ = self.loss_.init_estimator()

        self.estimators_ = np.empty((self.n_estimators, self.loss_.K),
                                    dtype=np.object)
        self.train_score_ = np.zeros((self.n_estimators,), dtype=np.float64)
        # do oob?
        if self.subsample < 1.0:
            self.oob_improvement_ = np.zeros((self.n_estimators),
                                             dtype=np.float64)

    def _clear_state(self):
        """Clear the state of the gradient boosting model. """
        if hasattr(self, 'estimators_'):
            self.estimators_ = np.empty((0, 0), dtype=np.object)
        if hasattr(self, 'train_score_'):
            del self.train_score_
        if hasattr(self, 'oob_improvement_'):
            del self.oob_improvement_
        if hasattr(self, 'init_'):
            del self.init_
        if hasattr(self, '_rng'):
            del self._rng

    def _resize_state(self):
        """Add additional ``n_estimators`` entries to all attributes. """
        # self.n_estimators is the number of additional est to fit
        total_n_estimators = self.n_estimators
        if total_n_estimators < self.estimators_.shape[0]:
            raise ValueError('resize with smaller n_estimators %d < %d' %
                             (total_n_estimators, self.estimators_[0]))

        self.estimators_ = np.resize(self.estimators_,
                                     (total_n_estimators, self.loss_.K))
        self.train_score_ = np.resize(self.train_score_, total_n_estimators)
        if (self.subsample < 1 or hasattr(self, 'oob_improvement_')):
            # if do oob resize arrays or create new if not available
            if hasattr(self, 'oob_improvement_'):
                self.oob_improvement_ = np.resize(self.oob_improvement_,
                                                  total_n_estimators)
            else:
                self.oob_improvement_ = np.zeros((total_n_estimators,),
                                                 dtype=np.float64)

    def _is_initialized(self):
        return len(getattr(self, 'estimators_', [])) > 0

    def _check_initialized(self):
        """Check that the estimator is initialized, raising an error if not."""
        check_is_fitted(self)

    def fit(self, X, y, sample_weight=None, monitor=None):
        """Fit the gradient boosting model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        y : array-like, shape (n_samples,)
            Target values (strings or integers in classification, real numbers
            in regression)
            For classification, labels must correspond to classes.

        sample_weight : array-like, shape (n_samples,) or None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. In the case of
            classification, splits are also ignored if they would result in any
            single class carrying a negative weight in either child node.

        monitor : callable, optional
            The monitor is called after each iteration with the current
            iteration, a reference to the estimator and the local variables of
            ``_fit_stages`` as keyword arguments ``callable(i, self,
            locals())``. If the callable returns ``True`` the fitting procedure
            is stopped. The monitor can be used for various things such as
            computing held-out estimates, early stopping, model introspect, and
            snapshoting.

        Returns
        -------
        self : object
        """
        # if not warmstart - clear the estimator state
        if not self.warm_start:
            self._clear_state()

        # Check input
        # Since check_array converts both X and y to the same dtype, but the
        # trees use different types for X and y, checking them separately.
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'], dtype=DTYPE)
        n_samples, self.n_features_ = X.shape

        sample_weight_is_none = sample_weight is None
        if sample_weight_is_none:
            sample_weight = np.ones(n_samples, dtype=np.float32)
        else:
            sample_weight = column_or_1d(sample_weight, warn=True)
            sample_weight_is_none = False

        check_consistent_length(X, y, sample_weight)

        y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)
        y = column_or_1d(y, warn=True)
        y = self._validate_y(y, sample_weight)

        if self.n_iter_no_change is not None:
            stratify = y if is_classifier(self) else None
            X, X_val, y, y_val, sample_weight, sample_weight_val = (
                train_test_split(X, y, sample_weight,
                                 random_state=self.random_state,
                                 test_size=self.validation_fraction,
                                 stratify=stratify))
            if is_classifier(self):
                if self.n_classes_ != np.unique(y).shape[0]:
                    # We choose to error here. The problem is that the init
                    # estimator would be trained on y, which has some missing
                    # classes now, so its predictions would not have the
                    # correct shape.
                    raise ValueError(
                        'The training data after the early stopping split '
                        'is missing some classes. Try using another random '
                        'seed.'
                    )
        else:
            X_val = y_val = sample_weight_val = None

        self._check_params()

        if not self._is_initialized():
            # init state
            self._init_state()

            # fit initial model and initialize raw predictions
            if self.init_ == 'zero':
                raw_predictions = np.zeros(shape=(X.shape[0], self.loss_.K),
                                           dtype=np.float64)
            else:
                # XXX clean this once we have a support_sample_weight tag
                if sample_weight_is_none:
                    self.init_.fit(X, y)
                else:
                    msg = ("The initial estimator {} does not support sample "
                           "weights.".format(self.init_.__class__.__name__))
                    try:
                        self.init_.fit(X, y, sample_weight=sample_weight)
                    except TypeError:  # regular estimator without SW support
                        raise ValueError(msg)
                    except ValueError as e:
                        if "pass parameters to specific steps of "\
                           "your pipeline using the "\
                           "stepname__parameter" in str(e):  # pipeline
                            raise ValueError(msg) from e
                        else:  # regular estimator whose input checking failed
                            raise

                raw_predictions = \
                    self.loss_.get_init_raw_predictions(X, self.init_)

            begin_at_stage = 0

            # The rng state must be preserved if warm_start is True
            self._rng = check_random_state(self.random_state)

        else:
            # add more estimators to fitted model
            # invariant: warm_start = True
            if self.n_estimators < self.estimators_.shape[0]:
                raise ValueError('n_estimators=%d must be larger or equal to '
                                 'estimators_.shape[0]=%d when '
                                 'warm_start==True'
                                 % (self.n_estimators,
                                    self.estimators_.shape[0]))
            begin_at_stage = self.estimators_.shape[0]
            # The requirements of _decision_function (called in two lines
            # below) are more constrained than fit. It accepts only CSR
            # matrices.
            X = check_array(X, dtype=DTYPE, order="C", accept_sparse='csr')
            raw_predictions = self._raw_predict(X)
            self._resize_state()

        if self.presort is True and issparse(X):
            raise ValueError(
                "Presorting is not supported for sparse matrices.")

        presort = self.presort
        # Allow presort to be 'auto', which means True if the dataset is dense,
        # otherwise it will be False.
        if presort == 'auto':
            presort = not issparse(X)

        X_idx_sorted = None
        if presort:
            X_idx_sorted = np.asfortranarray(np.argsort(X, axis=0),
                                             dtype=np.int32)

        # fit the boosting stages
        n_stages = self._fit_stages(
            X, y, raw_predictions, sample_weight, self._rng, X_val, y_val,
            sample_weight_val, begin_at_stage, monitor, X_idx_sorted)

        # change shape of arrays after fit (early-stopping or additional ests)
        if n_stages != self.estimators_.shape[0]:
            self.estimators_ = self.estimators_[:n_stages]
            self.train_score_ = self.train_score_[:n_stages]
            if hasattr(self, 'oob_improvement_'):
                self.oob_improvement_ = self.oob_improvement_[:n_stages]

        self.n_estimators_ = n_stages
        return self

    def _fit_stages(self, X, y, raw_predictions, sample_weight, random_state,
                    X_val, y_val, sample_weight_val,
                    begin_at_stage=0, monitor=None, X_idx_sorted=None):
        """Iteratively fits the stages.

        For each stage it computes the progress (OOB, train score)
        and delegates to ``_fit_stage``.
        Returns the number of stages fit; might differ from ``n_estimators``
        due to early stopping.
        """
        n_samples = X.shape[0]
        do_oob = self.subsample < 1.0
        sample_mask = np.ones((n_samples, ), dtype=np.bool)
        n_inbag = max(1, int(self.subsample * n_samples))
        loss_ = self.loss_

        if self.verbose:
            verbose_reporter = VerboseReporter(self.verbose)
            verbose_reporter.init(self, begin_at_stage)

        X_csc = csc_matrix(X) if issparse(X) else None
        X_csr = csr_matrix(X) if issparse(X) else None

        if self.n_iter_no_change is not None:
            loss_history = np.full(self.n_iter_no_change, np.inf)
            # We create a generator to get the predictions for X_val after
            # the addition of each successive stage
            y_val_pred_iter = self._staged_raw_predict(X_val)

        # perform boosting iterations
        i = begin_at_stage
        for i in range(begin_at_stage, self.n_estimators):

            # subsampling
            if do_oob:
                sample_mask = _random_sample_mask(n_samples, n_inbag,
                                                  random_state)
                # OOB score before adding this stage
                old_oob_score = loss_(y[~sample_mask],
                                      raw_predictions[~sample_mask],
                                      sample_weight[~sample_mask])

            # fit next stage of trees
            raw_predictions = self._fit_stage(
                i, X, y, raw_predictions, sample_weight, sample_mask,
                random_state, X_idx_sorted, X_csc, X_csr)

            # track deviance (= loss)
            if do_oob:
                self.train_score_[i] = loss_(y[sample_mask],
                                             raw_predictions[sample_mask],
                                             sample_weight[sample_mask])
                self.oob_improvement_[i] = (
                    old_oob_score - loss_(y[~sample_mask],
                                          raw_predictions[~sample_mask],
                                          sample_weight[~sample_mask]))
            else:
                # no need to fancy index w/ no subsampling
                self.train_score_[i] = loss_(y, raw_predictions, sample_weight)

            if self.verbose > 0:
                verbose_reporter.update(i, self)

            if monitor is not None:
                early_stopping = monitor(i, self, locals())
                if early_stopping:
                    break

            # We also provide an early stopping based on the score from
            # validation set (X_val, y_val), if n_iter_no_change is set
            if self.n_iter_no_change is not None:
                # By calling next(y_val_pred_iter), we get the predictions
                # for X_val after the addition of the current stage
                validation_loss = loss_(y_val, next(y_val_pred_iter),
                                        sample_weight_val)

                # Require validation_score to be better (less) than at least
                # one of the last n_iter_no_change evaluations
                if np.any(validation_loss + self.tol < loss_history):
                    loss_history[i % len(loss_history)] = validation_loss
                else:
                    break

        return i + 1

    def _make_estimator(self, append=True):
        # we don't need _make_estimator
        raise NotImplementedError()

    def _raw_predict_init(self, X):
        """Check input and compute raw predictions of the init estimtor."""
        self._check_initialized()
        X = self.estimators_[0, 0]._validate_X_predict(X, check_input=True)
        if X.shape[1] != self.n_features_:
            raise ValueError("X.shape[1] should be {0:d}, not {1:d}.".format(
                self.n_features_, X.shape[1]))
        if self.init_ == 'zero':
            raw_predictions = np.zeros(shape=(X.shape[0], self.loss_.K),
                                       dtype=np.float64)
        else:
            raw_predictions = self.loss_.get_init_raw_predictions(
                X, self.init_).astype(np.float64)
        return raw_predictions

    def _raw_predict(self, X):
        """Return the sum of the trees raw predictions (+ init estimator)."""
        raw_predictions = self._raw_predict_init(X)
        predict_stages(self.estimators_, X, self.learning_rate,
                       raw_predictions)
        return raw_predictions

    def _staged_raw_predict(self, X):
        """Compute raw predictions of ``X`` for each iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        raw_predictions : generator of array, shape (n_samples, k)
            The raw predictions of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
            Regression and binary classification are special cases with
            ``k == 1``, otherwise ``k==n_classes``.
        """
        X = check_array(X, dtype=DTYPE, order="C", accept_sparse='csr')
        raw_predictions = self._raw_predict_init(X)
        for i in range(self.estimators_.shape[0]):
            predict_stage(self.estimators_, i, X, self.learning_rate,
                          raw_predictions)
            yield raw_predictions.copy()

    @property
    def feature_importances_(self):
        """Return the feature importances (the higher, the more important the
           feature).

        Returns
        -------
        feature_importances_ : array, shape (n_features,)
            The values of this array sum to 1, unless all trees are single node
            trees consisting of only the root node, in which case it will be an
            array of zeros.
        """
        self._check_initialized()

        relevant_trees = [tree
                          for stage in self.estimators_ for tree in stage
                          if tree.tree_.node_count > 1]
        if not relevant_trees:
            # degenerate case where all trees have only one node
            return np.zeros(shape=self.n_features_, dtype=np.float64)

        relevant_feature_importances = [
            tree.tree_.compute_feature_importances(normalize=False)
            for tree in relevant_trees
        ]
        avg_feature_importances = np.mean(relevant_feature_importances,
                                          axis=0, dtype=np.float64)
        return avg_feature_importances / np.sum(avg_feature_importances)

    def _compute_partial_dependence_recursion(self, grid, target_features):
        """Fast partial dependence computation.

        Parameters
        ----------
        grid : ndarray, shape (n_samples, n_target_features)
            The grid points on which the partial dependence should be
            evaluated.
        target_features : ndarray, shape (n_target_features)
            The set of target features for which the partial dependence
            should be evaluated.

        Returns
        -------
        averaged_predictions : ndarray, shape \
                (n_trees_per_iteration, n_samples)
            The value of the partial dependence function on each grid point.
        """
        check_is_fitted(self,
                        msg="'estimator' parameter must be a fitted estimator")
        if self.init is not None:
            warnings.warn(
                'Using recursion method with a non-constant init predictor '
                'will lead to incorrect partial dependence values. '
                'Got init=%s.' % self.init,
                UserWarning
            )
        grid = np.asarray(grid, dtype=DTYPE, order='C')
        n_estimators, n_trees_per_stage = self.estimators_.shape
        averaged_predictions = np.zeros((n_trees_per_stage, grid.shape[0]),
                                        dtype=np.float64, order='C')
        for stage in range(n_estimators):
            for k in range(n_trees_per_stage):
                tree = self.estimators_[stage, k].tree_
                tree.compute_partial_dependence(grid, target_features,
                                                averaged_predictions[k])
        averaged_predictions *= self.learning_rate

        return averaged_predictions

    def _validate_y(self, y, sample_weight):
        # 'sample_weight' is not utilised but is used for
        # consistency with similar method _validate_y of GBC
        self.n_classes_ = 1
        if y.dtype.kind == 'O':
            y = y.astype(DOUBLE)
        # Default implementation
        return y

    def apply(self, X):
        """Apply trees in the ensemble to X, return leaf indices.

        .. versionadded:: 0.17

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will
            be converted to a sparse ``csr_matrix``.

        Returns
        -------
        X_leaves : array-like, shape (n_samples, n_estimators, n_classes)
            For each datapoint x in X and for each tree in the ensemble,
            return the index of the leaf x ends up in each estimator.
            In the case of binary classification n_classes is 1.
        """

        self._check_initialized()
        X = self.estimators_[0, 0]._validate_X_predict(X, check_input=True)

        # n_classes will be equal to 1 in the binary classification or the
        # regression case.
        n_estimators, n_classes = self.estimators_.shape
        leaves = np.zeros((X.shape[0], n_estimators, n_classes))

        for i in range(n_estimators):
            for j in range(n_classes):
                estimator = self.estimators_[i, j]
                leaves[:, i, j] = estimator.apply(X, check_input=False)

        return leaves


class GradientBoostingClassifier(ClassifierMixin, BaseGradientBoosting):
    """Gradient Boosting for classification.

    GB builds an additive model in a
    forward stage-wise fashion; it allows for the optimization of
    arbitrary differentiable loss functions. In each stage ``n_classes_``
    regression trees are fit on the negative gradient of the
    binomial or multinomial deviance loss function. Binary classification
    is a special case where only a single regression tree is induced.

    Read more in the :ref:`User Guide <gradient_boosting>`.

    Parameters
    ----------
    loss : {'deviance', 'exponential'}, optional (default='deviance')
        loss function to be optimized. 'deviance' refers to
        deviance (= logistic regression) for classification
        with probabilistic outputs. For loss 'exponential' gradient
        boosting recovers the AdaBoost algorithm.

    learning_rate : float, optional (default=0.1)
        learning rate shrinks the contribution of each tree by `learning_rate`.
        There is a trade-off between learning_rate and n_estimators.

    n_estimators : int (default=100)
        The number of boosting stages to perform. Gradient boosting
        is fairly robust to over-fitting so a large number usually
        results in better performance.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_estimators`.
        Choosing `subsample < 1.0` leads to a reduction of variance
        and an increase in bias.

    criterion : string, optional (default="friedman_mse")
        The function to measure the quality of a split. Supported criteria
        are "friedman_mse" for the mean squared error with improvement
        score by Friedman, "mse" for mean squared error, and "mae" for
        the mean absolute error. The default value of "friedman_mse" is
        generally the best as it can provide a better approximation in
        some cases.

        .. versionadded:: 0.18

    min_samples_split : int, float, optional (default=2)
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

    min_impurity_decrease : float, optional (default=0.)
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    min_impurity_split : float, (default=1e-7)
        Threshold for early stopping in tree growth. A node will split
        if its impurity is above the threshold, otherwise it is a leaf.

        .. deprecated:: 0.19
           ``min_impurity_split`` has been deprecated in favor of
           ``min_impurity_decrease`` in 0.19. The default value of
           ``min_impurity_split`` will change from 1e-7 to 0 in 0.23 and it
           will be removed in 0.25. Use ``min_impurity_decrease`` instead.

    init : estimator or 'zero', optional (default=None)
        An estimator object that is used to compute the initial predictions.
        ``init`` has to provide :meth:`fit` and :meth:`predict_proba`. If
        'zero', the initial raw predictions are set to zero. By default, a
        ``DummyEstimator`` predicting the classes priors is used.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `int(max_features * n_features)` features are considered at each
          split.
        - If "auto", then `max_features=sqrt(n_features)`.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.

        Choosing `max_features < n_features` leads to a reduction of variance
        and an increase in bias.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    verbose : int, default: 0
        Enable verbose output. If 1 then it prints progress and performance
        once in a while (the more trees the lower the frequency). If greater
        than 1 then it prints progress and performance for every tree.

    max_leaf_nodes : int or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    warm_start : bool, default: False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just erase the
        previous solution. See :term:`the Glossary <warm_start>`.

    presort : bool or 'auto', optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error.

        .. versionadded:: 0.17
           *presort* parameter.

    validation_fraction : float, optional, default 0.1
        The proportion of training data to set aside as validation set for
        early stopping. Must be between 0 and 1.
        Only used if ``n_iter_no_change`` is set to an integer.

        .. versionadded:: 0.20

    n_iter_no_change : int, default None
        ``n_iter_no_change`` is used to decide if early stopping will be used
        to terminate training when validation score is not improving. By
        default it is set to None to disable early stopping. If set to a
        number, it will set aside ``validation_fraction`` size of the training
        data as validation and terminate training when validation score is not
        improving in all of the previous ``n_iter_no_change`` numbers of
        iterations. The split is stratified.

        .. versionadded:: 0.20

    tol : float, optional, default 1e-4
        Tolerance for the early stopping. When the loss is not improving
        by at least tol for ``n_iter_no_change`` iterations (if set to a
        number), the training stops.

        .. versionadded:: 0.20

    ccp_alpha : non-negative float, optional (default=0.0)
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    Attributes
    ----------
    n_estimators_ : int
        The number of estimators as selected by early stopping (if
        ``n_iter_no_change`` is specified). Otherwise it is set to
        ``n_estimators``.

        .. versionadded:: 0.20

    feature_importances_ : array, shape (n_features,)
        The feature importances (the higher, the more important the feature).

    oob_improvement_ : array, shape (n_estimators,)
        The improvement in loss (= deviance) on the out-of-bag samples
        relative to the previous iteration.
        ``oob_improvement_[0]`` is the improvement in
        loss of the first stage over the ``init`` estimator.

    train_score_ : array, shape (n_estimators,)
        The i-th score ``train_score_[i]`` is the deviance (= loss) of the
        model at iteration ``i`` on the in-bag sample.
        If ``subsample == 1`` this is the deviance on the training data.

    loss_ : LossFunction
        The concrete ``LossFunction`` object.

    init_ : estimator
        The estimator that provides the initial predictions.
        Set via the ``init`` argument or ``loss.init_estimator``.

    estimators_ : ndarray of DecisionTreeRegressor,\
shape (n_estimators, ``loss_.K``)
        The collection of fitted sub-estimators. ``loss_.K`` is 1 for binary
        classification, otherwise n_classes.

    classes_ : array of shape = [n_classes]
        The classes labels.

    Notes
    -----
    The features are always randomly permuted at each split. Therefore,
    the best found split may vary, even with the same training data and
    ``max_features=n_features``, if the improvement of the criterion is
    identical for several splits enumerated during the search of the best
    split. To obtain a deterministic behaviour during fitting,
    ``random_state`` has to be fixed.

    See also
    --------
    sklearn.ensemble.HistGradientBoostingClassifier,
    sklearn.tree.DecisionTreeClassifier, RandomForestClassifier
    AdaBoostClassifier

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.

    J. Friedman, Stochastic Gradient Boosting, 1999

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning Ed. 2, Springer, 2009.
    """

    _SUPPORTED_LOSS = ('deviance', 'exponential')

    def __init__(self, loss='deviance', learning_rate=0.1, n_estimators=100,
                 subsample=1.0, criterion='friedman_mse', min_samples_split=2,
                 min_samples_leaf=1, min_weight_fraction_leaf=0.,
                 max_depth=3, min_impurity_decrease=0.,
                 min_impurity_split=None, init=None,
                 random_state=None, max_features=None, verbose=0,
                 max_leaf_nodes=None, warm_start=False,
                 presort='auto', validation_fraction=0.1,
                 n_iter_no_change=None, tol=1e-4, ccp_alpha=0.0):

        super().__init__(
            loss=loss, learning_rate=learning_rate, n_estimators=n_estimators,
            criterion=criterion, min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_depth=max_depth, init=init, subsample=subsample,
            max_features=max_features,
            random_state=random_state, verbose=verbose,
            max_leaf_nodes=max_leaf_nodes,
            min_impurity_decrease=min_impurity_decrease,
            min_impurity_split=min_impurity_split,
            warm_start=warm_start, presort=presort,
            validation_fraction=validation_fraction,
            n_iter_no_change=n_iter_no_change, tol=tol, ccp_alpha=ccp_alpha)

    def _validate_y(self, y, sample_weight):
        check_classification_targets(y)
        self.classes_, y = np.unique(y, return_inverse=True)
        n_trim_classes = np.count_nonzero(np.bincount(y, sample_weight))
        if n_trim_classes < 2:
            raise ValueError("y contains %d class after sample_weight "
                             "trimmed classes with zero weights, while a "
                             "minimum of 2 classes are required."
                             % n_trim_classes)
        self.n_classes_ = len(self.classes_)
        return y

    def decision_function(self, X):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        score : array, shape (n_samples, n_classes) or (n_samples,)
            The decision function of the input samples, which corresponds to
            the raw values predicted from the trees of the ensemble . The
            order of the classes corresponds to that in the attribute
            :term:`classes_`. Regression and binary classification produce an
            array of shape [n_samples].
        """
        X = check_array(X, dtype=DTYPE, order="C", accept_sparse='csr')
        raw_predictions = self._raw_predict(X)
        if raw_predictions.shape[1] == 1:
            return raw_predictions.ravel()
        return raw_predictions

    def staged_decision_function(self, X):
        """Compute decision function of ``X`` for each iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        score : generator of array, shape (n_samples, k)
            The decision function of the input samples, which corresponds to
            the raw values predicted from the trees of the ensemble . The
            classes corresponds to that in the attribute :term:`classes_`.
            Regression and binary classification are special cases with
            ``k == 1``, otherwise ``k==n_classes``.
        """
        yield from self._staged_raw_predict(X)

    def predict(self, X):
        """Predict class for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : array, shape (n_samples,)
            The predicted values.
        """
        raw_predictions = self.decision_function(X)
        encoded_labels = \
            self.loss_._raw_prediction_to_decision(raw_predictions)
        return self.classes_.take(encoded_labels, axis=0)

    def staged_predict(self, X):
        """Predict class at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : generator of array of shape (n_samples,)
            The predicted value of the input samples.
        """
        for raw_predictions in self._staged_raw_predict(X):
            encoded_labels = \
                self.loss_._raw_prediction_to_decision(raw_predictions)
            yield self.classes_.take(encoded_labels, axis=0)

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array, shape (n_samples, n_classes)
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        raw_predictions = self.decision_function(X)
        try:
            return self.loss_._raw_prediction_to_proba(raw_predictions)
        except NotFittedError:
            raise
        except AttributeError:
            raise AttributeError('loss=%r does not support predict_proba' %
                                 self.loss)

    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array, shape (n_samples, n_classes)
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        proba = self.predict_proba(X)
        return np.log(proba)

    def staged_predict_proba(self, X):
        """Predict class probabilities at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : generator of array of shape (n_samples,)
            The predicted value of the input samples.
        """
        try:
            for raw_predictions in self._staged_raw_predict(X):
                yield self.loss_._raw_prediction_to_proba(raw_predictions)
        except NotFittedError:
            raise
        except AttributeError:
            raise AttributeError('loss=%r does not support predict_proba' %
                                 self.loss)


class GradientBoostingRegressor(RegressorMixin, BaseGradientBoosting):
    """Gradient Boosting for regression.

    GB builds an additive model in a forward stage-wise fashion;
    it allows for the optimization of arbitrary differentiable loss functions.
    In each stage a regression tree is fit on the negative gradient of the
    given loss function.

    Read more in the :ref:`User Guide <gradient_boosting>`.

    Parameters
    ----------
    loss : {'ls', 'lad', 'huber', 'quantile'}, optional (default='ls')
        loss function to be optimized. 'ls' refers to least squares
        regression. 'lad' (least absolute deviation) is a highly robust
        loss function solely based on order information of the input
        variables. 'huber' is a combination of the two. 'quantile'
        allows quantile regression (use `alpha` to specify the quantile).

    learning_rate : float, optional (default=0.1)
        learning rate shrinks the contribution of each tree by `learning_rate`.
        There is a trade-off between learning_rate and n_estimators.

    n_estimators : int (default=100)
        The number of boosting stages to perform. Gradient boosting
        is fairly robust to over-fitting so a large number usually
        results in better performance.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_estimators`.
        Choosing `subsample < 1.0` leads to a reduction of variance
        and an increase in bias.

    criterion : string, optional (default="friedman_mse")
        The function to measure the quality of a split. Supported criteria
        are "friedman_mse" for the mean squared error with improvement
        score by Friedman, "mse" for mean squared error, and "mae" for
        the mean absolute error. The default value of "friedman_mse" is
        generally the best as it can provide a better approximation in
        some cases.

        .. versionadded:: 0.18

    min_samples_split : int, float, optional (default=2)
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

    min_impurity_decrease : float, optional (default=0.)
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    min_impurity_split : float, (default=1e-7)
        Threshold for early stopping in tree growth. A node will split
        if its impurity is above the threshold, otherwise it is a leaf.

        .. deprecated:: 0.19
           ``min_impurity_split`` has been deprecated in favor of
           ``min_impurity_decrease`` in 0.19. The default value of
           ``min_impurity_split`` will change from 1e-7 to 0 in 0.23 and it
           will be removed in 0.25. Use ``min_impurity_decrease`` instead.

    init : estimator or 'zero', optional (default=None)
        An estimator object that is used to compute the initial predictions.
        ``init`` has to provide :term:`fit` and :term:`predict`. If 'zero', the
        initial raw predictions are set to zero. By default a
        ``DummyEstimator`` is used, predicting either the average target value
        (for loss='ls'), or a quantile for the other losses.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `int(max_features * n_features)` features are considered at each
          split.
        - If "auto", then `max_features=n_features`.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.

        Choosing `max_features < n_features` leads to a reduction of variance
        and an increase in bias.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    alpha : float (default=0.9)
        The alpha-quantile of the huber loss function and the quantile
        loss function. Only if ``loss='huber'`` or ``loss='quantile'``.

    verbose : int, default: 0
        Enable verbose output. If 1 then it prints progress and performance
        once in a while (the more trees the lower the frequency). If greater
        than 1 then it prints progress and performance for every tree.

    max_leaf_nodes : int or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    warm_start : bool, default: False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just erase the
        previous solution. See :term:`the Glossary <warm_start>`.

    presort : bool or 'auto', optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error.

        .. versionadded:: 0.17
           optional parameter *presort*.

    validation_fraction : float, optional, default 0.1
        The proportion of training data to set aside as validation set for
        early stopping. Must be between 0 and 1.
        Only used if ``n_iter_no_change`` is set to an integer.

        .. versionadded:: 0.20

    n_iter_no_change : int, default None
        ``n_iter_no_change`` is used to decide if early stopping will be used
        to terminate training when validation score is not improving. By
        default it is set to None to disable early stopping. If set to a
        number, it will set aside ``validation_fraction`` size of the training
        data as validation and terminate training when validation score is not
        improving in all of the previous ``n_iter_no_change`` numbers of
        iterations.

        .. versionadded:: 0.20

    tol : float, optional, default 1e-4
        Tolerance for the early stopping. When the loss is not improving
        by at least tol for ``n_iter_no_change`` iterations (if set to a
        number), the training stops.

        .. versionadded:: 0.20

    ccp_alpha : non-negative float, optional (default=0.0)
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    Attributes
    ----------
    feature_importances_ : array, shape (n_features,)
        The feature importances (the higher, the more important the feature).

    oob_improvement_ : array, shape (n_estimators,)
        The improvement in loss (= deviance) on the out-of-bag samples
        relative to the previous iteration.
        ``oob_improvement_[0]`` is the improvement in
        loss of the first stage over the ``init`` estimator.

    train_score_ : array, shape (n_estimators,)
        The i-th score ``train_score_[i]`` is the deviance (= loss) of the
        model at iteration ``i`` on the in-bag sample.
        If ``subsample == 1`` this is the deviance on the training data.

    loss_ : LossFunction
        The concrete ``LossFunction`` object.

    init_ : estimator
        The estimator that provides the initial predictions.
        Set via the ``init`` argument or ``loss.init_estimator``.

    estimators_ : array of DecisionTreeRegressor, shape (n_estimators, 1)
        The collection of fitted sub-estimators.

    Notes
    -----
    The features are always randomly permuted at each split. Therefore,
    the best found split may vary, even with the same training data and
    ``max_features=n_features``, if the improvement of the criterion is
    identical for several splits enumerated during the search of the best
    split. To obtain a deterministic behaviour during fitting,
    ``random_state`` has to be fixed.

    See also
    --------
    sklearn.ensemble.HistGradientBoostingRegressor,
    sklearn.tree.DecisionTreeRegressor, RandomForestRegressor

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.

    J. Friedman, Stochastic Gradient Boosting, 1999

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning Ed. 2, Springer, 2009.
    """

    _SUPPORTED_LOSS = ('ls', 'lad', 'huber', 'quantile')

    def __init__(self, loss='ls', learning_rate=0.1, n_estimators=100,
                 subsample=1.0, criterion='friedman_mse', min_samples_split=2,
                 min_samples_leaf=1, min_weight_fraction_leaf=0.,
                 max_depth=3, min_impurity_decrease=0.,
                 min_impurity_split=None, init=None, random_state=None,
                 max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None,
                 warm_start=False, presort='auto', validation_fraction=0.1,
                 n_iter_no_change=None, tol=1e-4, ccp_alpha=0.0):

        super().__init__(
            loss=loss, learning_rate=learning_rate, n_estimators=n_estimators,
            criterion=criterion, min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_depth=max_depth, init=init, subsample=subsample,
            max_features=max_features,
            min_impurity_decrease=min_impurity_decrease,
            min_impurity_split=min_impurity_split,
            random_state=random_state, alpha=alpha, verbose=verbose,
            max_leaf_nodes=max_leaf_nodes, warm_start=warm_start,
            presort=presort, validation_fraction=validation_fraction,
            n_iter_no_change=n_iter_no_change, tol=tol, ccp_alpha=ccp_alpha)

    def predict(self, X):
        """Predict regression target for X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : array, shape (n_samples,)
            The predicted values.
        """
        X = check_array(X, dtype=DTYPE, order="C", accept_sparse='csr')
        # In regression we can directly return the raw value from the trees.
        return self._raw_predict(X).ravel()

    def staged_predict(self, X):
        """Predict regression target at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : generator of array of shape (n_samples,)
            The predicted value of the input samples.
        """
        for raw_predictions in self._staged_raw_predict(X):
            yield raw_predictions.ravel()

    def apply(self, X):
        """Apply trees in the ensemble to X, return leaf indices.

        .. versionadded:: 0.17

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will
            be converted to a sparse ``csr_matrix``.

        Returns
        -------
        X_leaves : array-like, shape (n_samples, n_estimators)
            For each datapoint x in X and for each tree in the ensemble,
            return the index of the leaf x ends up in each estimator.
        """

        leaves = super().apply(X)
        leaves = leaves.reshape(X.shape[0], self.estimators_.shape[0])
        return leaves
