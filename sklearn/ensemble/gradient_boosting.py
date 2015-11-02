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

from __future__ import print_function
from __future__ import division

from abc import ABCMeta
from abc import abstractmethod

from .base import BaseEnsemble
from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..base import RegressorMixin
from ..externals import six
from ..feature_selection.from_model import _LearntSelectorMixin

from ._gradient_boosting import predict_stages
from ._gradient_boosting import predict_stage
from ._gradient_boosting import _random_sample_mask

import numbers
import numpy as np
import warnings

from scipy import stats
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

from time import time
from ..tree.tree import DecisionTreeRegressor
from ..tree._tree import DTYPE
from ..tree._tree import TREE_LEAF

from ..utils import check_random_state
from ..utils import check_array
from ..utils import check_X_y
from ..utils import column_or_1d
from ..utils import check_consistent_length
from ..utils import deprecated
from ..utils.extmath import logsumexp
from ..utils.fixes import expit
from ..utils.fixes import bincount
from ..utils.stats import _weighted_percentile
from ..utils.validation import check_is_fitted
from ..utils.multiclass import check_classification_targets
from ..exceptions import NotFittedError

from ..metrics.scorer import check_scoring
from ..model_selection import check_cv
from ..model_selection import ParameterGrid
from ..externals.joblib import Parallel, delayed
from ..model_selection._search import _CVScoreTuple


class QuantileEstimator(BaseEstimator):
    """An estimator predicting the alpha-quantile of the training targets."""
    def __init__(self, alpha=0.9):
        if not 0 < alpha < 1.0:
            raise ValueError("`alpha` must be in (0, 1.0) but was %r" % alpha)
        self.alpha = alpha

    def fit(self, X, y, sample_weight=None):
        if sample_weight is None:
            self.quantile = stats.scoreatpercentile(y, self.alpha * 100.0)
        else:
            self.quantile = _weighted_percentile(y, sample_weight,
                                                 self.alpha * 100.0)

    def predict(self, X):
        check_is_fitted(self, 'quantile')

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.quantile)
        return y


class MeanEstimator(BaseEstimator):
    """An estimator predicting the mean of the training targets."""
    def fit(self, X, y, sample_weight=None):
        if sample_weight is None:
            self.mean = np.mean(y)
        else:
            self.mean = np.average(y, weights=sample_weight)

    def predict(self, X):
        check_is_fitted(self, 'mean')

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.mean)
        return y


class LogOddsEstimator(BaseEstimator):
    """An estimator predicting the log odds ratio."""
    scale = 1.0

    def fit(self, X, y, sample_weight=None):
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
        check_is_fitted(self, 'prior')

        y = np.empty((X.shape[0], 1), dtype=np.float64)
        y.fill(self.prior)
        return y


class ScaledLogOddsEstimator(LogOddsEstimator):
    """Log odds ratio scaled by 0.5 -- for exponential loss. """
    scale = 0.5


class PriorProbabilityEstimator(BaseEstimator):
    """An estimator predicting the probability of each
    class in the training data.
    """
    def fit(self, X, y, sample_weight=None):
        if sample_weight is None:
            sample_weight = np.ones_like(y, dtype=np.float64)
        class_counts = bincount(y, weights=sample_weight)
        self.priors = class_counts / class_counts.sum()

    def predict(self, X):
        check_is_fitted(self, 'priors')

        y = np.empty((X.shape[0], self.priors.shape[0]), dtype=np.float64)
        y[:] = self.priors
        return y


class ZeroEstimator(BaseEstimator):
    """An estimator that simply predicts zero. """

    def fit(self, X, y, sample_weight=None):
        if np.issubdtype(y.dtype, int):
            # classification
            self.n_classes = np.unique(y).shape[0]
            if self.n_classes == 2:
                self.n_classes = 1
        else:
            # regression
            self.n_classes = 1

    def predict(self, X):
        check_is_fitted(self, 'n_classes')

        y = np.empty((X.shape[0], self.n_classes), dtype=np.float64)
        y.fill(0.0)
        return y


class LossFunction(six.with_metaclass(ABCMeta, object)):
    """Abstract base class for various loss functions.

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
        """Compute the loss of prediction ``pred`` and ``y``. """

    @abstractmethod
    def negative_gradient(self, y, y_pred, **kargs):
        """Compute the negative gradient.

        Parameters
        ---------
        y : np.ndarray, shape=(n,)
            The target labels.
        y_pred : np.ndarray, shape=(n,):
            The predictions.
        """

    def update_terminal_regions(self, tree, X, y, residual, y_pred,
                                sample_weight, sample_mask,
                                learning_rate=1.0, k=0):
        """Update the terminal regions (=leaves) of the given tree and
        updates the current predictions of the model. Traverses tree
        and invokes template method `_update_terminal_region`.

        Parameters
        ----------
        tree : tree.Tree
            The tree object.
        X : ndarray, shape=(n, m)
            The data array.
        y : ndarray, shape=(n,)
            The target labels.
        residual : ndarray, shape=(n,)
            The residuals (usually the negative gradient).
        y_pred : ndarray, shape=(n,)
            The predictions.
        sample_weight : ndarray, shape=(n,)
            The weight of each sample.
        sample_mask : ndarray, shape=(n,)
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


class RegressionLossFunction(six.with_metaclass(ABCMeta, LossFunction)):
    """Base class for regression loss functions. """

    def __init__(self, n_classes):
        if n_classes != 1:
            raise ValueError("``n_classes`` must be 1 for regression but "
                             "was %r" % n_classes)
        super(RegressionLossFunction, self).__init__(n_classes)


class LeastSquaresError(RegressionLossFunction):
    """Loss function for least squares (LS) estimation.
    Terminal regions need not to be updated for least squares. """
    def init_estimator(self):
        return MeanEstimator()

    def __call__(self, y, pred, sample_weight=None):
        if sample_weight is None:
            return np.mean((y - pred.ravel()) ** 2.0)
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * ((y - pred.ravel()) ** 2.0)))

    def negative_gradient(self, y, pred, **kargs):
        return y - pred.ravel()

    def update_terminal_regions(self, tree, X, y, residual, y_pred,
                                sample_weight, sample_mask,
                                learning_rate=1.0, k=0):
        """Least squares does not need to update terminal regions.

        But it has to update the predictions.
        """
        # update predictions
        y_pred[:, k] += learning_rate * tree.predict(X).ravel()

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        pass


class LeastAbsoluteError(RegressionLossFunction):
    """Loss function for least absolute deviation (LAD) regression. """
    def init_estimator(self):
        return QuantileEstimator(alpha=0.5)

    def __call__(self, y, pred, sample_weight=None):
        if sample_weight is None:
            return np.abs(y - pred.ravel()).mean()
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * np.abs(y - pred.ravel())))

    def negative_gradient(self, y, pred, **kargs):
        """1.0 if y - pred > 0.0 else -1.0"""
        pred = pred.ravel()
        return 2.0 * (y - pred > 0.0) - 1.0

    def _update_terminal_region(self, tree, terminal_regions, leaf, X, y,
                                residual, pred, sample_weight):
        """LAD updates terminal regions to median estimates. """
        terminal_region = np.where(terminal_regions == leaf)[0]
        sample_weight = sample_weight.take(terminal_region, axis=0)
        diff = y.take(terminal_region, axis=0) - pred.take(terminal_region, axis=0)
        tree.value[leaf, 0, 0] = _weighted_percentile(diff, sample_weight, percentile=50)


class HuberLossFunction(RegressionLossFunction):
    """Huber loss function for robust regression.

    M-Regression proposed in Friedman 2001.

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.
    """

    def __init__(self, n_classes, alpha=0.9):
        super(HuberLossFunction, self).__init__(n_classes)
        self.alpha = alpha
        self.gamma = None

    def init_estimator(self):
        return QuantileEstimator(alpha=0.5)

    def __call__(self, y, pred, sample_weight=None):
        pred = pred.ravel()
        diff = y - pred
        gamma = self.gamma
        if gamma is None:
            if sample_weight is None:
                gamma = stats.scoreatpercentile(np.abs(diff), self.alpha * 100)
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
        pred = pred.ravel()
        diff = y - pred
        if sample_weight is None:
            gamma = stats.scoreatpercentile(np.abs(diff), self.alpha * 100)
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


class QuantileLossFunction(RegressionLossFunction):
    """Loss function for quantile regression.

    Quantile regression allows to estimate the percentiles
    of the conditional distribution of the target.
    """

    def __init__(self, n_classes, alpha=0.9):
        super(QuantileLossFunction, self).__init__(n_classes)
        assert 0 < alpha < 1.0
        self.alpha = alpha
        self.percentile = alpha * 100.0

    def init_estimator(self):
        return QuantileEstimator(self.alpha)

    def __call__(self, y, pred, sample_weight=None):
        pred = pred.ravel()
        diff = y - pred
        alpha = self.alpha

        mask = y > pred
        if sample_weight is None:
            loss = (alpha * diff[mask].sum() +
                    (1.0 - alpha) * diff[~mask].sum()) / y.shape[0]
        else:
            loss = ((alpha * np.sum(sample_weight[mask] * diff[mask]) +
                    (1.0 - alpha) * np.sum(sample_weight[~mask] * diff[~mask])) /
                    sample_weight.sum())
        return loss

    def negative_gradient(self, y, pred, **kargs):
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


class ClassificationLossFunction(six.with_metaclass(ABCMeta, LossFunction)):
    """Base class for classification loss functions. """

    def _score_to_proba(self, score):
        """Template method to convert scores to probabilities.

         the does not support probabilites raises AttributeError.
        """
        raise TypeError('%s does not support predict_proba' % type(self).__name__)

    @abstractmethod
    def _score_to_decision(self, score):
        """Template method to convert scores to decisions.

        Returns int arrays.
        """


class BinomialDeviance(ClassificationLossFunction):
    """Binomial deviance loss function for binary classification.

    Binary classification is a special case; here, we only need to
    fit one tree instead of ``n_classes`` trees.
    """
    def __init__(self, n_classes):
        if n_classes != 2:
            raise ValueError("{0:s} requires 2 classes.".format(
                self.__class__.__name__))
        # we only need to fit one tree for binary clf.
        super(BinomialDeviance, self).__init__(1)

    def init_estimator(self):
        return LogOddsEstimator()

    def __call__(self, y, pred, sample_weight=None):
        """Compute the deviance (= 2 * negative log-likelihood). """
        # logaddexp(0, v) == log(1.0 + exp(v))
        pred = pred.ravel()
        if sample_weight is None:
            return -2.0 * np.mean((y * pred) - np.logaddexp(0.0, pred))
        else:
            return (-2.0 / sample_weight.sum() *
                    np.sum(sample_weight * ((y * pred) - np.logaddexp(0.0, pred))))

    def negative_gradient(self, y, pred, **kargs):
        """Compute the residual (= negative gradient). """
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

        if denominator == 0.0:
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


class MultinomialDeviance(ClassificationLossFunction):
    """Multinomial deviance loss function for multi-class classification.

    For multi-class classification we need to fit ``n_classes`` trees at
    each stage.
    """

    is_multi_class = True

    def __init__(self, n_classes):
        if n_classes < 3:
            raise ValueError("{0:s} requires more than 2 classes.".format(
                self.__class__.__name__))
        super(MultinomialDeviance, self).__init__(n_classes)

    def init_estimator(self):
        return PriorProbabilityEstimator()

    def __call__(self, y, pred, sample_weight=None):
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
        """Compute negative gradient for the ``k``-th class. """
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

        if denominator == 0.0:
            tree.value[leaf, 0, 0] = 0.0
        else:
            tree.value[leaf, 0, 0] = numerator / denominator

    def _score_to_proba(self, score):
        return np.nan_to_num(
            np.exp(score - (logsumexp(score, axis=1)[:, np.newaxis])))

    def _score_to_decision(self, score):
        proba = self._score_to_proba(score)
        return np.argmax(proba, axis=1)


class ExponentialLoss(ClassificationLossFunction):
    """Exponential loss function for binary classification.

    Same loss as AdaBoost.

    References
    ----------
    Greg Ridgeway, Generalized Boosted Models: A guide to the gbm package, 2007
    """
    def __init__(self, n_classes):
        if n_classes != 2:
            raise ValueError("{0:s} requires 2 classes.".format(
                self.__class__.__name__))
        # we only need to fit one tree for binary clf.
        super(ExponentialLoss, self).__init__(1)

    def init_estimator(self):
        return ScaledLogOddsEstimator()

    def __call__(self, y, pred, sample_weight=None):
        pred = pred.ravel()
        if sample_weight is None:
            return np.mean(np.exp(-(2. * y - 1.) * pred))
        else:
            return (1.0 / sample_weight.sum() *
                    np.sum(sample_weight * np.exp(-(2 * y - 1) * pred)))

    def negative_gradient(self, y, pred, **kargs):
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

        if denominator == 0.0:
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


LOSS_FUNCTIONS = {'ls': LeastSquaresError,
                  'lad': LeastAbsoluteError,
                  'huber': HuberLossFunction,
                  'quantile': QuantileLossFunction,
                  'deviance': None,    # for both, multinomial and binomial
                  'exponential': ExponentialLoss,
                  }


INIT_ESTIMATORS = {'zero': ZeroEstimator}


class VerboseReporter(object):
    """Reports verbose output to stdout.

    If ``verbose==1`` output is printed once in a while (when iteration mod
    verbose_mod is zero).; if larger than 1 then output is printed for
    each update.
    """

    def __init__(self, verbose):
        self.verbose = verbose

    def init(self, est, begin_at_stage=0):
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
        """Update reporter with new iteration. """
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


class BaseGradientBoosting(six.with_metaclass(ABCMeta, BaseEnsemble,
                                              _LearntSelectorMixin)):
    """Abstract base class for Gradient Boosting. """

    @abstractmethod
    def __init__(self, loss, learning_rate, n_estimators, criterion,
                 min_samples_split, min_samples_leaf, min_weight_fraction_leaf,
                 max_depth, min_impurity_split, init, subsample, max_features,
                 random_state, alpha=0.9, verbose=0, max_leaf_nodes=None,
                 warm_start=False, presort='auto'):

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
        self.min_impurity_split = min_impurity_split
        self.init = init
        self.random_state = random_state
        self.alpha = alpha
        self.verbose = verbose
        self.max_leaf_nodes = max_leaf_nodes
        self.warm_start = warm_start
        self.presort = presort

        self.estimators_ = np.empty((0, 0), dtype=np.object)

    def _fit_stage(self, i, X, y, y_pred, sample_weight, sample_mask,
                   random_state, X_idx_sorted, X_csc=None, X_csr=None):
        """Fit another stage of ``n_classes_`` trees to the boosting model. """

        assert sample_mask.dtype == np.bool
        loss = self.loss_
        original_y = y

        for k in range(loss.K):
            if loss.is_multi_class:
                y = np.array(original_y == k, dtype=np.float64)

            residual = loss.negative_gradient(y, y_pred, k=k,
                                              sample_weight=sample_weight)

            # induce regression tree on residuals
            tree = DecisionTreeRegressor(
                criterion=self.criterion,
                splitter='best',
                max_depth=self.max_depth,
                min_samples_split=self.min_samples_split,
                min_samples_leaf=self.min_samples_leaf,
                min_weight_fraction_leaf=self.min_weight_fraction_leaf,
                max_features=self.max_features,
                max_leaf_nodes=self.max_leaf_nodes,
                random_state=random_state,
                presort=self.presort)

            if self.subsample < 1.0:
                # no inplace multiplication!
                sample_weight = sample_weight * sample_mask.astype(np.float64)

            if X_csc is not None:
                tree.fit(X_csc, residual, sample_weight=sample_weight,
                         check_input=False, X_idx_sorted=X_idx_sorted)
            else:
                tree.fit(X, residual, sample_weight=sample_weight,
                         check_input=False, X_idx_sorted=X_idx_sorted)

            # update tree leaves
            if X_csr is not None:
                loss.update_terminal_regions(tree.tree_, X_csr, y, residual, y_pred,
                                             sample_weight, sample_mask,
                                             self.learning_rate, k=k)
            else:
                loss.update_terminal_regions(tree.tree_, X, y, residual, y_pred,
                                             sample_weight, sample_mask,
                                             self.learning_rate, k=k)

            # add tree to ensemble
            self.estimators_[i, k] = tree

        return y_pred

    def _check_params(self):
        """Check validity of parameters and raise ValueError if not valid. """
        if self.n_estimators <= 0:
            raise ValueError("n_estimators must be greater than 0 but "
                             "was %r" % self.n_estimators)

        if self.learning_rate <= 0.0:
            raise ValueError("learning_rate must be greater than 0 but "
                             "was %r" % self.learning_rate)

        if (self.loss not in self._SUPPORTED_LOSS
                or self.loss not in LOSS_FUNCTIONS):
            raise ValueError("Loss '{0:s}' not supported. ".format(self.loss))

        if self.loss == 'deviance':
            loss_class = (MultinomialDeviance
                          if len(self.classes_) > 2
                          else BinomialDeviance)
        else:
            loss_class = LOSS_FUNCTIONS[self.loss]

        if self.loss in ('huber', 'quantile'):
            self.loss_ = loss_class(self.n_classes_, self.alpha)
        else:
            self.loss_ = loss_class(self.n_classes_)

        if not (0.0 < self.subsample <= 1.0):
            raise ValueError("subsample must be in (0,1] but "
                             "was %r" % self.subsample)

        if self.init is not None:
            if isinstance(self.init, six.string_types):
                if self.init not in INIT_ESTIMATORS:
                    raise ValueError('init="%s" is not supported' % self.init)
            else:
                if (not hasattr(self.init, 'fit')
                        or not hasattr(self.init, 'predict')):
                    raise ValueError("init=%r must be valid BaseEstimator "
                                     "and support both fit and "
                                     "predict" % self.init)

        if not (0.0 < self.alpha < 1.0):
            raise ValueError("alpha must be in (0.0, 1.0) but "
                             "was %r" % self.alpha)

        if isinstance(self.max_features, six.string_types):
            if self.max_features == "auto":
                # if is_classification
                if self.n_classes_ > 1:
                    max_features = max(1, int(np.sqrt(self.n_features)))
                else:
                    # is regression
                    max_features = self.n_features
            elif self.max_features == "sqrt":
                max_features = max(1, int(np.sqrt(self.n_features)))
            elif self.max_features == "log2":
                max_features = max(1, int(np.log2(self.n_features)))
            else:
                raise ValueError("Invalid value for max_features: %r. "
                                 "Allowed string values are 'auto', 'sqrt' "
                                 "or 'log2'." % self.max_features)
        elif self.max_features is None:
            max_features = self.n_features
        elif isinstance(self.max_features, (numbers.Integral, np.integer)):
            max_features = self.max_features
        else:  # float
            if 0. < self.max_features <= 1.:
                max_features = max(int(self.max_features * self.n_features), 1)
            else:
                raise ValueError("max_features must be in (0, n_features]")

        self.max_features_ = max_features

    def _init_state(self):
        """Initialize model state and allocate model state data structures. """

        if self.init is None:
            self.init_ = self.loss_.init_estimator()
        elif isinstance(self.init, six.string_types):
            self.init_ = INIT_ESTIMATORS[self.init]()
        else:
            self.init_ = self.init

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

    def _resize_state(self):
        """Add additional ``n_estimators`` entries to all attributes. """
        # self.n_estimators is the number of additional est to fit
        total_n_estimators = self.n_estimators
        if total_n_estimators < self.estimators_.shape[0]:
            raise ValueError('resize with smaller n_estimators %d < %d' %
                             (total_n_estimators, self.estimators_[0]))

        self.estimators_.resize((total_n_estimators, self.loss_.K))
        self.train_score_.resize(total_n_estimators)
        if (self.subsample < 1 or hasattr(self, 'oob_improvement_')):
            # if do oob resize arrays or create new if not available
            if hasattr(self, 'oob_improvement_'):
                self.oob_improvement_.resize(total_n_estimators)
            else:
                self.oob_improvement_ = np.zeros((total_n_estimators,),
                                                 dtype=np.float64)

    def _is_initialized(self):
        return len(getattr(self, 'estimators_', [])) > 0

    def _check_initialized(self):
        """Check that the estimator is initialized, raising an error if not."""
        if self.estimators_ is None or len(self.estimators_) == 0:
            raise NotFittedError("Estimator not fitted, call `fit`"
                                 " before making predictions`.")

    def fit(self, X, y, sample_weight=None, monitor=None):
        """Fit the gradient boosting model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes.

        sample_weight : array-like, shape = [n_samples] or None
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
            Returns self.
        """
        # if not warmstart - clear the estimator state
        if not self.warm_start:
            self._clear_state()

        # Check input
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'], dtype=DTYPE)
        n_samples, self.n_features = X.shape
        if sample_weight is None:
            sample_weight = np.ones(n_samples, dtype=np.float32)
        else:
            sample_weight = column_or_1d(sample_weight, warn=True)

        check_consistent_length(X, y, sample_weight)

        y = self._validate_y(y)

        random_state = check_random_state(self.random_state)
        self._check_params()

        if not self._is_initialized():
            # init state
            self._init_state()

            # fit initial model - FIXME make sample_weight optional
            self.init_.fit(X, y, sample_weight)

            # init predictions
            y_pred = self.init_.predict(X)
            begin_at_stage = 0
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
            y_pred = self._decision_function(X)
            self._resize_state()

        X_idx_sorted = None
        presort = self.presort
        # Allow presort to be 'auto', which means True if the dataset is dense,
        # otherwise it will be False.
        if presort == 'auto' and issparse(X):
            presort = False
        elif presort == 'auto':
            presort = True

        if presort == True:
            if issparse(X):
                raise ValueError("Presorting is not supported for sparse matrices.")
            else:
                X_idx_sorted = np.asfortranarray(np.argsort(X, axis=0),
                                                 dtype=np.int32)

        # fit the boosting stages
        n_stages = self._fit_stages(X, y, y_pred, sample_weight, random_state,
                                    begin_at_stage, monitor, X_idx_sorted)
        # change shape of arrays after fit (early-stopping or additional ests)
        if n_stages != self.estimators_.shape[0]:
            self.estimators_ = self.estimators_[:n_stages]
            self.train_score_ = self.train_score_[:n_stages]
            if hasattr(self, 'oob_improvement_'):
                self.oob_improvement_ = self.oob_improvement_[:n_stages]

        return self

    def _fit_stages(self, X, y, y_pred, sample_weight, random_state,
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

        # Set min_weight_leaf from min_weight_fraction_leaf
        if self.min_weight_fraction_leaf != 0. and sample_weight is not None:
            min_weight_leaf = (self.min_weight_fraction_leaf *
                               np.sum(sample_weight))
        else:
            min_weight_leaf = 0.

        if self.verbose:
            verbose_reporter = VerboseReporter(self.verbose)
            verbose_reporter.init(self, begin_at_stage)

        X_csc = csc_matrix(X) if issparse(X) else None
        X_csr = csr_matrix(X) if issparse(X) else None

        # perform boosting iterations
        i = begin_at_stage
        for i in range(begin_at_stage, self.n_estimators):

            # subsampling
            if do_oob:
                sample_mask = _random_sample_mask(n_samples, n_inbag,
                                                  random_state)
                # OOB score before adding this stage
                old_oob_score = loss_(y[~sample_mask],
                                      y_pred[~sample_mask],
                                      sample_weight[~sample_mask])

            # fit next stage of trees
            y_pred = self._fit_stage(i, X, y, y_pred, sample_weight,
                                     sample_mask, random_state, X_idx_sorted,
                                     X_csc, X_csr)

            # track deviance (= loss)
            if do_oob:
                self.train_score_[i] = loss_(y[sample_mask],
                                             y_pred[sample_mask],
                                             sample_weight[sample_mask])
                self.oob_improvement_[i] = (
                    old_oob_score - loss_(y[~sample_mask],
                                          y_pred[~sample_mask],
                                          sample_weight[~sample_mask]))
            else:
                # no need to fancy index w/ no subsampling
                self.train_score_[i] = loss_(y, y_pred, sample_weight)

            if self.verbose > 0:
                verbose_reporter.update(i, self)

            if monitor is not None:
                early_stopping = monitor(i, self, locals())
                if early_stopping:
                    break
        return i + 1

    def _make_estimator(self, append=True):
        # we don't need _make_estimator
        raise NotImplementedError()

    def _init_decision_function(self, X):
        """Check input and compute prediction of ``init``. """
        self._check_initialized()
        X = self.estimators_[0, 0]._validate_X_predict(X, check_input=True)
        if X.shape[1] != self.n_features:
            raise ValueError("X.shape[1] should be {0:d}, not {1:d}.".format(
                self.n_features, X.shape[1]))
        score = self.init_.predict(X).astype(np.float64)
        return score

    def _decision_function(self, X):
        # for use in inner loop, not raveling the output in single-class case,
        # not doing input validation.
        score = self._init_decision_function(X)
        predict_stages(self.estimators_, X, self.learning_rate, score)
        return score

    @deprecated(" and will be removed in 0.19")
    def decision_function(self, X):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : array, shape = [n_samples, n_classes] or [n_samples]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification produce an array of shape
            [n_samples].
        """

        self._check_initialized()
        X = self.estimators_[0, 0]._validate_X_predict(X, check_input=True)
        score = self._decision_function(X)
        if score.shape[1] == 1:
            return score.ravel()
        return score

    def _staged_decision_function(self, X):
        """Compute decision function of ``X`` for each iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : generator of array, shape = [n_samples, k]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification are special cases with
            ``k == 1``, otherwise ``k==n_classes``.
        """
        X = check_array(X, dtype=DTYPE, order="C")
        score = self._init_decision_function(X)
        for i in range(self.estimators_.shape[0]):
            predict_stage(self.estimators_, i, X, self.learning_rate, score)
            yield score.copy()

    @deprecated(" and will be removed in 0.19")
    def staged_decision_function(self, X):
        """Compute decision function of ``X`` for each iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : generator of array, shape = [n_samples, k]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification are special cases with
            ``k == 1``, otherwise ``k==n_classes``.
        """
        for dec in self._staged_decision_function(X):
            # no yield from in Python2.X
            yield dec

    @property
    def feature_importances_(self):
        """Return the feature importances (the higher, the more important the
           feature).

        Returns
        -------
        feature_importances_ : array, shape = [n_features]
        """
        self._check_initialized()

        total_sum = np.zeros((self.n_features, ), dtype=np.float64)
        for stage in self.estimators_:
            stage_sum = sum(tree.feature_importances_
                            for tree in stage) / len(stage)
            total_sum += stage_sum

        importances = total_sum / len(self.estimators_)
        return importances

    def _validate_y(self, y):
        self.n_classes_ = 1
        if y.dtype.kind == 'O':
            y = y.astype(np.float64)
        # Default implementation
        return y

    def apply(self, X):
        """Apply trees in the ensemble to X, return leaf indices.

        .. versionadded:: 0.17

        Parameters
        ----------
        X : array-like or sparse matrix, shape = [n_samples, n_features]
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will
            be converted to a sparse ``csr_matrix``.

        Returns
        -------
        X_leaves : array_like, shape = [n_samples, n_estimators, n_classes]
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


class GradientBoostingClassifier(BaseGradientBoosting, ClassifierMixin):
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

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

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
        - If float, then `min_samples_split` is a percentage and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node:

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a percentage and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.


    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the input samples required to be at a
        leaf node.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_estimators`.
        Choosing `subsample < 1.0` leads to a reduction of variance
        and an increase in bias.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a percentage and
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

    max_leaf_nodes : int or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_split : float, optional (default=1e-7)
        Threshold for early stopping in tree growth. A node will split
        if its impurity is above the threshold, otherwise it is a leaf.

        .. versionadded:: 0.18

    init : BaseEstimator, None, optional (default=None)
        An estimator object that is used to compute the initial
        predictions. ``init`` has to provide ``fit`` and ``predict``.
        If None it uses ``loss.init_estimator``.

    verbose : int, default: 0
        Enable verbose output. If 1 then it prints progress and performance
        once in a while (the more trees the lower the frequency). If greater
        than 1 then it prints progress and performance for every tree.

    warm_start : bool, default: False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just erase the
        previous solution.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    presort : bool or 'auto', optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error.

        .. versionadded:: 0.17
           *presort* parameter.

    Attributes
    ----------
    feature_importances_ : array, shape = [n_features]
        The feature importances (the higher, the more important the feature).

    oob_improvement_ : array, shape = [n_estimators]
        The improvement in loss (= deviance) on the out-of-bag samples
        relative to the previous iteration.
        ``oob_improvement_[0]`` is the improvement in
        loss of the first stage over the ``init`` estimator.

    train_score_ : array, shape = [n_estimators]
        The i-th score ``train_score_[i]`` is the deviance (= loss) of the
        model at iteration ``i`` on the in-bag sample.
        If ``subsample == 1`` this is the deviance on the training data.

    loss_ : LossFunction
        The concrete ``LossFunction`` object.

    init : BaseEstimator
        The estimator that provides the initial predictions.
        Set via the ``init`` argument or ``loss.init_estimator``.

    estimators_ : ndarray of DecisionTreeRegressor, shape = [n_estimators, ``loss_.K``]
        The collection of fitted sub-estimators. ``loss_.K`` is 1 for binary
        classification, otherwise n_classes.


    See also
    --------
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
                 max_depth=3, min_impurity_split=1e-7, init=None,
                 random_state=None, max_features=None, verbose=0,
                 max_leaf_nodes=None, warm_start=False,
                 presort='auto'):

        super(GradientBoostingClassifier, self).__init__(
            loss=loss, learning_rate=learning_rate, n_estimators=n_estimators,
            criterion=criterion, min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_depth=max_depth, init=init, subsample=subsample,
            max_features=max_features,
            random_state=random_state, verbose=verbose,
            max_leaf_nodes=max_leaf_nodes,
            min_impurity_split=min_impurity_split,
            warm_start=warm_start,
            presort=presort)

    def _validate_y(self, y):
        check_classification_targets(y)
        self.classes_, y = np.unique(y, return_inverse=True)
        self.n_classes_ = len(self.classes_)
        return y

    def decision_function(self, X):
        """Compute the decision function of ``X``.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : array, shape = [n_samples, n_classes] or [n_samples]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification produce an array of shape
            [n_samples].
        """
        X = check_array(X, dtype=DTYPE, order="C")
        score = self._decision_function(X)
        if score.shape[1] == 1:
            return score.ravel()
        return score

    def staged_decision_function(self, X):
        """Compute decision function of ``X`` for each iteration.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        score : generator of array, shape = [n_samples, k]
            The decision function of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
            Regression and binary classification are special cases with
            ``k == 1``, otherwise ``k==n_classes``.
        """
        for dec in self._staged_decision_function(X):
            # no yield from in Python2.X
            yield dec

    def predict(self, X):
        """Predict class for X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y: array of shape = ["n_samples]
            The predicted values.
        """
        score = self.decision_function(X)
        decisions = self.loss_._score_to_decision(score)
        return self.classes_.take(decisions, axis=0)

    def staged_predict(self, X):
        """Predict class at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : generator of array of shape = [n_samples]
            The predicted value of the input samples.
        """
        for score in self._staged_decision_function(X):
            decisions = self.loss_._score_to_decision(score)
            yield self.classes_.take(decisions, axis=0)

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        score = self.decision_function(X)
        try:
            return self.loss_._score_to_proba(score)
        except NotFittedError:
            raise
        except AttributeError:
            raise AttributeError('loss=%r does not support predict_proba' %
                                 self.loss)

    def predict_log_proba(self, X):
        """Predict class log-probabilities for X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        proba = self.predict_proba(X)
        return np.log(proba)

    def staged_predict_proba(self, X):
        """Predict class probabilities at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : generator of array of shape = [n_samples]
            The predicted value of the input samples.
        """
        try:
            for score in self._staged_decision_function(X):
                yield self.loss_._score_to_proba(score)
        except NotFittedError:
            raise
        except AttributeError:
            raise AttributeError('loss=%r does not support predict_proba' %
                                 self.loss)


class GradientBoostingRegressor(BaseGradientBoosting, RegressorMixin):
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

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

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
        - If float, then `min_samples_split` is a percentage and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

    min_samples_leaf : int, float, optional (default=1)
        The minimum number of samples required to be at a leaf node:

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a percentage and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

    min_weight_fraction_leaf : float, optional (default=0.)
        The minimum weighted fraction of the input samples required to be at a
        leaf node.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_estimators`.
        Choosing `subsample < 1.0` leads to a reduction of variance
        and an increase in bias.

    max_features : int, float, string or None, optional (default=None)
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a percentage and
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

    max_leaf_nodes : int or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_split : float, optional (default=1e-7)
        Threshold for early stopping in tree growth. A node will split
        if its impurity is above the threshold, otherwise it is a leaf.

        .. versionadded:: 0.18

    alpha : float (default=0.9)
        The alpha-quantile of the huber loss function and the quantile
        loss function. Only if ``loss='huber'`` or ``loss='quantile'``.

    init : BaseEstimator, None, optional (default=None)
        An estimator object that is used to compute the initial
        predictions. ``init`` has to provide ``fit`` and ``predict``.
        If None it uses ``loss.init_estimator``.

    verbose : int, default: 0
        Enable verbose output. If 1 then it prints progress and performance
        once in a while (the more trees the lower the frequency). If greater
        than 1 then it prints progress and performance for every tree.

    warm_start : bool, default: False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just erase the
        previous solution.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    presort : bool or 'auto', optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error.

        .. versionadded:: 0.17
           optional parameter *presort*.

    Attributes
    ----------
    feature_importances_ : array, shape = [n_features]
        The feature importances (the higher, the more important the feature).

    oob_improvement_ : array, shape = [n_estimators]
        The improvement in loss (= deviance) on the out-of-bag samples
        relative to the previous iteration.
        ``oob_improvement_[0]`` is the improvement in
        loss of the first stage over the ``init`` estimator.

    train_score_ : array, shape = [n_estimators]
        The i-th score ``train_score_[i]`` is the deviance (= loss) of the
        model at iteration ``i`` on the in-bag sample.
        If ``subsample == 1`` this is the deviance on the training data.

    loss_ : LossFunction
        The concrete ``LossFunction`` object.

    `init` : BaseEstimator
        The estimator that provides the initial predictions.
        Set via the ``init`` argument or ``loss.init_estimator``.

    estimators_ : ndarray of DecisionTreeRegressor, shape = [n_estimators, 1]
        The collection of fitted sub-estimators.

    See also
    --------
    DecisionTreeRegressor, RandomForestRegressor

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
                 max_depth=3, min_impurity_split=1e-7, init=None, random_state=None,
                 max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None,
                 warm_start=False, presort='auto'):

        super(GradientBoostingRegressor, self).__init__(
            loss=loss, learning_rate=learning_rate, n_estimators=n_estimators,
            criterion=criterion, min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_depth=max_depth, init=init, subsample=subsample,
            max_features=max_features, min_impurity_split=min_impurity_split,
            random_state=random_state, alpha=alpha, verbose=verbose,
            max_leaf_nodes=max_leaf_nodes, warm_start=warm_start,
            presort=presort)

    def predict(self, X):
        """Predict regression target for X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples]
            The predicted values.
        """
        X = check_array(X, dtype=DTYPE, order="C")
        return self._decision_function(X).ravel()

    def staged_predict(self, X):
        """Predict regression target at each stage for X.

        This method allows monitoring (i.e. determine error on testing set)
        after each stage.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : generator of array of shape = [n_samples]
            The predicted value of the input samples.
        """
        for y in self._staged_decision_function(X):
            yield y.ravel()

    def apply(self, X):
        """Apply trees in the ensemble to X, return leaf indices.

        .. versionadded:: 0.17

        Parameters
        ----------
        X : array-like or sparse matrix, shape = [n_samples, n_features]
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will
            be converted to a sparse ``csr_matrix``.

        Returns
        -------
        X_leaves : array_like, shape = [n_samples, n_estimators]
            For each datapoint x in X and for each tree in the ensemble,
            return the index of the leaf x ends up in each estimator.
        """

        leaves = super(GradientBoostingRegressor, self).apply(X)
        leaves = leaves.reshape(X.shape[0], self.estimators_.shape[0])
        return leaves


class GradientBoostingClassifierCV(BaseEstimator):
    """ Gradient Boosting Classification with Cross Validation

    This class keeps adding estimators to the underlying GB estimator till
    adding a new estimators does not improve the score by much. The goal of
    this class is to determine the best choice of `n_estimators` for each set
    of parameters. In most cases, this will be much faster than using
    :class:`sklearn.model_selection.GridSearchCV` due to early stopping.
    All parameters which are not parsed by this class are passed directly to
    the underlying Gradient Boosting estimator. If one or more parameters is a
    list, the Gradient Boosting model is evaluated for all elements of the
    list.

    Parameters
    ----------
    n_stop_rounds : int, optional, default=10
        If the score on the test set rounded off to `score_precision` decimal
        places does not change for `n_stop_rounds` iterations, the gradient
        boosting is halted.

    score_precision : int, optional, default=2
        The number of decimal places to round the score off by before
        comparing it with the scores of previous iterations.

    max_iterations : int, optinal, default=10000
        The maximum number of estimators that will be added to the Gradient
        Boosting class. If the model converges, the number of estimators will
        be lower than this.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a `(Stratified)KFold`,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, test splits.

        For integer/None inputs, ``StratifiedKFold`` is used for classification
        tasks, when ``y`` is binary or multiclass.

        See the :mod:`sklearn.model_selection` module for the list of
        cross-validation strategies that can be used here.

        Also refer :ref:`cross-validation documentation <cross_validation>`

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    refit : boolean, default=True
        If set, the data will be fit on a new Gradient Boosting model with
        the parameters that gave the best mean score. If not set, it is not
        possbile to make predections later with this class.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    loss : string, list of strings, optional (default='deviance')
        loss function to be optimized. 'deviance' refers to
        deviance (= logistic regression) for classification
        with probabilistic outputs. For loss 'exponential' gradient
        boosting recovers the AdaBoost algorithm. If a list is given, all
        loss functions are evaluated to choose the best one.

    learning_rate : float, list of floats, optional (default=0.1)
        learning rate shrinks the contribution of each tree by `learning_rate`.
        There is a trade-off between learning_rate and n_estimators. If a list
        is given, all learning rates will be evaluated.

    max_depth : integer, list of integers, optional (default=3)
        maximum depth of the individual regression estimators. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.
        Ignored if ``max_leaf_nodes`` is not None.
        If a list is given, each item in the list will be evaluated.

    min_samples_split : int, float, list of int or float, optional (default=2)
        The minimum number of samples required to split an internal node:
        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a percentage and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.
        If a list is given, each item in the list will be evaluated.

    min_samples_leaf : int, float, list of int or floats, optional (default=1)
        The minimum number of samples required to be at a leaf node:
        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a percentage and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.
        If a list is given, each item in the list will be evaluated.

    min_weight_fraction_leaf : float, list of floats, optional (default=0.)
        The minimum weighted fraction of the input samples required to be at a
        leaf node. If a list is given, each item in the list will be evaluated.

    subsample : float, list of floats, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_estimators`.
        Choosing `subsample < 1.0` leads to a reduction of variance
        and an increase in bias. If a list is given, each item in the list
        will be evaluated.

    max_features : int, float, string or None, list, optional (default=None)
        The number of features to consider when looking for the best split:
        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a percentage and
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
        effectively inspect more than ``max_features`` features. If a list is
        given, each item in the list will be evaluated.

    max_leaf_nodes : int, list of int, or None, optional (default=None)
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.
        If not None then ``max_depth`` will be ignored.
        If a list is given, each item in the list will be evaluated.

    init : BaseEstimator, None, list, optional (default=None)
        An estimator object that is used to compute the initial
        predictions. ``init`` has to provide ``fit`` and ``predict``.
        If None it uses ``loss.init_estimator``.
        If a list is given, each item in the list will be evaluated.

    random_state : int, RandomState instance or None,
                   list of int or RandomState, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    presort : bool or 'auto', list of string or bool optional (default='auto')
        Whether to presort the data to speed up the finding of best splits in
        fitting. Auto mode by default will use presorting on dense data and
        default to normal sorting on sparse data. Setting presort to true on
        sparse data will raise an error. If a list is given, each item in the
        list will be evaluated.

    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    best_score_tuple_ : named tuple
        The attributes of the model which gave the best score.

        * ``parameters``, a dict of parameter settings, which includes the
          `n_estimators` parameter
        * ``mean_validation_score``, the mean score over the
          cross-validation folds
        * ``cv_validation_scores``, the list of scores for each fold

    best_params_ : dict
        The parameters which gave the best score, along with the `n_estimators`
        parameter.

    best_estimator_ : GradientBoostingClassifier
        If `refit` was `True`, the model with the best score fit to the entire
        dataset.

    """

    _estimator_class = GradientBoostingClassifier

    def __init__(self, n_stop_rounds=10, score_precision=2,
                 max_iterations=10000, cv=3, scoring=None, refit=True,
                 n_jobs=1, pre_dispatch='2*n_jobs', verbose=0, loss='deviance',
                 learning_rate=0.1, subsample=1.0, min_samples_split=2,
                 min_samples_leaf=1, min_weight_fraction_leaf=0., max_depth=3,
                 init=None, random_state=None, max_features=None,
                 max_leaf_nodes=None, presort='auto'):

        self.n_stop_rounds = n_stop_rounds
        self.score_precision = score_precision
        self.max_iterations = max_iterations
        self.cv = cv
        self.scoring = scoring
        self.refit = refit
        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self.verbose = verbose

        self.loss = loss
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_depth = max_depth
        self.init = init
        self.random_state = random_state
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.presort = presort

    def fit(self, X, y):
        """Run fit with all sets of parameters till convergence.
    
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes.
        """
    
        params = self._get_params()
        X, y = check_X_y(X, y, dtype=DTYPE)
    
        cv = check_cv(self.cv, y)
        param_iter = ParameterGrid(params)
    
        parallel = Parallel(n_jobs=self.n_jobs, verbose=self.verbose,
                            pre_dispatch=self.pre_dispatch)
    
        out = parallel(delayed(_fit_single_param)
                       (self._estimator_class, X, y, train, test, params,
                        self.n_stop_rounds, self.max_iterations, self.scoring,
                        self.score_precision)
                       for train, test in cv.split(X)
                       for params in param_iter)
    
        n_splits = int(len(out)/len(param_iter))
    
        grid_scores = []
        i = 0
        best_mean = -np.inf
        for params in param_iter:
        
            score_list = []
            n_est_list = []
            for idx in range(n_splits):
                score, n_estimators = out[i]
                n_est_list.append(n_estimators)
                score_list.append(score)
                i += 1
        
            scores = np.array(score_list)
            mean = np.mean(scores)
            params['n_estimators'] = int(np.mean(n_est_list))
            score_tuple = _CVScoreTuple(params, mean, scores)
            grid_scores.append(_CVScoreTuple(params, mean, scores))
        
            if mean > best_mean:
                best_tuple = score_tuple
                best_params = params
        
        self.grid_scores_ = grid_scores
        self.best_score_tuple_ = best_tuple
        self.best_params_ = best_params
        
        if self.refit:
            self.best_params_['random_state'] = self.random_state
            gb = self._estimator_class(**self.best_params_)
            gb.fit(X, y)
            self.best_estimator_ = gb
    
        return self
    
    def predict(self, X):
        """Call `predict` on the best estimator.
    
        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.
        Returns
        -------
        y: array of shape = ["n_samples]
            The predicted values.
        """
        return self.best_estimator_.predict(X)

    def predict_proba(self, X):
        """Call `predict_proba` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Raises
        ------
        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        return self.best_estimator_.predict_proba(X)

    def predict_log_proba(self, X):
        """Call `predict_log_proba` on the best estimator.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.
        Raises
        ------

        AttributeError
            If the ``loss`` does not support probabilities.

        Returns
        -------
        p : array of shape = [n_samples]
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute `classes_`.
        """
        return self.best_estimator_.predict_log_proba(X)

    def _get_params(self):
        "Return parametrs to iterate over for underlying gradient boosting"


        params = {
            'loss': self.loss,
            'learning_rate': self.learning_rate,
            'subsample': self.subsample,
            'min_samples_split': self.min_samples_split,
            'min_samples_leaf': self.min_samples_leaf,
            'min_weight_fraction_leaf': self.min_weight_fraction_leaf,
            'max_depth': self.max_depth,
            'init': self.init,
            'max_features': self.max_features,
            'max_leaf_nodes': self.max_leaf_nodes,
            'presort': self.presort,
        }

        # Pre processing to ensure every parameter is a list
        # which `ParameterGrid` can iterate over
        for key, value in params.items():
            params[key] = np.atleast_1d(value)

        return params


def _fit_single_param(estimator, X, y, train, test, params, stop_rounds,
                      max_iter, scoring, score_precision):
    """ Fit a single estimator till stopping criteria is reached.

    Parameters
    ---------
    estimator : GradientBoostingRegressor or GradientBoostingClassifier
        The class of the model to instantiate.

    X : array-like, shape = [n_samples, n_features]
        Training vectors, where n_samples is the number of samples
        and n_features is the number of features.

    y : array-like, shape = [n_samples]
        Target values (integers in classification, real numbers in
        regression)
        For classification, labels must correspond to classes.

    train : array-like
        The indices of the training samples.

    test : array-like
        The indices of the testing samples.

    params : dict
        The parameters to instantiate the Gradient Boosting model with

    stop_rounds : int
        The number of iterations for which the score should remain constant
        over to stop training.

    max_iter : int
        The maximum number of estimators to fit during gradient boosting.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    score_precision : int
        The number of decimal places the score during each iterations is
        rounded off to before comparing with the previous ones.

    Returns
    -------
    score: float
        The score when the model converged.

    n_iterations:
        The number of iterations it took to converge. This is the same as
        parameter `n_estimators` for the gradient boosting model
    """

    params['warm_start'] = True
    gb = estimator(**params)
    scorer = check_scoring(estimator, scoring=scoring)
    scores = np.ones((stop_rounds,))

    X_train = X[train]
    y_train = y[train]
    X_test = X[test]
    y_test = y[test]

    for i in range(1, max_iter):
        gb.n_estimators = i
        gb.fit(X_train, y_train)

        scores = np.roll(scores, shift=-1)
        scores[-1] = scorer(gb, X_test, y_test)
        rscores = np.round(scores, score_precision)

        if np.allclose(rscores, rscores[-1]):
            break

    if i == max_iter:
        warnings.warn(str(gb) + ' failed to converge')

    return scores[-1], i
