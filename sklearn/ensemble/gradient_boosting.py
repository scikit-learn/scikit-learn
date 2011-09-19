# Authors: Peter Prettenhofer
#
# License: BSD Style.

from __future__ import division
import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..utils import check_random_state

from ..tree.tree import _build_tree
from ..tree.tree import _apply_tree
from ..tree._tree import MSE


def _predict_tree(tree, X):
    predictions = np.zeros(X.shape[0], dtype=np.float64)
    for idx, sample in enumerate(X):
        predictions[idx] = _apply_tree(tree, sample)
    return predictions


class MedianPredictor(object):
    """A simple initial estimator that predicts the median
    of the training targets.
    """

    median = None

    def fit(self, X, y):
        self.median = np.median(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.median)
        return y


class MeanPredictor(object):
    """A simple initial estimator that predicts the mean
    of the training targets.
    """

    mean = None

    def fit(self, X, y):
        self.mean = np.mean(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.mean)
        return y


class ClassPriorPredictor(object):
    """A simple initial estimator that predicts the mean
    of the training targets.
    """

    prior = None

    def fit(self, X, y):
        y_bar = y.mean()
        self.prior = 0.5 * np.log((1.0 + y_bar) / (1.0 - y_bar))
        print "prior", self.prior

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.prior)
        return y


class LossFunction(object):
    """Abstract base class for various loss functions."""

    def init_estimator(self, X, y):
        pass

    def loss(self, y, pred):
        pass

    def negative_gradient(self, y, pred):
        """Compute the negative gradient."""
        pass

    def update_terminal_regions(self, tree, X, y, residual, pred):
        """Update the terminal regions (=leafs) of the given
        tree. Traverses tree and invokes template method
        `_update_terminal_region`. """
        if tree.is_leaf:
            self._update_terminal_region(tree, X, y, residual, pred)
        else:
            self.update_terminal_regions(tree.left, X, y, residual, pred)
            self.update_terminal_regions(tree.right, X, y, residual, pred)

    def _update_terminal_region(self, node, X, y, residual, pred):
        """Template method for updating terminal regions (=leafs). """
        pass


class LeastSquaresError(LossFunction):
    """Loss function for least squares (LS) estimation. """

    def init_estimator(self):
        return MeanPredictor()

    def loss(self, y, pred):
        return 0.5 * np.sum((y - pred) ** 2.0)

    def negative_gradient(self, y, pred):
        return y - pred

    def update_terminal_regions(self, tree, X, y, residual, pred):
        """Terminal regions need not to be updated for least squares. """
        pass


class LeastAbsoluteError(LossFunction):
    """Loss function for least absolute deviation (LAD) regression. """

    def init_estimator(self):
        return MedianPredictor()

    def loss(self, y, pred):
        return np.abs(y - pred)

    def negative_gradient(self, y, pred):
        return np.sign(y - pred)

    def _update_terminal_region(self, node, X, y, residual, pred):
        """LAD updates terminal regions to median estimates. """
        node.value = np.median(y[node.sample_mask] - pred[node.sample_mask])


class BinomialDeviance(LossFunction):

    def init_estimator(self):
        return ClassPriorPredictor()

    def loss(self, y, pred):
        return np.log(1 + np.exp(-2.0 * y * pred))

    def negative_gradient(self, y, pred):
        return (2.0 * y) / (1.0 + np.exp(2.0 * y * pred))

    def _update_terminal_region(self, node, X, y, residual, pred):
        """Make a single Newton-Raphson step. """
        targets = residual[node.sample_mask]
        abs_targets = np.abs(targets)
        node.value = targets.sum() / np.sum(abs_targets * \
                                                 (2.0 - abs_targets))


LOSS_FUNCTIONS = {'ls': LeastSquaresError,
                  'lad': LeastAbsoluteError,
                  'deviance': BinomialDeviance}


class BaseGradientBoosting(BaseEstimator):

    trees = []

    def __init__(self, loss, learn_rate, n_iter, min_split, max_depth, init,
                 subsample, random_state):
        if n_iter <= 0:
            raise ValueError("n_iter must be greater than 0")
        self.n_iter = n_iter

        if learn_rate <= 0.0:
            raise ValueError("learn_rate must be greater than 0")
        self.learn_rate = learn_rate

        if loss not in LOSS_FUNCTIONS:
            raise ValueError("loss not supported")
        self.loss = LOSS_FUNCTIONS[loss]()

        if min_split <= 0:
            raise ValueError("min_split must be larger than 0.")
        self.min_split = min_split

        if subsample <= 0.0 or subsample > 1:
            raise ValueError("subsample must be in (0,1]")
        self.subsample = subsample

        if max_depth <= 0:
            raise ValueError("max_depth must be larger than 0.")
        self.max_depth = max_depth

        if init is None:
            self.init = self.loss.init_estimator()
        else:
            if not hasattr(init, 'fit') or not hasattr(init, 'predict'):
                raise ValueError("init must be valid estimator")
            self.init = init

        self.random_state = check_random_state(random_state)

    def fit(self, X, y):
        """Fit the gradient boosting model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes
            0, 1, ..., n_classes-1

        Returns
        -------
        self : object
            Returns self.
        """
        X = np.asanyarray(X, dtype=np.float32, order='F')
        y = np.asanyarray(y, order='C')
        n_samples, n_features = X.shape
        if y.shape[0] != n_samples:
            raise ValueError("Number of labels does not match " \
                             "number of samples.")

        # create argsorted X for fast tree induction
        X_argsorted = np.asfortranarray(
            np.argsort(X.T, axis=1).astype(np.int32).T)

        # fit initial model
        self.init.fit(X, y)

        # init predictions
        y_pred = self.init.predict(X)

        self.trees = []
        loss = self.loss

        # perform boosting iterations
        for i in xrange(self.n_iter):

            # subsampling
            sample_mask = np.random.rand(n_samples) > (1.0 - self.subsample)

            residual = loss.negative_gradient(y, y_pred)

            # induce regression tree on residuals
            tree = _build_tree(False, X, residual, MSE(), self.max_depth,
                               self.min_split, None, 1, self.random_state,
                               0.0, sample_mask, X_argsorted, True)

            loss.update_terminal_regions(tree, X, y, residual, y_pred)
            self.trees.append(tree)

            y_pred = self._predict(X, old_pred=y_pred,
                                   learn_rate=self.learn_rate)

        return self

    def _predict(self, X, old_pred=None, learn_rate=1.0):
        """Predict targets with current model. Re-uses predictions
        from previous iteration if available.
        """
        if old_pred is not None:
            return old_pred + learn_rate * _predict_tree(self.trees[-1], X)
        else:
            y = self.init.predict(X)
            for tree in self.trees:
                y += learn_rate * _predict_tree(tree, X)
            return y


class GradientBoostingClassifier(BaseGradientBoosting, ClassifierMixin):

    def __init__(self, loss='deviance', learn_rate=0.1, n_iter=100,
                 subsample=1.0, min_split=5, max_depth=4,
                 init=None, random_state=None):

        super(GradientBoostingClassifier, self).__init__(
            loss, learn_rate, n_iter, min_split, max_depth, init, subsample,
            random_state)

    def fit(self, X, y):
        self.classes = np.unique(y)
        if self.classes.shape[0] != 2:
            raise ValueError("only binary classification supported")
        y = np.searchsorted(self.classes, y)
        y[y == 0] = -1
        print "classes", self.classes
        super(GradientBoostingClassifier, self).fit(X, y)

    def predict(self, X):
        proba = self.predict_proba(X)
        return self.classes[(proba > 0.0).astype(np.int32)]

    def predict_proba(self, X):
        X = np.atleast_2d(X)
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        y = self._predict(X)
        return 1.0 / (1.0 + np.exp(-2.0 * y))


class GradientBoostingRegressor(BaseGradientBoosting, RegressorMixin):

    def __init__(self, loss='ls', learn_rate=0.1, n_iter=100, subsample=1.0,
                 min_split=5, max_depth=4, init=None, random_state=None):

        super(GradientBoostingRegressor, self).__init__(
            loss, learn_rate, n_iter, min_split, max_depth, init, subsample,
            random_state)

    def predict(self, X):
        X = np.atleast_2d(X)
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        y = self._predict(X, learn_rate=self.learn_rate)
        return y
