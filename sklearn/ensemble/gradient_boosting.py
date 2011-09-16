# Authors: Peter Prettenhofer
#
# License: BSD Style.

from __future__ import division
import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..utils import check_random_state

from ..tree import DecisionTreeRegressor


class MedianPredictor(object):
    """A simple initial estimator that predicts the median
    of the training targets.
    """

    median = None

    def fit(self, X, y):
        y = np.asanyarray(y)
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
        y = np.asanyarray(y)
        self.mean = np.mean(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.mean)
        return y


class GradientBoostingBase(BaseEstimator):

    trees = []

    def __init__(self, learn_rate, n_iter, min_split, max_depth, init,
                 random_state):
        if n_iter <= 0:
            raise ValueError("n_iter must be greater than 0")
        self.n_iter = n_iter

        if learn_rate <= 0.0:
            raise ValueError("learn_rate must be greater than 0")
        self.learn_rate = learn_rate

        self.min_split = min_split
        self.max_depth = max_depth

        if init == 'median':
            self.init = MedianPredictor()
        elif init == 'mean':
            self.init = MeanPredictor()
        elif init == None:
            raise ValueError("init must not be None")
        else:
            if not hasattr(init, 'fit') or not hasattr(init, 'predict'):
                raise ValueError("init must be valid estimator")
            self.init = init

        self.random_state = check_random_state(random_state)

    def fit(self, X, y):
        X = np.asanyarray(X, dtype=np.float32, order='F')
        y = np.asanyarray(y, order='C')
        n_samples, n_features = X.shape
        if y.shape[0] != n_samples:
            raise ValueError("Number of labels does not match " \
                             "number of samples.")

        self.init.fit(X, y)
        self.trees = [self.init]
        y_pred = self._predict(X)
        for i in xrange(self.n_iter):
            #print "Boosting iteration %d" % i
            residual = y - y_pred
            tree = DecisionTreeRegressor(min_split=self.min_split,
                                         max_depth=self.max_depth)
            tree.fit(X, residual)
            self.trees.append(tree)
            y_pred = self._predict(X, old_pred=y_pred,
                                   learn_rate=self.learn_rate)

    def _predict(self, X, old_pred=None, learn_rate=1.0):
        if old_pred is not None:
            return old_pred + learn_rate * self.trees[-1].predict(X)
        else:
            y = np.zeros((X.shape[0],), dtype=np.float64)
            for i, tree in enumerate(self.trees):
                if i == 0:
                    y += tree.predict(X)
                else:
                    y += learn_rate * tree.predict(X)
            return y


class GradientBoostingClassifier(GradientBoostingBase, ClassifierMixin):

    def __init__(self, learn_rate=0.1, n_iter=100, min_split=5,
                 max_depth=4, init='median', random_state=None):

        super(GradientBoostingClassifier, self).__init__(
            learn_rate, n_iter, min_split, max_depth, init, random_state)

    def predict(self, X):
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        y = self._predict(X)
        y[y > 0.0] = 1.0
        y[y <= 0.0] = -1.0
        return y


class GradientBoostingRegressor(GradientBoostingBase, RegressorMixin):

    def __init__(self, learn_rate=0.1, n_iter=100, min_split=5,
                 max_depth=4, init='median', random_state=None):

        super(GradientBoostingRegressor, self).__init__(
            learn_rate, n_iter, min_split, max_depth, init, random_state)

    def predict(self, X):
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        y = self._predict(X, learn_rate=self.learn_rate)
        return y
