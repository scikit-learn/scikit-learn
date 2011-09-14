# Authors: Brian Holt, Peter Prettenhofer, Satrajit Ghosh
#
# License: BSD Style.

from __future__ import division
import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..utils import check_random_state

from ..tree import DecisionTreeRegressor


class MedianPredictor(object):

    median = None

    def fit(self, X, y):
        y = np.asanyarray(y)
        self.median = np.median(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.median)
        return y


class GradientBoostingClassifier(BaseEstimator, ClassifierMixin):

    def __init__(self, learn_rate=0.1, n_iter=100, min_split=5,
                 max_depth=4, init=None, random_state=None):
        if n_iter <= 0:
            raise ValueError("n_iter must be greater than 0")
        self.n_iter = n_iter

        if learn_rate <= 0.0:
            raise ValueError("learn_rate must be greater than 0")
        self.learn_rate = learn_rate

        self.min_split = min_split
        self.max_depth = max_depth

        if init == None:
            self.init = MedianPredictor()
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
            raise ValueError("Number of labels does not match number of samples.")

        self.init.fit(X, y)
        self.trees = [self.init]
        
        for i in xrange(self.n_iter):
            print "Boosting iteration %d" % i
            residual = y - self.predict(X)
            tree = DecisionTreeRegressor(min_split=self.min_split,
                                         max_depth=self.max_depth)
            tree.fit(X, residual)
            self.trees.append(tree)

    def predict(self, X):
        y = np.zeros((X.shape[0],), dtype=np.float64)
        for tree in self.trees:
            y += tree.predict(X)
        y[y > 0.0] = 1
        y[y <= 0.0] = -1
        return y
    
