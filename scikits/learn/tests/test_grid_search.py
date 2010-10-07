"""
Testing for grid search module (scikits.learn.grid_search)

"""
from nose.tools import assert_equal

import numpy as np
from scikits.learn.base import BaseEstimator
from scikits.learn.grid_search import GridSearchCV

class MockClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation
    
    """
    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y, **params):
        self._set_params(**params)
        return self

    def predict(self, T):
        return T.shape[0]

    def score(self, X=None, Y=None):
        if self.foo_param == 2:
            score = 1.
        else:
            score = 0.
        return score

X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
y = np.array([1, 1, 2, 2])


def test_GridSearch():
    clf = MockClassifier()
    cross_validation = GridSearchCV(clf, {'foo_param':[1, 2, 3]})
    assert_equal(cross_validation.fit(X, y).best_estimator.foo_param, 2)
