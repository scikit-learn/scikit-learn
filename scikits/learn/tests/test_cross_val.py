""" Test the cross_val module
"""

import numpy as np

import nose

from ..base import BaseEstimator
from .. import cross_val

class MockClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation
    
    """
    def __init__(self, a=0):
        self.a = a

    def fit(self, X, Y, **params):
        self._set_params(**params)
        return self

    def predict(self, T):
        return T.shape[0]

    def score(self, X=None, Y=None):
        return 1./(1+np.abs(self.a))


X = np.ones((10, 2))
y = np.arange(10)/2

################################################################################
# Tests

def test_kfold():
    # Check that errors are raise if there is not enough samples
    nose.tools.assert_raises(AssertionError, cross_val.KFold, 3, 3)
    y = [0, 0, 1, 1, 2]
    nose.tools.assert_raises(AssertionError, cross_val.StratifiedKFold, y, 3)


def test_cross_val_score():
    clf = MockClassifier()
    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        score = cross_val.cross_val_score(clf, X, y)
        np.testing.assert_array_equal(score,  clf.score(X, y))

