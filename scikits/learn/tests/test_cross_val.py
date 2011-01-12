""" Test the cross_val module
"""

import numpy as np

import nose
from nose.tools import assert_true

from ..base import BaseEstimator
from ..datasets import load_iris
from ..metrics import zero_one_score
from ..cross_val import StratifiedKFold
from ..svm import SVC
from .. import cross_val
from ..cross_val import permutation_score


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


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    y = iris.target
    svm = SVC(kernel='linear')
    cv = StratifiedKFold(y, 2)

    score, scores, pvalue = permutation_score(svm, X, y, zero_one_score, cv)
    assert_true(score > 0.9)
    np.testing.assert_equal(pvalue,  0.0)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = permutation_score(svm, X, y, zero_one_score, cv)
    assert_true(score < 0.5)
    assert_true(pvalue > 0.4)
