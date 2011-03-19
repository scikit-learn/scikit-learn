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
from ..cross_val import permutation_test_score


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

##############################################################################
# Tests

def test_kfold():
    # Check that errors are raise if there is not enough samples
    nose.tools.assert_raises(AssertionError, cross_val.KFold, 3, 4)
    y = [0, 0, 1, 1, 2]
    nose.tools.assert_raises(AssertionError, cross_val.StratifiedKFold, y, 3)


def test_cross_val_score():
    clf = MockClassifier()
    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        score = cross_val.cross_val_score(clf, X, y)
        np.testing.assert_array_equal(score, clf.score(X, y))


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    y = iris.target
    svm = SVC(kernel='linear')
    cv = StratifiedKFold(y, 2)

    score, scores, pvalue = permutation_test_score(svm, X, y,
                                                   zero_one_score, cv)
    assert_true(score > 0.9)
    np.testing.assert_almost_equal(pvalue, 0.0, 1)

    score_label, _, pvalue_label = permutation_test_score(svm, X, y,
                                                    zero_one_score,
                                                    cv, labels=np.ones(y.size),
                                                    rng=0)
    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = permutation_test_score(svm, X, y,
                                                   zero_one_score, cv)
    assert_true(score < 0.5)
    assert_true(pvalue > 0.4)


def test_cross_val_generator_with_indices():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cross_val.LeaveOneOut(4, indices=True)
    lpo = cross_val.LeavePOut(4, 2, indices=True)
    kf = cross_val.KFold(4, 2, indices=True)
    skf = cross_val.StratifiedKFold(y, 2, indices=True)
    lolo = cross_val.LeaveOneLabelOut(labels, indices=True)
    lopo = cross_val.LeavePLabelOut(labels, 2, indices=True)
    for cv in [loo, lpo, kf, skf, lolo, lopo]:
        for train, test in cv:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]
