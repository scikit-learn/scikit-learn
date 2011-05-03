"""
Testing for grid search module (scikits.learn.grid_search)

"""
from nose.tools import assert_equal
from numpy.testing import assert_array_equal

import numpy as np
import scipy.sparse as sp

from scikits.learn.base import BaseEstimator
from scikits.learn.grid_search import GridSearchCV
from scikits.learn.datasets.samples_generator import test_dataset_classif
from scikits.learn.svm import LinearSVC
from scikits.learn.svm.sparse import LinearSVC as SparseLinearSVC
from scikits.learn.metrics import f1_score

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
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score

X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
y = np.array([1, 1, 2, 2])


def test_grid_search():
    """Test that the best estimator contains the right value for foo_param"""
    clf = MockClassifier()
    cross_validation = GridSearchCV(clf, {'foo_param': [1, 2, 3]})
    # make sure it selects the smallest parameter in case of ties
    assert_equal(cross_validation.fit(X, y).best_estimator.foo_param, 2)

    for i, foo_i in enumerate([1, 2, 3]):
        assert cross_validation.grid_scores_[i][0] == {'foo_param' : foo_i}


def test_grid_search_sparse():
    """Test that grid search works with both dense and sparse matrices"""
    X_, y_ = test_dataset_classif(n_samples=200, n_features=100, seed=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C':[0.1, 1.0]})
    cv.fit(X_[:180], y_[:180])
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator.C

    X_ = sp.csr_matrix(X_)
    clf = SparseLinearSVC()
    cv = GridSearchCV(clf, {'C':[0.1, 1.0]})
    cv.fit(X_[:180], y_[:180])
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator.C

    assert np.mean(y_pred == y_pred2) >= .9
    assert_equal(C, C2)


def test_grid_search_sparse_score_func():
    X_, y_ = test_dataset_classif(n_samples=200, n_features=100, seed=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, score_func=f1_score)
    # XXX: set refit to False due to a random bug when True (default)
    cv.fit(X_[:180], y_[:180], refit=False)
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator.C

    X_ = sp.csr_matrix(X_)
    clf = SparseLinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, score_func=f1_score)
    # XXX: set refit to False due to a random bug when True (default)
    cv.fit(X_[:180], y_[:180], refit=False)
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator.C

    assert_array_equal(y_pred, y_pred2)
    assert_equal(C, C2)


