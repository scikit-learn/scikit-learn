"""Test the cross_validation module"""

import numpy as np
from scipy.sparse import coo_matrix

from nose.tools import assert_true
from nose.tools import assert_raises

from ..base import BaseEstimator
from ..datasets import make_regression
from ..datasets import load_iris
from ..metrics import zero_one_score
from ..metrics import f1_score
from ..metrics import mean_square_error
from ..metrics import r2_score
from ..metrics import explained_variance_score
from ..cross_validation import StratifiedKFold
from ..svm import SVC
from ..linear_model import Ridge
from ..svm.sparse import SVC as SparseSVC
from .. import cross_validation
from ..cross_validation import permutation_test_score

from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal


class MockClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation"""

    def __init__(self, a=0):
        self.a = a

    def fit(self, X, Y):
        return self

    def predict(self, T):
        return T.shape[0]

    def score(self, X=None, Y=None):
        return 1. / (1 + np.abs(self.a))


X = np.ones((10, 2))
X_sparse = coo_matrix(X)
y = np.arange(10) / 2

##############################################################################
# Tests


def test_kfold():
    # Check that errors are raise if there is not enough samples
    assert_raises(AssertionError, cross_validation.KFold, 3, 4)
    y = [0, 0, 1, 1, 2]
    assert_raises(AssertionError, cross_validation.StratifiedKFold, y, 3)


def test_cross_val_score():
    clf = MockClassifier()
    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        scores = cross_validation.cross_val_score(clf, X, y)
        assert_array_equal(scores, clf.score(X, y))

        scores = cross_validation.cross_val_score(clf, X_sparse, y)
        assert_array_equal(scores, clf.score(X_sparse, y))


def test_cross_val_score_with_score_func_classification():
    iris = load_iris()
    clf = SVC(kernel='linear')

    # Default score (should be the accuracy score)
    scores = cross_validation.cross_val_score(clf, iris.data, iris.target,
                                              cv=5)
    assert_array_almost_equal(scores, [1., 0.97, 0.90, 0.97, 1.], 2)

    # Correct classification score (aka. zero / one score) - should be the
    # same as the default estimator score
    zo_scores = cross_validation.cross_val_score(clf, iris.data, iris.target,
                                          score_func=zero_one_score, cv=5)
    assert_array_almost_equal(zo_scores, [1., 0.97, 0.90, 0.97, 1.], 2)

    # F1 score (class are balanced so f1_score should be equal to zero/one
    # score
    f1_scores = cross_validation.cross_val_score(clf, iris.data, iris.target,
                                          score_func=f1_score, cv=5)
    assert_array_almost_equal(f1_scores, [1., 0.97, 0.90, 0.97, 1.], 2)


def test_cross_val_score_with_score_func_regression():
    X, y = make_regression(n_samples=30, n_features=20, n_informative=5,
                           random_state=0)
    reg = Ridge()

    # Default score of the Ridge regression estimator
    scores = cross_validation.cross_val_score(reg, X, y, cv=5)
    assert_array_almost_equal(scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # R2 score (aka. determination coefficient) - should be the
    # same as the default estimator score
    r2_scores = cross_validation.cross_val_score(reg, X, y,
                                                 score_func=r2_score, cv=5)
    assert_array_almost_equal(r2_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # Mean squared error
    mse_scores = cross_validation.cross_val_score(reg, X, y, cv=5,
                                           score_func=mean_square_error)
    expected_mse = [4578.47, 3319.02, 1646.29, 1639.58, 10092.00]
    assert_array_almost_equal(mse_scores, expected_mse, 2)

    # Explained variance
    ev_scores = cross_validation.cross_val_score(reg, X, y, cv=5,
                                          score_func=explained_variance_score)
    assert_array_almost_equal(ev_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    X_sparse = coo_matrix(X)
    y = iris.target
    svm = SVC(kernel='linear')
    cv = StratifiedKFold(y, 2)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, zero_one_score, cv)

    assert_true(score > 0.9)
    np.testing.assert_almost_equal(pvalue, 0.0, 1)

    score_label, _, pvalue_label = permutation_test_score(
        svm, X, y, zero_one_score, cv, labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # check that we obtain the same results with a sparse representation
    svm_sparse = SparseSVC(kernel='linear')
    cv_sparse = StratifiedKFold(y, 2, indices=True)
    score_label, _, pvalue_label = permutation_test_score(
        svm_sparse, X_sparse, y, zero_one_score, cv_sparse,
        labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, zero_one_score, cv)

    assert_true(score < 0.5)
    assert_true(pvalue > 0.4)


def test_cross_val_generator_with_mask():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cross_validation.LeaveOneOut(4, indices=False)
    lpo = cross_validation.LeavePOut(4, 2, indices=False)
    kf = cross_validation.KFold(4, 2, indices=False)
    skf = cross_validation.StratifiedKFold(y, 2, indices=False)
    lolo = cross_validation.LeaveOneLabelOut(labels, indices=False)
    lopo = cross_validation.LeavePLabelOut(labels, 2, indices=False)
    ss = cross_validation.ShuffleSplit(4, indices=False)
    for cv in [loo, lpo, kf, skf, lolo, lopo, ss]:
        for train, test in cv:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


def test_cross_val_generator_with_indices():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cross_validation.LeaveOneOut(4, indices=True)
    lpo = cross_validation.LeavePOut(4, 2, indices=True)
    kf = cross_validation.KFold(4, 2, indices=True)
    skf = cross_validation.StratifiedKFold(y, 2, indices=True)
    lolo = cross_validation.LeaveOneLabelOut(labels, indices=True)
    lopo = cross_validation.LeavePLabelOut(labels, 2, indices=True)
    b = cross_validation.Bootstrap(2)  # only in index mode
    ss = cross_validation.ShuffleSplit(2, indices=True)
    for cv in [loo, lpo, kf, skf, lolo, lopo, b, ss]:
        for train, test in cv:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


def test_bootstrap_errors():
    assert_raises(ValueError, cross_validation.Bootstrap, 10, n_train=100)
    assert_raises(ValueError, cross_validation.Bootstrap, 10, n_test=100)
    assert_raises(ValueError, cross_validation.Bootstrap, 10, n_train=1.1)
    assert_raises(ValueError, cross_validation.Bootstrap, 10, n_test=1.1)


def test_shufflesplit_errors():
    assert_raises(ValueError, cross_validation.ShuffleSplit, 10,
                  test_fraction=2.0)
    assert_raises(ValueError, cross_validation.ShuffleSplit, 10,
                  test_fraction=1.0)
    assert_raises(ValueError, cross_validation.ShuffleSplit, 10,
                  test_fraction=0.1, train_fraction=0.95)


def test_cross_indices_exception():
    X = coo_matrix(np.array([[1, 2], [3, 4], [5, 6], [7, 8]]))
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cross_validation.LeaveOneOut(4, indices=False)
    lpo = cross_validation.LeavePOut(4, 2, indices=False)
    kf = cross_validation.KFold(4, 2, indices=False)
    skf = cross_validation.StratifiedKFold(y, 2, indices=False)
    lolo = cross_validation.LeaveOneLabelOut(labels, indices=False)
    lopo = cross_validation.LeavePLabelOut(labels, 2, indices=False)

    assert_raises(ValueError, cross_validation.check_cv, loo, X, y)
    assert_raises(ValueError, cross_validation.check_cv, lpo, X, y)
    assert_raises(ValueError, cross_validation.check_cv, kf, X, y)
    assert_raises(ValueError, cross_validation.check_cv, skf, X, y)
    assert_raises(ValueError, cross_validation.check_cv, lolo, X, y)
    assert_raises(ValueError, cross_validation.check_cv, lopo, X, y)
