"""Test the cross_validation module"""

import warnings
import numpy as np
from scipy.sparse import coo_matrix

from nose.tools import assert_true, assert_equal
from nose.tools import assert_raises
from sklearn.utils.testing import assert_greater, assert_less
from sklearn.utils.fixes import unique

from sklearn import cross_validation as cval
from sklearn.base import BaseEstimator
from sklearn.datasets import make_regression
from sklearn.datasets import load_iris
from sklearn.metrics import zero_one_score
from sklearn.metrics import f1_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import explained_variance_score
from sklearn.svm import SVC
from sklearn.linear_model import Ridge
from sklearn.svm.sparse import SVC as SparseSVC

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
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, cval.KFold, 3, 4)
    y = [0, 0, 1, 1, 2]
    assert_raises(ValueError, cval.StratifiedKFold, y, 3)

    # Check all indices are returned in the test folds
    kf = cval.KFold(300, 3)
    all_folds = None
    for train, test in kf:
        if all_folds is None:
            all_folds = test.copy()
        else:
            all_folds = np.concatenate((all_folds, test))

    all_folds.sort()
    assert_array_equal(all_folds, np.arange(300))


def test_shuffle_kfold():
    # Check the indices are shuffled properly, and that all indices are
    # returned in the different test folds
    kf1 = cval.KFold(300, 3, shuffle=True, random_state=0, indices=True)
    kf2 = cval.KFold(300, 3, shuffle=True, random_state=0, indices=False)
    ind = np.arange(300)

    for kf in (kf1, kf2):
        all_folds = None
        for train, test in kf:
            sorted_array = np.arange(100)
            assert np.any(sorted_array != ind[train])
            sorted_array = np.arange(101, 200)
            assert np.any(sorted_array != ind[train])
            sorted_array = np.arange(201, 300)
            assert np.any(sorted_array != ind[train])
            if all_folds is None:
                all_folds = ind[test].copy()
            else:
                all_folds = np.concatenate((all_folds, ind[test]))

        all_folds.sort()
        assert_array_equal(all_folds, ind)


def test_shuffle_split():
    ss1 = cval.ShuffleSplit(10, test_size=0.2, random_state=0)
    ss2 = cval.ShuffleSplit(10, test_size=2, random_state=0)
    ss3 = cval.ShuffleSplit(10, test_size=np.int32(2), random_state=0)
    ss4 = cval.ShuffleSplit(10, test_size=long(2), random_state=0)
    for t1, t2, t3, t4 in zip(ss1, ss2, ss3, ss4):
        assert_array_equal(t1[0], t2[0])
        assert_array_equal(t2[0], t3[0])
        assert_array_equal(t3[0], t4[0])
        assert_array_equal(t1[1], t2[1])
        assert_array_equal(t2[1], t3[1])
        assert_array_equal(t3[1], t4[1])


def test_stratified_shuffle_split():
    y = np.asarray([0, 1, 1, 1, 2, 2, 2])
    # Check that error is raised if there is a class with only one sample
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.2)

    # Check that error is raised if the test set size is smaller than n_classes
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 2)
    # Check that error is raised if the train set size is smaller than
    # n_classes
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 3, 2)

    y = np.asarray([0, 0, 0, 1, 1, 1, 2, 2, 2])
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.5, 0.6)
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 8, 0.6)
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.6, 8)

    ys = [
        np.array([1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3]),
        np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]),
        np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]),
        np.array([1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4]),
        np.array([-1] * 800 + [1] * 50)
        ]

    for y in ys:
        sss = cval.StratifiedShuffleSplit(y, 6, test_size=0.33,
                                          random_state=0, indices=True)
        for train, test in sss:
            assert_array_equal(unique(y[train]), unique(y[test]))
            # Checks if folds keep classes proportions
            p_train = np.bincount(
                unique(y[train], return_inverse=True)[1]
                ) / float(len(y[train]))
            p_test = np.bincount(
                unique(y[test], return_inverse=True)[1]
                ) / float(len(y[test]))
            assert_array_almost_equal(p_train, p_test, 1)
            assert_equal(y[train].size + y[test].size, y.size)
            assert_array_equal(np.lib.arraysetops.intersect1d(train, test), [])


def test_cross_val_score():
    clf = MockClassifier()
    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        scores = cval.cross_val_score(clf, X, y)
        assert_array_equal(scores, clf.score(X, y))

        scores = cval.cross_val_score(clf, X_sparse, y)
        assert_array_equal(scores, clf.score(X_sparse, y))


def test_train_test_split_errors():
    assert_raises(ValueError, cval.train_test_split)
    assert_raises(ValueError, cval.train_test_split, range(3),
            train_size=1.1)
    assert_raises(ValueError, cval.train_test_split, range(3),
            test_size=0.6, train_size=0.6)
    assert_raises(ValueError, cval.train_test_split, range(3),
            test_size=np.float32(0.6), train_size=np.float32(0.6))
    assert_raises(ValueError, cval.train_test_split, range(3),
            test_size="wrong_type")
    assert_raises(ValueError, cval.train_test_split, range(3),
            test_size=2, train_size=4)
    assert_raises(TypeError, cval.train_test_split, range(3),
            some_argument=1.1)
    assert_raises(ValueError, cval.train_test_split, range(3), range(42))


def test_shuffle_split_warnings():
    expected_message = ("test_fraction is deprecated in 0.11 and scheduled "
                        "for removal in 0.13, use test_size instead",
                        "train_fraction is deprecated in 0.11 and scheduled "
                        "for removal in 0.13, use train_size instead")

    with warnings.catch_warnings(record=True) as warn_queue:
        cval.ShuffleSplit(10, 3, test_fraction=0.1)
        cval.ShuffleSplit(10, 3, train_fraction=0.1)
        cval.train_test_split(range(3), test_fraction=0.1)
        cval.train_test_split(range(3), train_fraction=0.1)

    assert_equal(len(warn_queue), 4)
    assert_equal(str(warn_queue[0].message), expected_message[0])
    assert_equal(str(warn_queue[1].message), expected_message[1])
    assert_equal(str(warn_queue[2].message), expected_message[0])
    assert_equal(str(warn_queue[3].message), expected_message[1])


def test_train_test_split():
    X = np.arange(100).reshape((10, 10))
    X_s = coo_matrix(X)
    y = range(10)
    X_train, X_test, X_s_train, X_s_test, y_train, y_test = \
            cval.train_test_split(X, X_s, y)
    assert_array_equal(X_train, X_s_train.toarray())
    assert_array_equal(X_test, X_s_test.toarray())
    assert_array_equal(X_train[:, 0], y_train * 10)
    assert_array_equal(X_test[:, 0], y_test * 10)


def test_cross_val_score_with_score_func_classification():
    iris = load_iris()
    clf = SVC(kernel='linear')

    # Default score (should be the accuracy score)
    scores = cval.cross_val_score(clf, iris.data, iris.target, cv=5)
    assert_array_almost_equal(scores, [1., 0.97, 0.90, 0.97, 1.], 2)

    # Correct classification score (aka. zero / one score) - should be the
    # same as the default estimator score
    zo_scores = cval.cross_val_score(clf, iris.data, iris.target,
            score_func=zero_one_score, cv=5)
    assert_array_almost_equal(zo_scores, [1., 0.97, 0.90, 0.97, 1.], 2)

    # F1 score (class are balanced so f1_score should be equal to zero/one
    # score
    f1_scores = cval.cross_val_score(clf, iris.data, iris.target,
            score_func=f1_score, cv=5)
    assert_array_almost_equal(f1_scores, [1., 0.97, 0.90, 0.97, 1.], 2)


def test_cross_val_score_with_score_func_regression():
    X, y = make_regression(n_samples=30, n_features=20, n_informative=5,
                           random_state=0)
    reg = Ridge()

    # Default score of the Ridge regression estimator
    scores = cval.cross_val_score(reg, X, y, cv=5)
    assert_array_almost_equal(scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # R2 score (aka. determination coefficient) - should be the
    # same as the default estimator score
    r2_scores = cval.cross_val_score(reg, X, y, score_func=r2_score, cv=5)
    assert_array_almost_equal(r2_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # Mean squared error
    mse_scores = cval.cross_val_score(reg, X, y, cv=5,
            score_func=mean_squared_error)
    expected_mse = np.array([763.07, 553.16, 274.38, 273.26, 1681.99])
    assert_array_almost_equal(mse_scores, expected_mse, 2)

    # Explained variance
    ev_scores = cval.cross_val_score(reg, X, y, cv=5,
            score_func=explained_variance_score)
    assert_array_almost_equal(ev_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    X_sparse = coo_matrix(X)
    y = iris.target
    svm = SVC(kernel='linear')
    cv = cval.StratifiedKFold(y, 2)

    score, scores, pvalue = cval.permutation_test_score(
        svm, X, y, zero_one_score, cv)

    assert_greater(score, 0.9)
    np.testing.assert_almost_equal(pvalue, 0.0, 1)

    score_label, _, pvalue_label = cval.permutation_test_score(
        svm, X, y, zero_one_score, cv, labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # check that we obtain the same results with a sparse representation
    svm_sparse = SparseSVC(kernel='linear')
    cv_sparse = cval.StratifiedKFold(y, 2, indices=True)
    score_label, _, pvalue_label = cval.permutation_test_score(
        svm_sparse, X_sparse, y, zero_one_score, cv_sparse,
        labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = cval.permutation_test_score(svm, X, y,
            zero_one_score, cv)

    assert_less(score, 0.5)
    assert_greater(pvalue, 0.4)


def test_cross_val_generator_with_mask():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cval.LeaveOneOut(4, indices=False)
    lpo = cval.LeavePOut(4, 2, indices=False)
    kf = cval.KFold(4, 2, indices=False)
    skf = cval.StratifiedKFold(y, 2, indices=False)
    lolo = cval.LeaveOneLabelOut(labels, indices=False)
    lopo = cval.LeavePLabelOut(labels, 2, indices=False)
    ss = cval.ShuffleSplit(4, indices=False)
    for cv in [loo, lpo, kf, skf, lolo, lopo, ss]:
        for train, test in cv:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


def test_cross_val_generator_with_indices():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cval.LeaveOneOut(4, indices=True)
    lpo = cval.LeavePOut(4, 2, indices=True)
    kf = cval.KFold(4, 2, indices=True)
    skf = cval.StratifiedKFold(y, 2, indices=True)
    lolo = cval.LeaveOneLabelOut(labels, indices=True)
    lopo = cval.LeavePLabelOut(labels, 2, indices=True)
    b = cval.Bootstrap(2)  # only in index mode
    ss = cval.ShuffleSplit(2, indices=True)
    for cv in [loo, lpo, kf, skf, lolo, lopo, b, ss]:
        for train, test in cv:
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


def test_bootstrap_errors():
    assert_raises(ValueError, cval.Bootstrap, 10, train_size=100)
    assert_raises(ValueError, cval.Bootstrap, 10, test_size=100)
    assert_raises(ValueError, cval.Bootstrap, 10, train_size=1.1)
    assert_raises(ValueError, cval.Bootstrap, 10, test_size=1.1)


def test_shufflesplit_errors():
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=2.0)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=1.0)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=0.1,
            train_size=0.95)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=11)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=10)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=8,
            train_size=3)


def test_shufflesplit_reproducible():
    # Check that iterating twice on the ShuffleSplit gives the same
    # sequence of train-test when the random_state is given
    ss = cval.ShuffleSplit(10, random_state=21)
    assert_array_equal(list(a for a, b in ss), list(a for a, b in ss))


def test_cross_indices_exception():
    X = coo_matrix(np.array([[1, 2], [3, 4], [5, 6], [7, 8]]))
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cval.LeaveOneOut(4, indices=False)
    lpo = cval.LeavePOut(4, 2, indices=False)
    kf = cval.KFold(4, 2, indices=False)
    skf = cval.StratifiedKFold(y, 2, indices=False)
    lolo = cval.LeaveOneLabelOut(labels, indices=False)
    lopo = cval.LeavePLabelOut(labels, 2, indices=False)

    assert_raises(ValueError, cval.check_cv, loo, X, y)
    assert_raises(ValueError, cval.check_cv, lpo, X, y)
    assert_raises(ValueError, cval.check_cv, kf, X, y)
    assert_raises(ValueError, cval.check_cv, skf, X, y)
    assert_raises(ValueError, cval.check_cv, lolo, X, y)
    assert_raises(ValueError, cval.check_cv, lopo, X, y)
