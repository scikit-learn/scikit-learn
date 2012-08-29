import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_equal
import nose
from nose.tools import assert_raises, raises

from sklearn.utils.testing import assert_greater

from sklearn.linear_model import logistic
from sklearn import datasets

X = [[-1, 0], [0, 1], [1, 1]]
X_sp = sp.csr_matrix(X)
Y1 = [0, 1, 1]
Y2 = [2, 1, 0]
iris = datasets.load_iris()


def test_predict_2_classes():
    """Simple sanity check on a 2 classes dataset

    Make sure it predicts the correct result on simple datasets.
    """
    clf = logistic.LogisticRegression().fit(X, Y1)
    assert_array_equal(clf.predict(X), Y1)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression().fit(X_sp, Y1)
    assert_array_equal(clf.predict(X_sp), Y1)
    assert_array_equal(clf.predict_proba(X_sp).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression(C=100).fit(X, Y1)
    assert_array_equal(clf.predict(X), Y1)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression(C=100).fit(X_sp, Y1)
    assert_array_equal(clf.predict(X_sp), Y1)
    assert_array_equal(clf.predict_proba(X_sp).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression(fit_intercept=False).fit(X, Y1)
    assert_array_equal(clf.predict(X), Y1)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression(fit_intercept=False).fit(X_sp, Y1)
    assert_array_equal(clf.predict(X_sp), Y1)
    assert_array_equal(clf.predict_proba(X_sp).argmax(axis=1), Y1)


def test_error():
    """Test for appropriate exception on errors"""
    assert_raises(ValueError, logistic.LogisticRegression(C=-1).fit, X, Y1)


def test_predict_3_classes():
    clf = logistic.LogisticRegression(C=10).fit(X, Y2)
    assert_array_equal(clf.predict(X), Y2)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y2)

    clf = logistic.LogisticRegression(C=10).fit(X_sp, Y2)
    assert_array_equal(clf.predict(X_sp), Y2)
    assert_array_equal(clf.predict_proba(X_sp).argmax(axis=1), Y2)


def test_predict_iris():
    """Test logisic regression with the iris dataset"""

    clf = logistic.LogisticRegression(C=len(iris.data)).fit(iris.data,
                                                            iris.target)

    pred = clf.predict(iris.data)
    assert_greater(np.mean(pred == iris.target), .95)

    pred = clf.predict_proba(iris.data).argmax(axis=1)
    assert_greater(np.mean(pred == iris.target), .95)


def test_inconsistent_input():
    """Test that an exception is raised on inconsistent input"""
    rng = np.random.RandomState(0)
    X_ = rng.random_sample((5, 10))
    y_ = np.ones(X_.shape[0])
    y_[0] = 0

    clf = logistic.LogisticRegression()

    # Wrong dimensions for training data
    y_wrong = y_[:-1]
    assert_raises(ValueError, clf.fit, X, y_wrong)

    # Wrong dimensions for test data
    assert_raises(ValueError, clf.fit(X_, y_).predict,
            rng.random_sample((3, 12)))


@raises(ValueError)
def test_nan():
    """Test proper NaN handling.

    Regression test for Issue #252: fit used to go into an infinite loop.
    """
    Xnan = np.array(X, dtype=np.float64)
    Xnan[0, 1] = np.nan
    logistic.LogisticRegression().fit(Xnan, Y1)


if __name__ == '__main__':
    nose.runmodule()
