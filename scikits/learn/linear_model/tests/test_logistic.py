import numpy as np

from numpy.testing import assert_array_equal, \
     assert_array_almost_equal, assert_almost_equal
import nose
from nose.tools import assert_raises

from scikits.learn.linear_model import logistic
from scikits.learn import datasets

X = [[-1, 0], [0, 1], [1, 1]]
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

    clf = logistic.LogisticRegression(C=100).fit(X, Y1)
    assert_array_equal(clf.predict(X), Y1)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y1)

    clf = logistic.LogisticRegression(fit_intercept=False).fit(X, Y1)
    assert_array_equal(clf.predict(X), Y1)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y1)


def test_error():
    """Test for appropriate exception on errors"""
    assert_raises (ValueError, logistic.LogisticRegression(C=-1).fit, X, Y1)


def test_predict_3_classes():
    clf = logistic.LogisticRegression(C=10).fit(X, Y2)
    assert_array_equal(clf.predict(X), Y2)
    assert_array_equal(clf.predict_proba(X).argmax(axis=1), Y2)

def test_predict_iris():
    """Test logisic regression with the iris dataset"""

    clf = logistic.LogisticRegression().fit(iris.data, iris.target)

    pred = clf.predict(iris.data)
    assert np.mean(pred == iris.target) > .95

    pred = clf.predict_proba(iris.data).argmax(axis=1)
    assert np.mean(pred == iris.target) > .95

def test_inconsistent_input():
    """Test that an exception is raised when input to predict is inconsistent"""
    X_ = np.random.random((5, 10))
    y_ = np.ones(X_.shape[0])
    assert_raises(ValueError,
                  logistic.LogisticRegression().fit(X_, y_).predict,
                  np.random.random((3,12)))

def test_transform():
    clf = logistic.LogisticRegression(penalty="l1")
    clf.fit(iris.data, iris.target)
    X_new = clf.transform(iris.data)
    clf = logistic.LogisticRegression()
    clf.fit(X_new, iris.target)
    pred = clf.predict(X_new)
    assert np.mean(pred == iris.target) >= 0.75

if __name__ == '__main__':
    import nose
    nose.runmodule()

