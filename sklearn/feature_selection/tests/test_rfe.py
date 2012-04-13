"""
Testing Recursive feature elimination
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_true

from sklearn.feature_selection.rfe import RFE, RFECV
from sklearn.datasets import load_iris, make_classification, \
make_regression
from sklearn.metrics import zero_one
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier, SGDRegressor
from sklearn.utils import check_random_state


def test_rfe():
    generator = check_random_state(0)

    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    clf = SVC(kernel="linear")
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    X_r = rfe.transform(X)

    assert_true(X_r.shape == iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    assert_array_almost_equal(rfe.predict(X), clf.predict(iris.data))
    assert_true(rfe.score(X, y) == clf.score(iris.data, iris.target))


def test_rfecv():
    generator = check_random_state(0)

    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    # Test using the score function
    rfecv = RFECV(estimator=SVC(kernel="linear", C=100), step=1, cv=3)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert_true(X_r.shape == iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    # Test using a customized loss function
    rfecv = RFECV(estimator=SVC(kernel="linear", C=100), step=1, cv=3,
            loss_func=zero_one)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert_true(X_r.shape == iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])

# Test binary classification, multi-class classification, and regression
# with estimators that support warm_start

def test_binary_rfecv_warmstart():
    X, y = make_classification(n_samples=5000, n_features=4, n_informative=3, n_redundant=0, n_repeated=0, n_classes=2)
    clf = SGDClassifier(alpha=1, warm_start=False)

    rfecv = RFECV(estimator=clf, cv=10, step=1)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    print X_r.shape
    assert_true(X_r.shape[1] == 3)


def test_multiclass_rfecv_warmstart():
    X, y = make_classification(n_samples=3000, n_features=4, n_informative=3, n_redundant=0, n_repeated=0, n_classes=3)
    clf = SGDClassifier(alpha=1, warm_start=False)

    rfecv = RFECV(estimator=clf, cv=10, step=1)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    print X_r.shape
    assert_true(X_r.shape[1] == 3)



def test_regression_rfecv_warmstart():
    X, y = make_regression(n_samples=1000, n_features=4, n_informative=3)
    generator = check_random_state(0)
    clf = SGDRegressor(alpha=1, warm_start=False)

    rfecv = RFECV(estimator=clf, cv=10, step=0.1)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert_true(X_r.shape[1] == 3)
