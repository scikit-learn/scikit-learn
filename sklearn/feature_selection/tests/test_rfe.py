"""
Testing Recursive feature elimination
"""

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from nose.tools import assert_equal
from scipy import sparse

from sklearn.feature_selection.rfe import RFE, RFECV
from sklearn.datasets import load_iris
from sklearn.metrics import zero_one
from sklearn.svm import SVC
from sklearn.utils import check_random_state


def test_rfe_set_params():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target
    clf = SVC(kernel="linear")
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    y_pred = rfe.fit(X, y).predict(X)

    clf = SVC()
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1,
              estimator_params={'kernel': 'linear'})
    y_pred2 = rfe.fit(X, y).predict(X)
    assert_array_equal(y_pred, y_pred2)


def test_rfe():
    generator = check_random_state(0)
    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    X_sparse = sparse.csr_matrix(X)
    y = iris.target

    # dense model
    clf = SVC(kernel="linear")
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    X_r = rfe.transform(X)
    clf.fit(X_r, y)
    assert_equal(len(rfe.ranking_), X.shape[1])

    # sparse model
    clf_sparse = SVC(kernel="linear")
    rfe_sparse = RFE(estimator=clf_sparse, n_features_to_select=4, step=0.1)
    rfe_sparse.fit(X_sparse, y)
    X_r_sparse = rfe_sparse.transform(X_sparse)

    assert_equal(X_r.shape, iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    assert_array_almost_equal(rfe.predict(X), clf.predict(iris.data))
    assert_equal(rfe.score(X, y), clf.score(iris.data, iris.target))
    assert_array_almost_equal(X_r, X_r_sparse.toarray())


def test_rfecv():
    generator = check_random_state(0)

    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = list(iris.target)   # regression test: list should be supported

    # Test using the score function
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=3)
    rfecv.fit(X, y)
    # non-regression test for missing worst feature:
    assert_equal(len(rfecv.cv_scores_), X.shape[1])
    assert_equal(len(rfecv.ranking_), X.shape[1])
    X_r = rfecv.transform(X)

    # same in sparse
    rfecv_sparse = RFECV(estimator=SVC(kernel="linear"), step=1, cv=3)
    X_sparse = sparse.csr_matrix(X)
    rfecv_sparse.fit(X_sparse, y)
    X_r_sparse = rfecv_sparse.transform(X_sparse)

    assert_equal(X_r.shape, iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])
    assert_array_almost_equal(X_r_sparse.toarray(), X_r)

    # Test using a customized loss function
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=3,
            loss_func=zero_one)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert_equal(X_r.shape, iris.data.shape)
    assert_array_almost_equal(X_r[:10], iris.data[:10])
