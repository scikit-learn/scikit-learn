import numpy as np
from scipy import sparse
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_

from ..coordinate_descent import ElasticNet


def test_sparse_predict():
    """Check that the predict method works with dense coef_ and sparse X"""
    X = sparse.lil_matrix((3, 2))
    X[0, 0] = 1
    X[0, 1] = 0.5
    X[1, 0] = -1
    coef_ = np.array([1, -1])

    predicted = ElasticNet(coef_=coef_).predict(X)
    np.testing.assert_array_equal([0.5, -1.0, 0.0], predicted)


def test_enet_toy_list_input():
    """Test ElasticNet for various parameters of alpha and rho with list X"""

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line
    T = [[2], [3], [4]]  # test sample

    # this should be the same as lasso
    clf = ElasticNet(alpha=0, rho=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    #assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.3)
    clf.fit(X, Y, maxit=1000)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163,  1.5245,  2.0327], decimal=3)
    #assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090,  1.3636,  1.8181], 3)
    #assert_almost_equal(clf.dual_gap_, 0)


def test_enet_toy_explicit_sparse_input():
    """Test ElasticNet for various parameters of alpha and rho with sparse X"""

    # training samples
    X = sparse.lil_matrix((3, 1))
    X[0, 0] = -1
    # X[1, 0] = 0
    X[2, 0] = 1
    Y = [-1, 0, 1]       # just a straight line (the identity function)

    # test samples
    T = sparse.lil_matrix((3, 1))
    T[0, 0] = 2
    T[1, 0] = 3
    T[2, 0] = 4

    # this should be the same as lasso
    clf = ElasticNet(alpha=0, rho=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    #assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.3)
    clf.fit(X, Y, maxit=1000)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163,  1.5245,  2.0327], decimal=3)
    #assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090,  1.3636,  1.8181], 3)
    #assert_almost_equal(clf.dual_gap_, 0)


def test_sparse_enet_not_as_toy_dataset():

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples
    n_samples, n_features, maxit = 100, 100, 50
    np.random.seed(42)

    # generate a ground truth model
    w = np.random.randn(n_features)
    w[10:] = 0.0 # only the top 10 features are impacting the model

    X = np.random.randn(n_samples, n_features)
    rnd = np.random.uniform(size=(n_samples, n_features))
    X[rnd > 0.5] = 0.0 # 50% for zeros in input signal

    # generate training ground truth labels
    y = np.dot(X, w)

    X_train, X_test = X[:n_samples / 2], X[n_samples / 2:]
    y_train, y_test = y[:n_samples / 2], y[n_samples / 2:]

    clf = ElasticNet(alpha=0.1, rho=0.8).fit(X_train, y_train, maxit=maxit)
    assert_almost_equal(clf.score(X_test, y_test), 0.75, 1)


