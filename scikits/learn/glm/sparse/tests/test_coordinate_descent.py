import numpy as np
from scipy import sparse
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_

from scikits.learn.glm.sparse.coordinate_descent import Lasso as SparseLasso
from scikits.learn.glm.sparse.coordinate_descent import ElasticNet as SparseENet
from scikits.learn.glm.coordinate_descent import Lasso as DenseLasso
from scikits.learn.glm.coordinate_descent import ElasticNet as DenseENet


def test_sparse_predict():
    """Check that the predict method works with dense coef_ and sparse X"""
    X = sparse.lil_matrix((3, 2))
    X[0, 0] = 1
    X[0, 1] = 0.5
    X[1, 0] = -1
    coef_ = np.array([1, -1])

    predicted = SparseENet(coef_=coef_).predict(X)
    np.testing.assert_array_equal([0.5, -1.0, 0.0], predicted)


def test_lasso_zero():
    """Check that the sparse lasso can handle zero data without crashing"""
    X = sparse.csc_matrix((3, 1))
    y = [0, 0, 0]
    clf = SparseLasso().fit(X, y)
    pred = clf.predict([[1], [2], [3]])
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_,  0)


def test_enet_toy_list_input():
    """Test ElasticNet for various parameters of alpha and rho with list X"""

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line
    T = [[2], [3], [4]]  # test sample

    # this should be the same as unregularized least squares
    clf = SparseENet(alpha=0, rho=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = SparseENet(alpha=0.5, rho=0.3)
    clf.fit(X, Y, maxit=1000)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163,  1.5245,  2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf = SparseENet(alpha=0.5, rho=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090,  1.3636,  1.8181], 3)
    assert_almost_equal(clf.dual_gap_, 0)


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
    clf = SparseENet(alpha=0, rho=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = SparseENet(alpha=0.5, rho=0.3)
    clf.fit(X, Y, maxit=1000)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163,  1.5245,  2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf = SparseENet(alpha=0.5, rho=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090,  1.3636,  1.8181], 3)
    assert_almost_equal(clf.dual_gap_, 0)


def make_sparse_data(n_samples, n_features, n_informative, seed=42):
    rng = np.random.RandomState(seed)

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples

    # generate a ground truth model
    w = rng.randn(n_features)
    w[n_informative:] = 0.0 # only the top features are impacting the model

    X = rng.randn(n_samples, n_features)
    rnd = rng.uniform(size=(n_samples, n_features))
    X[rnd > 0.5] = 0.0 # 50% of zeros in input signal

    # generate training ground truth labels
    y = np.dot(X, w)
    return X, y


def test_sparse_enet_not_as_toy_dataset():
    n_samples, n_features, maxit = 100, 100, 1000
    n_informative = 10

    X, y = make_sparse_data(n_samples, n_features, n_informative)

    X_train, X_test = X[n_samples / 2:], X[:n_samples / 2]
    y_train, y_test = y[n_samples / 2:], y[:n_samples / 2]

    s_clf = SparseENet(alpha=0.1, rho=0.8, fit_intercept=False)
    s_clf.fit(X_train, y_train, maxit=maxit, tol=1e-7)
    assert_almost_equal(s_clf.dual_gap_, 0, 4)
    assert_(s_clf.score(X_test, y_test) > 0.85)

    # check the convergence is the same as the dense version
    d_clf = DenseENet(alpha=0.1, rho=0.8, fit_intercept=False)
    d_clf.fit(X_train, y_train, maxit=maxit, tol=1e-7)
    assert_almost_equal(d_clf.dual_gap_, 0, 4)
    assert_(d_clf.score(X_test, y_test) > 0.85)

    assert_almost_equal(s_clf.coef_, d_clf.coef_, 5)

    # check that the coefs are sparse
    assert_(np.sum(s_clf.coef_ != 0.0) < 2 * n_informative)


def test_sparse_lasso_not_as_toy_dataset():
    n_samples, n_features, maxit = 100, 100, 1000
    n_informative = 10

    X, y = make_sparse_data(n_samples, n_features, n_informative)

    X_train, X_test = X[n_samples / 2:], X[:n_samples / 2]
    y_train, y_test = y[n_samples / 2:], y[:n_samples / 2]

    s_clf = SparseLasso(alpha=0.1, fit_intercept=False)
    s_clf.fit(X_train, y_train, maxit=maxit, tol=1e-7)
    assert_almost_equal(s_clf.dual_gap_, 0, 4)
    assert_(s_clf.score(X_test, y_test) > 0.85)

    # check the convergence is the same as the dense version
    d_clf = DenseLasso(alpha=0.1, fit_intercept=False)
    d_clf.fit(X_train, y_train, maxit=maxit, tol=1e-7)
    assert_almost_equal(d_clf.dual_gap_, 0, 4)
    assert_(d_clf.score(X_test, y_test) > 0.85)

    # check that the coefs are sparse
    assert_equal(np.sum(s_clf.coef_ != 0.0), n_informative)


