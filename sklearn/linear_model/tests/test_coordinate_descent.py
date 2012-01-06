# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_almost_equal, \
                          assert_equal

from sklearn.linear_model.coordinate_descent import Lasso, \
    LassoCV, ElasticNet, ElasticNetCV


def test_lasso_zero():
    """Check that the lasso can handle zero data without crashing"""
    X = [[0], [0], [0]]
    y = [0, 0, 0]
    clf = Lasso(alpha=0).fit(X, y)
    pred = clf.predict([[1], [2], [3]])
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_, 0)


def test_lasso_toy():
    """
    Test Lasso on a toy example for various values of alpha.

    When validating this against glmnet notice that glmnet divides it
    against nobs.
    """

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line
    T = [[2], [3], [4]]  # test sample

    clf = Lasso(alpha=0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=0.1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.85])
    assert_array_almost_equal(pred, [1.7, 2.55, 3.4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.25])
    assert_array_almost_equal(pred, [0.5, 0.75, 1.])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_, 0)


def test_enet_toy():
    """
    Test ElasticNet for various parameters of alpha and rho.

    Actualy, the parameters alpha = 0 should not be alowed. However,
    we test it as a border case.

    ElasticNet is tested with and without precomputed Gram matrix
    """

    X = np.array([[-1.], [0.], [1.]])
    Y = [-1, 0, 1]       # just a straight line
    T = [[2.], [3.], [4.]]  # test sample

    # this should be the same as lasso
    clf = ElasticNet(alpha=0, rho=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.3, max_iter=100,
                     precompute=False)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=True)
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=np.dot(X.T, X))
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, rho=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090, 1.3636, 1.8181], 3)
    assert_almost_equal(clf.dual_gap_, 0)


def test_lasso_path():

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples
    n_samples, n_features, max_iter = 50, 200, 30
    random_state = np.random.RandomState(0)
    w = random_state.randn(n_features)
    w[10:] = 0.0  # only the top 10 features are impacting the model
    X = random_state.randn(n_samples, n_features)
    y = np.dot(X, w)

    clf = LassoCV(n_alphas=10, eps=1e-3, max_iter=max_iter).fit(X, y)
    assert_almost_equal(clf.alpha, 0.011, 2)

    clf = LassoCV(n_alphas=10, eps=1e-3, max_iter=max_iter, precompute=True)
    clf.fit(X, y)
    assert_almost_equal(clf.alpha, 0.011, 2)

    # test set
    X_test = random_state.randn(n_samples, n_features)
    y_test = np.dot(X_test, w)
    assert clf.score(X_test, y_test) > 0.85


def test_enet_path():

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples
    n_samples, n_features, max_iter = 50, 200, 50
    random_state = np.random.RandomState(0)
    w = random_state.randn(n_features)
    w[10:] = 0.0  # only the top 10 features are impacting the model
    X = random_state.randn(n_samples, n_features)
    y = np.dot(X, w)

    clf = ElasticNetCV(n_alphas=10, eps=1e-3, rho=0.95, cv=5,
            max_iter=max_iter)
    clf.fit(X, y)
    assert_almost_equal(clf.alpha, 0.002, 2)

    clf = ElasticNetCV(n_alphas=10, eps=1e-3, rho=0.95, cv=5,
                       max_iter=max_iter, precompute=True)
    clf.fit(X, y)
    assert_almost_equal(clf.alpha, 0.002, 2)

    # test set
    X_test = random_state.randn(n_samples, n_features)
    y_test = np.dot(X_test, w)
    assert clf.score(X_test, y_test) > 0.99


def test_path_parameters():

    # build an ill-posed linear regression problem with many noisy features and
    # comparatively few samples
    n_samples, n_features, max_iter = 50, 200, 50
    random_state = np.random.RandomState(0)
    w = random_state.randn(n_features)
    w[10:] = 0.0  # only the top 10 features are impacting the model
    X = random_state.randn(n_samples, n_features)
    y = np.dot(X, w)

    clf = ElasticNetCV(n_alphas=50, eps=1e-3, max_iter=max_iter,
                       rho=0.5)
    clf.fit(X, y)  # new params
    assert_almost_equal(0.5, clf.rho)
    assert_equal(50, clf.n_alphas)
    assert_equal(50, len(clf.alphas))


if __name__ == '__main__':
    import nose
    nose.runmodule()
