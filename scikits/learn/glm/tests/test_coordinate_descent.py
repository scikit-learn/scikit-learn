# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
from numpy.testing import assert_array_almost_equal

from ..coordinate_descent import Lasso, lasso_path
from ..coordinate_descent import ElasticNet, enet_path

def test_Lasso_toy():
    """
    test predict on a toy example.

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
    assert clf.dual_gap_ == 0.

    clf = Lasso(alpha=0.1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.95])
    assert_array_almost_equal(pred, [1.9, 2.85, 3.8])
    assert clf.dual_gap_ == 0.

    clf = Lasso(alpha=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.75])
    assert_array_almost_equal(pred, [1.5, 2.25, 3.])
    assert clf.dual_gap_ == 0.

    clf = Lasso(alpha=1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.5])
    assert_array_almost_equal(pred, [1, 1.5, 2.])
    assert clf.dual_gap_ == 0.


def test_Enet_toy():
    """
    Test ElasticNet for various parameters of alpha and beta.

    Actualy, the parameters alpha = 0 should not be alowed. However,
    we test it as a border case.
    """

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line
    T = [[2], [3], [4]]  # test sample

    # this should be the same as lasso
    clf = ElasticNet(alpha=0, beta=0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert clf.dual_gap_ == 0.

    clf = ElasticNet(alpha=0.5, beta=1.)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.5])
    assert_array_almost_equal(pred, [1, 1.5, 2.])
    assert 0 < clf.dual_gap_ < 1e-10

def test_lasso_path():
    """
    Test for the complete lasso path.

    As the weigths_lasso array is quite big, we only test at the first
    & last index.
    """
    n_samples, n_features, maxit = 5, 10, 30
    np.random.seed(0)
    Y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    alphas_lasso, weights_lasso = lasso_path(X, Y, factor=0.97,
                                             n_alphas = 10, tol=1e-3)
    assert_array_almost_equal(alphas_lasso,
                              [ 4.498, 4.363, 4.232, 4.105, 3.982,
                              3.863, 3.747, 3.634, 3.525, 3.420],
                              decimal=3)

    assert weights_lasso.shape == (10, 10)

    assert_array_almost_equal(weights_lasso[0],
                              [0, 0, 0, 0, 0 , -0.016, 0, 0, 0, 0],
                              decimal=3)

    assert_array_almost_equal(weights_lasso[9],
                              [-0.038, 0, 0, 0, 0, -0.148, 0, -0.095, 0, 0],
                              decimal=3)

def test_enet_path():
    n_samples, n_features, maxit = 5, 10, 30
    np.random.seed(0)
    Y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    alphas_enet, weights_enet = enet_path(X, Y, factor=0.97,
                                             n_alphas = 10, tol=1e-3)
    assert_array_almost_equal(alphas_enet,
                              [ 4.498, 4.363, 4.232, 4.105, 3.982,
                              3.863, 3.747, 3.634, 3.525, 3.420],
                              decimal=3)

    assert weights_enet.shape == (10, 10)

    assert_array_almost_equal(weights_enet[0],
                              [0, 0, 0, 0, 0 , -0.016, 0, 0, 0, 0],
                              decimal=3)

    assert_array_almost_equal(weights_enet[9],
                              [-0.028, 0, 0, 0, 0, -0.131, 0, -0.081, 0, 0],
                              decimal=3)
