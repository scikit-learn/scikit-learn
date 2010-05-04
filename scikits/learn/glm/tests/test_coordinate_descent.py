# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
from numpy.testing import assert_array_almost_equal

from ..coordinate_descent import Lasso

# from ..coordinate_descent import ElasticNet
# from ..coordinate_descent import lasso_path
# from ..coordinate_descent import enet_path
# from ..coordinate_descent import enet_dual_gap, lasso_dual_gap


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

    clf = Lasso(alpha=0.1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.95])
    assert_array_almost_equal(pred, [1.9, 2.85, 3.8])

    clf = Lasso(alpha=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.75])
    assert_array_almost_equal(pred, [1.5, 2.25, 3.])

    clf = Lasso(alpha=1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.5])
    assert_array_almost_equal(pred, [1, 1.5, 2.])



# def test_lasso_cd_python_cython_sanity():
#     n_samples, n_features, maxit = 100, 50, 150
#     np.random.seed(0)
#     y = np.random.randn(n_samples)
#     X = np.random.randn(n_samples, n_features)

#     alpha = 1

#     model = Lasso(alpha=alpha)
#     model.fit(X, y, maxit=maxit)

    

#     model_fast_gap = lasso_dual_gap(X, y, model_fast.coef_, alpha)[0]

#     # check the convergence using the KKT condition
#     assert_array_almost_equal(model_fast_gap, 0, 4)

#     # check that python and cython implementations behave exactly the same
#     assert_array_almost_equal(model_slow.coef_, model_fast.coef_)
#     assert_array_almost_equal(model_slow_gap, model_fast_gap, 3)

#     # # check that the priori induces sparsity in the weights (feature selection)
#     assert_array_almost_equal(model.compute_density(), 0.88, 2)

# def test_enet_cd_python_cython_sanity():
#     n_samples, n_features, maxit = 100, 50, 150
#     np.random.seed(0)
#     y = np.random.randn(n_samples)
#     X = np.random.randn(n_samples, n_features)

#     alpha, beta = 1, 10

#     model_slow = ElasticNet(alpha=alpha, beta=beta)
#     model_slow.learner = enet_coordinate_descent_slow
#     model_slow.fit(X, y, maxit=maxit)

#     model_slow_gap = enet_dual_gap(X, y, model_slow.coef_, alpha, beta)[0]

#     # check the convergence using the KKT condition
#     assert_array_almost_equal(model_slow_gap, 0, 4)

#     model_fast = ElasticNet(alpha=alpha, beta=beta)
#     model_fast.learner = enet_coordinate_descent_fast
#     model_fast.fit(X, y, maxit=maxit)

#     model_fast_gap = enet_dual_gap(X, y, model_fast.coef_, alpha, beta)[0]

#     # check t convergence using the KKT condition
#     assert_array_almost_equal(model_fast_gap, 0, 4)

#     # check cython's sanity
#     assert_array_almost_equal(model_slow.coef_, model_fast.coef_)
#     assert_array_almost_equal(model_slow_gap, model_fast_gap, 3)

#     # check that the priori induces sparsity in the weights
#     # (feature selection) but not
#     assert_array_almost_equal(model_slow.compute_density(), 0.90, 2)


# def test_lasso_enet_cd_paths():
#     """Test Lasso and Elastic-Net path functions
#     """
#     n_samples, n_features, maxit = 5, 10, 30
#     np.random.seed(0)
#     y = np.random.randn(n_samples)
#     X = np.random.randn(n_samples, n_features)

#     alphas_lasso, weights_lasso = lasso_path(X, y, factor=0.97, n_alphas = 50,
#                                             tol=1e-2)
#     alphas_enet, weights_enet = enet_path(X, y, factor=0.97, n_alphas = 50,
#                                             beta=0.1, tol=1e-2)
