import numpy as np
from numpy.testing import assert_array_almost_equal

# from nose.tools import assert_equal, assert_true

from ..lars import lars_path, LassoLARS
from ..coordinate_descent import Lasso

from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target

# TODO: use another dataset that has multiple drops


def test_simple():
    """
    Principle of LARS is to keep covariances tied and decreasing
    """
    max_features = 10
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target,
                                        max_features=max_features, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < max_features:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == max_features


def test_simple_precomputed():
    """
    The same, with precomputed Gram matrix
    """
    max_features = 10
    G = np.dot (diabetes.data.T, diabetes.data)
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target,
                                Gram=G, max_features=max_features, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < max_features:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == max_features


def test_lars_lstsq():
    """
    Test that LARS gives least square solution at the end
    of the path
    """
    # test that it arrives to a least squares solution
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target,
                                                                method="lar")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_path_.T[-1], coef_lstsq)


def test_lasso_gives_lstsq_solution():
    """
    Test that LARS Lasso gives least square solution at the end
    of the path
    """

    alphas_, active, coef_path_ = lars_path(X, y, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq , coef_path_[:,-1])


def test_lasso_lars_vs_lasso_cd(verbose=False):
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results
    """
    alphas, _, lasso_path = lars_path(X, y, method='lasso')
    alphas /= X.shape[0]
    lasso_cd = Lasso(fit_intercept=False)
    for (c, a) in zip(lasso_path.T, alphas):
        lasso_cd.alpha = a
        lasso_cd.fit(X, y, tol=1e-8)
        error = np.linalg.norm(c - lasso_cd.coef_)
        assert error < 0.01


