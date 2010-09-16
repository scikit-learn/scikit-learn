import numpy as np
from numpy.testing import (assert_array_almost_equal,
                           assert_almost_equal)

from nose.tools import assert_equal, assert_true

from ..lars import lars_path, LassoLARS, LARS
from ..coordinate_descent import Lasso

from scikits.learn import datasets

n, m = 10, 10
np.random.seed (0)
diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target


def test_simple():
    """
    Principle of LARS is to keep covariances tied and decreasing
    """
    max_pred = 10
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target, max_iter=max_pred, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < max_pred:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == max_pred


def test_simple_precomputed():
    """
    The same, with precomputed Gram matrix
    """
    max_pred = 10
    G = np.dot (diabetes.data.T, diabetes.data)
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target, Gram=G, max_iter=max_pred, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < max_pred:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == max_pred



def test_lars_lstsq():
    """
    Test that LARS gives least square solution at the end
    of the path
    """
    # test that it arrives to a least squares solution
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target, method="lar")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_path_.T[-1], coef_lstsq)


def test_lasso_gives_lstsq_solution():
    """
    Test that LARS Lasso gives least square solution at the end
    of the path
    """

    alphas_, active, coef_path_ = lars_path(X, y, max_iter=12, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq , coef_path_[:,-1])

def test_lasso_lars_vs_lasso_cd(verbose=False):
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results
    """
    lasso_lars = LassoLARS(alpha=0.1, normalize=False)
    lasso = Lasso(alpha=0.1, fit_intercept=False)
    for alpha in [0.1, 0.01, 0.004]:
        lasso_lars.alpha = alpha
        lasso_lars.fit(X, y, max_iter=12)
        lasso.alpha = alpha
        lasso.fit(X, y, maxit=5000, tol=1e-13)

        # make sure results are the same than with Lasso Coordinate descent
        error = np.linalg.norm(lasso_lars.coef_ - lasso.coef_)
        if verbose:
            print 'Error : ', error
        assert error < 1e-5


