from scikits.learn import datasets
from scikits.learn import glm

import numpy as np
# from numpy.testing import *


n, m = 10, 10
np.random.seed (0)
diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target


#normalize data
_xmean = X.mean(0)
_ymean = y.mean(0)
X = X - _xmean
y = y - _ymean
_norms = np.apply_along_axis (np.linalg.norm, 0, X)
nonzeros = np.flatnonzero(_norms)
X[:, nonzeros] /= _norms[nonzeros]


def test_1():
    """
    Principle of LARS is to keep covariances tied and decreasing
    """
    
    alphas_, active, coef_path_ = glm.lars_path(X, y, 6, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        assert ocur == i + 1

def test_lasso_lars_vs_lasso_cd():
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results
    """
    lasso_lars = glm.LassoLARS(alpha=0.1)
    lasso_lars.fit(X, y)

    # make sure results are the same than with Lasso Coordinate descent
    lasso = glm.Lasso(alpha=0.1)
    lasso.fit(X, y, maxit=3000, tol=1e-10)

    error = np.linalg.norm(lasso_lars.coef_ - lasso.coef_)
    assert error < 1e-5
