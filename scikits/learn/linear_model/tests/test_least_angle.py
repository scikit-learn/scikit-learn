import numpy as np
from numpy.testing import assert_array_almost_equal

from scikits.learn import linear_model, datasets

diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target

# TODO: use another dataset that has multiple drops


def test_simple():
    """
    Principle of LARS is to keep covariances tied and decreasing
    """

    alphas_, active, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar")

    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < X.shape[1]:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == X.shape[1]


def test_simple_precomputed():
    """
    The same, with precomputed Gram matrix
    """

    G = np.dot (diabetes.data.T, diabetes.data)
    alphas_, active, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, Gram=G, method="lar")

    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < X.shape[1]:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == X.shape[1]


def test_lars_lstsq():
    """
    Test that LARS gives least square solution at the end
    of the path
    """
    # test that it arrives to a least squares solution
    alphas_, active, coef_path_ = linear_model.lars_path(diabetes.data, diabetes.target,
                                                                method="lar")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_path_.T[-1], coef_lstsq)


def test_lasso_gives_lstsq_solution():
    """
    Test that LARS Lasso gives least square solution at the end
    of the path
    """

    alphas_, active, coef_path_ = linear_model.lars_path(X, y, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq , coef_path_[:,-1])


def test_collinearity():
    """Check that lars_path is robust to collinearity in input"""

    X = np.array([[3., 3., 1.],
                  [2., 2., 0.],
                  [1., 1., 0]])
    y = np.array([1., 0., 0])

    _, _, coef_path_ = linear_model.lars_path(X, y)
    assert (not np.isnan(coef_path_).any())
    assert_array_almost_equal(np.dot(X, coef_path_[:,-1]), y)


def test_singular_matrix():
    """
    Test when input is a singular matrix
    """
    X1 = np.array([[1, 1.], [1., 1.]])
    y1 = np.array([1, 1])
    alphas, active, coef_path = linear_model.lars_path(X1, y1)
    assert_array_almost_equal(coef_path.T, [[0, 0], [1, 0], [1, 0]])


def test_lasso_lars_vs_lasso_cd(verbose=False):
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results
    """
    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
    lasso_cd = linear_model.Lasso(fit_intercept=False)
    for (c, a) in zip(lasso_path.T, alphas):
        lasso_cd.alpha = a
        lasso_cd.fit(X, y, tol=1e-8)
        error = np.linalg.norm(c - lasso_cd.coef_)
        assert error < 0.01

def test_lasso_lars_vs_lasso_cd_early_stopping(verbose=False):
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results when early stopping is used.
    (test : before, in the middle, and in the last part of the path)
    """
    alphas_min = [10, 0.9, 1e-4]
    for alphas_min in alphas_min:
        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                    alpha_min=0.9)
        lasso_cd = linear_model.Lasso(fit_intercept=False)
        lasso_cd.alpha = alphas[-1]
        lasso_cd.fit(X, y, tol=1e-8)
        error = np.linalg.norm(lasso_path[:,-1] - lasso_cd.coef_)
        assert error < 0.01

def test_lars_add_features(verbose=False):
    linear_model.LARS(verbose=verbose, fit_intercept=True).fit(
        np.array([[ 0.02863763,  0.88144085, -0.02052429, -0.10648066, -0.06396584, -0.18338974],
                  [ 0.02038287,  0.51463335, -0.31734681, -0.12830467,  0.16870657, 0.02169503],
                  [ 0.14411476,  0.37666599,  0.2764702 ,  0.0723859 , -0.03812009, 0.03663579],
                  [-0.29411448,  0.33321005,  0.09429278, -0.10635334,  0.02827505, -0.07307312],
                  [-0.40929514,  0.57692643, -0.12559217,  0.19001991,  0.07381565, -0.0072319 ],
                  [-0.01763028,  1.        ,  0.04437242,  0.11870747,  0.1235008 , -0.27375014],
                  [-0.06482493,  0.1233536 ,  0.15686536,  0.02059646, -0.31723546, 0.42050836],
                  [-0.18806577,  0.01970053,  0.02258482, -0.03216307,  0.17196751, 0.34123213],
                  [ 0.11277307,  0.15590351,  0.11231502,  0.22009306,  0.1811108 , 0.51456405],
                  [ 0.03228484, -0.12317732, -0.34223564,  0.08323492, -0.15770904, 0.39392212],
                  [-0.00586796,  0.04902901,  0.18020746,  0.04370165, -0.06686751, 0.50099547],
                  [-0.12951744,  0.21978613, -0.04762174, -0.27227304, -0.02722684, 0.57449581]]),
        np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]))

if __name__ == '__main__':
    import nose
    nose.runmodule()

