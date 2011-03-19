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


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
