import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn import linear_model, datasets

diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target

# TODO: use another dataset that has multiple drops


def test_simple():
    """
    Principle of Lars is to keep covariances tied and decreasing
    """

    alphas_, active, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar")

    for (i, coef_) in enumerate(coef_path_.T):
        res = y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[C - eps < abs(cov)])
        if i < X.shape[1]:
            assert ocur == i + 1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == X.shape[1]


def test_simple_precomputed():
    """
    The same, with precomputed Gram matrix
    """

    G = np.dot(diabetes.data.T, diabetes.data)
    alphas_, active, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, Gram=G, method="lar")

    for i, coef_ in enumerate(coef_path_.T):
        res = y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[C - eps < abs(cov)])
        if i < X.shape[1]:
            assert ocur == i + 1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == X.shape[1]


def test_lars_lstsq():
    """
    Test that Lars gives least square solution at the end
    of the path
    """
    X1 = 3 * diabetes.data  # use un-normalized dataset
    clf = linear_model.LassoLars(alpha=0.)
    clf.fit(X1, y)
    coef_lstsq = np.linalg.lstsq(X1, y)[0]
    assert_array_almost_equal(clf.coef_, coef_lstsq)


def test_lasso_gives_lstsq_solution():
    """
    Test that Lars Lasso gives least square solution at the end
    of the path
    """
    alphas_, active, coef_path_ = linear_model.lars_path(X, y, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq, coef_path_[:, -1])


def test_collinearity():
    """Check that lars_path is robust to collinearity in input"""
    X = np.array([[3., 3., 1.],
                  [2., 2., 0.],
                  [1., 1., 0]])
    y = np.array([1., 0., 0])

    _, _, coef_path_ = linear_model.lars_path(X, y)
    assert (not np.isnan(coef_path_).any())
    assert_array_almost_equal(np.dot(X, coef_path_[:, -1]), y)


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
    X = 3 * diabetes.data

    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
    lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
    for c, a in zip(lasso_path.T, alphas):
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
        error = np.linalg.norm(c - lasso_cd.coef_)
        assert error < 0.01

    # similar test, with the classifiers
    for alpha in np.linspace(1e-2, 1 - 1e-2):
        clf1 = linear_model.LassoLars(alpha=alpha, normalize=False).fit(X, y)
        clf2 = linear_model.Lasso(alpha=alpha, tol=1e-8,
                                  normalize=False).fit(X, y)
        err = np.linalg.norm(clf1.coef_ - clf2.coef_)
        assert err < 1e-3

    # same test, with normalized data
    X = diabetes.data
    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
    lasso_cd = linear_model.Lasso(fit_intercept=False, normalize=True,
                                  tol=1e-8)
    for c, a in zip(lasso_path.T, alphas):
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
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
        lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
        lasso_cd.alpha = alphas[-1]
        lasso_cd.fit(X, y)
        error = np.linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
        assert error < 0.01

    alphas_min = [10, 0.9, 1e-4]
    # same test, with normalization
    for alphas_min in alphas_min:
        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                    alpha_min=0.9)
        lasso_cd = linear_model.Lasso(fit_intercept=True, normalize=True,
                                      tol=1e-8)
        lasso_cd.alpha = alphas[-1]
        lasso_cd.fit(X, y)
        error = np.linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
        assert error < 0.01



def test_lars_add_features(verbose=False):
    """
    assure that at least some features get added if necessary

    test for 6d2b4c
    """
    linear_model.Lars(verbose=verbose, fit_intercept=True).fit(
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


def test_lars_n_nonzero_coefs(verbose=False):
    lars = linear_model.Lars(n_nonzero_coefs=6, verbose=verbose)
    lars.fit(X, y)
    assert len(lars.coef_.nonzero()[0]) == 6


def test_lars_cv():
    """ Test the LassoLarsCV object by checking that the optimal alpha
        increases as the number of samples increases.

        This property is not actualy garantied in general and is just a
        property of the given dataset, with the given steps chosen.
    """
    old_alpha = 0
    lars_cv = linear_model.LassoLarsCV()
    for length in (400, 200, 100):
        X = diabetes.data[:length]
        y = diabetes.target[:length]
        lars_cv.fit(X, y)
        np.testing.assert_array_less(old_alpha, lars_cv.alpha)
        old_alpha = lars_cv.alpha


def test_lasso_lars_ic():
    """ Test the LassoLarsIC object by checking that
        - some good features are selected.
        - alpha_bic > alpha_aic
        - n_nonzero_bic < n_nonzero_aic
    """
    lars_bic = linear_model.LassoLarsIC('bic')
    lars_aic = linear_model.LassoLarsIC('aic')
    rng = np.random.RandomState(42)
    X = diabetes.data
    y = diabetes.target
    X = np.c_[X, rng.randn(X.shape[0], 4)]  # add 4 bad features
    lars_bic.fit(X, y)
    lars_aic.fit(X, y)
    nonzero_bic = np.where(lars_bic.coef_)[0]
    nonzero_aic = np.where(lars_aic.coef_)[0]
    assert lars_bic.alpha_ > lars_aic.alpha_
    assert len(nonzero_bic) < len(nonzero_aic)
    assert np.max(nonzero_bic) < diabetes.data.shape[1]


if __name__ == '__main__':
    import nose
    nose.runmodule()
