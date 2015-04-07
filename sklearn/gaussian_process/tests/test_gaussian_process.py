"""
Testing for Gaussian Process module (sklearn.gaussian_process)
"""

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
# Licence: BSD 3 clause

from nose.tools import raises
from nose.tools import assert_true

import numpy as np

from sklearn.gaussian_process import GaussianProcess
from sklearn.gaussian_process import regression_models as regression
from sklearn.gaussian_process import correlation_models as correlation
from sklearn.datasets import make_regression
from sklearn.utils.testing import assert_greater


f = lambda x: x * np.sin(x)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = f(X).ravel()


def test_1d(regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    # MLE estimation of a one-dimensional Gaussian Process model.
    # Check random start optimization.
    # Test the interpolating property.
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
                         random_start=random_start, verbose=False).fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)
    y2_pred, MSE2 = gp.predict(X2, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.)
                and np.allclose(MSE2, 0., atol=10))


def test_2d(regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    # MLE estimation of a two-dimensional Gaussian Process model accounting for
    # anisotropy. Check random start optimization.
    # Test the interpolating property.
    b, kappa, e = 5., .5, .1
    g = lambda x: b - x[:, 1] - kappa * (x[:, 0] - e) ** 2.
    X = np.array([[-4.61611719, -6.00099547],
                  [4.10469096, 5.32782448],
                  [0.00000000, -0.50000000],
                  [-6.17289014, -4.6984743],
                  [1.3109306, -6.93271427],
                  [-5.03823144, 3.10584743],
                  [-2.87600388, 6.74310541],
                  [5.21301203, 4.26386883]])
    y = g(X).ravel()

    thetaL = [1e-4] * 2
    thetaU = [1e-1] * 2
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=[1e-2] * 2, thetaL=thetaL,
                         thetaU=thetaU,
                         random_start=random_start, verbose=False)
    gp.fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.))

    assert_true(np.all(gp.theta_ >= thetaL))  # Lower bounds of hyperparameters
    assert_true(np.all(gp.theta_ <= thetaU))  # Upper bounds of hyperparameters


def test_2d_2d(regr=regression.constant, corr=correlation.squared_exponential,
               random_start=10, beta0=None):
    # MLE estimation of a two-dimensional Gaussian Process model accounting for
    # anisotropy. Check random start optimization.
    # Test the GP interpolation for 2D output
    b, kappa, e = 5., .5, .1
    g = lambda x: b - x[:, 1] - kappa * (x[:, 0] - e) ** 2.
    f = lambda x: np.vstack((g(x), g(x))).T
    X = np.array([[-4.61611719, -6.00099547],
                  [4.10469096, 5.32782448],
                  [0.00000000, -0.50000000],
                  [-6.17289014, -4.6984743],
                  [1.3109306, -6.93271427],
                  [-5.03823144, 3.10584743],
                  [-2.87600388, 6.74310541],
                  [5.21301203, 4.26386883]])
    y = f(X)
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=[1e-2] * 2, thetaL=[1e-4] * 2,
                         thetaU=[1e-1] * 2,
                         random_start=random_start, verbose=False)
    gp.fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.))


@raises(ValueError)
def test_wrong_number_of_outputs():
    gp = GaussianProcess()
    gp.fit([[1, 2, 3], [4, 5, 6]], [1, 2, 3])


def test_more_builtin_correlation_models(random_start=1):
    # Repeat test_1d and test_2d for several built-in correlation
    # models specified as strings.
    all_corr = ['absolute_exponential', 'squared_exponential', 'cubic',
                'linear']

    for corr in all_corr:
        test_1d(regr='constant', corr=corr, random_start=random_start)
        test_2d(regr='constant', corr=corr, random_start=random_start)
        test_2d_2d(regr='constant', corr=corr, random_start=random_start)


def test_ordinary_kriging():
    # Repeat test_1d and test_2d with given regression weights (beta0) for
    # different regression models (Ordinary Kriging).
    test_1d(regr='linear', beta0=[0., 0.5])
    test_1d(regr='quadratic', beta0=[0., 0.5, 0.5])
    test_2d(regr='linear', beta0=[0., 0.5, 0.5])
    test_2d(regr='quadratic', beta0=[0., 0.5, 0.5, 0.5, 0.5, 0.5])
    test_2d_2d(regr='linear', beta0=[0., 0.5, 0.5])
    test_2d_2d(regr='quadratic', beta0=[0., 0.5, 0.5, 0.5, 0.5, 0.5])


def test_no_normalize():
    gp = GaussianProcess(normalize=False).fit(X, y)
    y_pred = gp.predict(X)
    assert_true(np.allclose(y_pred, y))


def test_random_starts():
    # Test that an increasing number of random-starts of GP fitting only
    # increases the reduced likelihood function of the optimal theta.
    n_samples, n_features = 50, 3
    np.random.seed(0)
    rng = np.random.RandomState(0)
    X = rng.randn(n_samples, n_features) * 2 - 1
    y = np.sin(X).sum(axis=1) + np.sin(3 * X).sum(axis=1)
    best_likelihood = -np.inf
    for random_start in range(1, 5):
        gp = GaussianProcess(regr="constant", corr="squared_exponential",
                             theta0=[1e-0] * n_features,
                             thetaL=[1e-4] * n_features,
                             thetaU=[1e+1] * n_features,
                             random_start=random_start, random_state=0,
                             verbose=False).fit(X, y)
        rlf = gp.reduced_likelihood_function()[0]
        assert_greater(rlf, best_likelihood - np.finfo(np.float32).eps)
        best_likelihood = rlf


def test_mse_solving():
    # test the MSE estimate to be sane.
    # non-regression test for ignoring off-diagonals of feature covariance,
    # testing with nugget that renders covariance useless, only
    # using the mean function, with low effective rank of data
    gp = GaussianProcess(corr='absolute_exponential', theta0=1e-4,
                         thetaL=1e-12, thetaU=1e-2, nugget=1e-2,
                         optimizer='Welch', regr="linear", random_state=0)

    X, y = make_regression(n_informative=3, n_features=60, noise=50,
                           random_state=0, effective_rank=1)

    gp.fit(X, y)
    assert_greater(1000, gp.predict(X, eval_MSE=True)[1].mean())
