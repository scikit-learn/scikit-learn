"""
Testing for Gaussian Process module (sklearn.gaussian_process)
"""

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
# Licence: BSD 3 clause
from macpath import norm_error

from nose.tools import raises
from nose.tools import assert_true

import numpy as np
import scipy.stats as ss

from sklearn.gaussian_process import GaussianProcess
from sklearn.gaussian_process import regression_models as regression
from sklearn.gaussian_process import correlation_models as correlation


f = lambda x: x * np.sin(x)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
X2 = np.atleast_2d([2., 4., 5.5, 6.5, 7.5]).T
y = f(X).ravel()


def test_1d(regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    """
    MLE estimation of a one-dimensional Gaussian Process model.
    Check random start optimization.

    Test the interpolating property.
    """
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
                         random_start=random_start, verbose=False).fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)
    y2_pred, MSE2 = gp.predict(X2, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.)
                and np.allclose(MSE2, 0., atol=10))


def check2d(X, regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    """
    MLE estimation of a two-dimensional Gaussian Process model accounting for
    anisotropy. Check random start optimization.

    Test the interpolating property.
    """
    b, kappa, e = 5., .5, .1
    g = lambda x: b - x[:, 1] - kappa * (x[:, 0] - e) ** 2.

    y = g(X).ravel()
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=[1e-2] * 2, thetaL=[1e-4] * 2,
                         thetaU=[1e-1] * 2,
                         random_start=random_start, verbose=False)
    gp.fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.))


def check2d_2d(X, regr=regression.constant, corr=correlation.squared_exponential,
               random_start=10, beta0=None, theta0=[1e-2] * 2, thetaL=[1e-4] * 2, thetaU=[1e-1] * 2):
    """
    MLE estimation of a two-dimensional Gaussian Process model accounting for
    anisotropy. Check random start optimization.

    Test the GP interpolation for 2D output
    """
    b, kappa, e = 5., .5, .1
    g = lambda x: b - x[:, 1] - kappa * (x[:, 0] - e) ** 2.
    f = lambda x: np.vstack((g(x), g(x))).T

    y = f(X)
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=theta0, thetaL=thetaL,
                         thetaU=thetaU,
                         random_start=random_start, verbose=False)
    gp.fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert_true(np.allclose(y_pred, y) and np.allclose(MSE, 0.))


@raises(ValueError)
def test_wrong_number_of_outputs():
    gp = GaussianProcess()
    gp.fit([[1, 2, 3], [4, 5, 6]], [1, 2, 3])


def test_more_builtin_correlation_models(random_start=1):
    """
    Repeat test_1d and test_2d for several built-in correlation
    models specified as strings.
    """
    all_corr = ['absolute_exponential', 'squared_exponential', 'cubic',
                'linear']

    X = np.array([[-4.61611719, -6.00099547],
                  [4.10469096, 5.32782448],
                  [0.00000000, -0.50000000],
                  [-6.17289014, -4.6984743],
                  [1.3109306, -6.93271427],
                  [-5.03823144, 3.10584743],
                  [-2.87600388, 6.74310541],
                  [5.21301203, 4.26386883]])

    for corr in all_corr:
        test_1d(regr='constant', corr=corr, random_start=random_start)
        check2d(X, regr='constant', corr=corr, random_start=random_start)
        check2d_2d(X, regr='constant', corr=corr, random_start=random_start)


def test_ordinary_kriging():
    """
    Repeat test_1d and test_2d with given regression weights (beta0) for
    different regression models (Ordinary Kriging).
    """

    X = np.array([[-4.61611719, -6.00099547],
                  [4.10469096, 5.32782448],
                  [0.00000000, -0.50000000],
                  [-6.17289014, -4.6984743],
                  [1.3109306, -6.93271427],
                  [-5.03823144, 3.10584743],
                  [-2.87600388, 6.74310541],
                  [5.21301203, 4.26386883]])

    test_1d(regr='linear', beta0=[0., 0.5])
    test_1d(regr='quadratic', beta0=[0., 0.5, 0.5])
    check2d(X, regr='linear', beta0=[0., 0.5, 0.5])
    check2d(X, regr='quadratic', beta0=[0., 0.5, 0.5, 0.5, 0.5, 0.5])
    check2d_2d(X, regr='linear', beta0=[0., 0.5, 0.5])
    check2d_2d(X, regr='quadratic', beta0=[0., 0.5, 0.5, 0.5, 0.5, 0.5])


def test_no_normalize():
    gp = GaussianProcess(normalize=False).fit(X, y)
    y_pred = gp.predict(X)
    assert_true(np.allclose(y_pred, y))


@raises(Exception)
def test_no_multiple_feature():
    '''
    check that multiple features are not allowed for non-noisy data
    '''
    X = np.array([[-4.61611719, -6.00099547],
                  [4.10469096, 5.32782448],
                  [0.00000000, -0.50000000],
                  [-6.17289014, -4.6984743],
                  [-6.17289014, -4.6984743],
                  [1.3109306, -6.93271427],
                  [-5.03823144, 3.10584743],
                  [-2.87600388, 6.74310541],
                  [5.21301203, 4.26386883]])

    check2d(X)


def test_1d_noisy(regr=regression.constant, corr=correlation.absolute_exponential,
                  random_start=10, beta0=None):
    """
    MLE estimation of a one-dimensional Gaussian Process model.
    Check random start optimization with noisy / duplicate inputs.

    Test the interpolating property.
    """

    X = np.atleast_2d([1., 3., 5., 6., 7., 8., 9., 10.] * 2).T
    x = np.atleast_2d(np.linspace(0, 10, 50)).T

    y = f(X).ravel() + np.random.normal(0, 0.1, len(X))

    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
                         random_start=random_start, verbose=False, nugget=0.01).fit(X, y)
    y_pred, MSE = gp.predict(x, eval_MSE=True)

    y = f(x).ravel()
    assert_true((np.abs(y_pred - y) <= np.abs(
        ss.norm.ppf(0.025, y_pred, np.sqrt(MSE)))  ).all())  #check that true value is within 95% conf. int.
