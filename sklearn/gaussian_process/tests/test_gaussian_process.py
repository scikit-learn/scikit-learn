"""
Testing for Gaussian Process module (sklearn.gaussian_process)
"""

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
# License: BSD style

from nose.tools import raises
import numpy as np

from .. import GaussianProcess
from .. import regression_models as regression
from .. import correlation_models as correlation


def test_1d(regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    """
    MLE estimation of a one-dimensional Gaussian Process model.
    Check random start optimization.

    Test the interpolating property.
    """
    f = lambda x: x * np.sin(x)
    X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
    y = f(X).ravel()
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=1e-2, thetaL=1e-4, thetaU=1e-1,
                         random_start=random_start, verbose=False).fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert np.allclose(y_pred, y) and np.allclose(MSE, 0.)


def test_2d(regr=regression.constant, corr=correlation.squared_exponential,
            random_start=10, beta0=None):
    """
    MLE estimation of a two-dimensional Gaussian Process model accounting for
    anisotropy. Check random start optimization.

    Test the interpolating property.
    """
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
    gp = GaussianProcess(regr=regr, corr=corr, beta0=beta0,
                         theta0=[1e-2] * 2, thetaL=[1e-4] * 2,
                         thetaU=[1e-1] * 2,
                         random_start=random_start, verbose=False)
    gp.fit(X, y)
    y_pred, MSE = gp.predict(X, eval_MSE=True)

    assert np.allclose(y_pred, y) and np.allclose(MSE, 0.)


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

    for corr in all_corr:
        test_1d(regr='constant', corr=corr, random_start=random_start)
        test_2d(regr='constant', corr=corr, random_start=random_start)


def test_ordinary_kriging():
    """
    Repeat test_1d and test_2d with given regression weights (beta0) for
    different regression models (Ordinary Kriging).
    """
    test_1d(regr='linear', beta0=[0., 0.5])
    test_1d(regr='quadratic', beta0=[0., 0.5, 0.5])
    test_2d(regr='linear', beta0=[0., 0.5, 0.5])
    test_2d(regr='quadratic', beta0=[0., 0.5, 0.5, 0.5, 0.5, 0.5])
