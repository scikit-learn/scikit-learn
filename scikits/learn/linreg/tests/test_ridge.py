# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id: test_cd.py 450 2010-03-03 14:21:06Z twigster $

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_raises

from ..regression import LinearRegression, bayesian_regression_ridge, Ridge, BayesianRidge



def test_toy_noprior():
    """
    Test BayesianRegression with no prior classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = LinearRegression()
    clf.fit(X, Y)
    T = [[1], [2], [3], [4]]
    assert_array_almost_equal(clf.predict(T), [1, 2, 3, 4]) # identity


def test_bayesian_ridge():
    """
    Test Ridge regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_regression_ridge(X, Y)
    assert np.abs(1-w)<1.e-3

    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = BayesianRidge()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_almost_equal(clf.predict(Test), [1, 2, 3, 4]) # identity

def test_ridge():
    alpha = 1.0

    # With more samples than features
    nsamples, nfeatures = 10, 5
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)

    # With more features than samples
    nsamples, nfeatures = 5, 10
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
