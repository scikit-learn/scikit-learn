# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

import numpy as np

from numpy.testing import assert_raises, assert_array_equal, \
     assert_array_almost_equal

from scikits.learn.glm import *

def test_LinearRegression():
    """
    Test LinearRegression on a simple dataset.
    """
    # a simple dataset
    X = [[1], [2]]
    Y = [1, 2]

    clf = LinearRegression()
    clf.fit(X, Y)

    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [1, 2])

    # test it also for degenerate input
    X = [[1]]
    Y = [0]

    clf = LinearRegression()
    clf.fit(X, Y)
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [0])

def test_bayesian_ridge():
    """
    Test Ridge regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_ridge_regression(X, Y)
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
    n_samples, n_features = 10, 5
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)

    # With more features than samples
    n_samples, n_features = 5, 10
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)



def test_toy_ridge_regression():
    """
    Test Ridge regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_ridge_regression(X, Y)
    assert(np.abs(1-w)<1.e-3)


def test_toy_ard_regression():
    """
    Test ARD regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_regression_ard(X, Y)
    assert(np.abs(1-w)<1.e-3)


def test_toy_ridge_object():
    """
    Test BayesianRegression ridge classifier
    TODO: test also n_samples > n_features
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = Ridge(alpha=0.0)
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_equal(clf.predict(Test), [1, 2, 3, 4]) # identity


def test_toy_bayesian_ridge_object():
    """
    Test BayesianRegression ridge classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = BayesianRidge()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_equal(clf.predict(Test), [1, 2, 3, 4]) # identity


def test_toy_ard_object():
    """
    Test BayesianRegression ARD classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = ARDRegression()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert(np.abs(clf.predict(Test)-[1, 2, 3, 4]).sum()<1.e-3) # identity
