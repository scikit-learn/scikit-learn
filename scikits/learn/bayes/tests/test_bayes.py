import numpy as np
from scikits.learn.bayes.regression import bayesian_regression_noprior, \
                                      bayesian_regression_ridge, \
                                      bayesian_regression_ard, \
                                      RidgeRegression, \
                                      BayesianRegression, \
                                      ARDRegression
from numpy.testing import assert_array_almost_equal

def test_toy_noprior_regression():
    """
    Test no prior regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, beta, log_likelihood = bayesian_regression_noprior(X, Y)
    assert_array_almost_equal(w, [1])


def test_toy_ridge_regression():
    """
    Test Ridge regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_regression_ridge(X, Y)
    assert_array_almost_equal(w, [1])


def test_toy_ard_regression():
    """
    Test ARD regression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w, alpha, beta, sigma, log_likelihood = bayesian_regression_ard(X, Y)
    assert(np.abs(1-w)<1.e-3)


def test_toy_noprior_object():
    """
    Test BayesianRegression with no prior classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = BayesianRegression()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_almost_equal(clf.predict(Test), [1, 2, 3, 4]) # identity


def test_toy_ridge_object():
    """
    Test BayesianRegression ridge classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = RidgeRegression()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_almost_equal(clf.predict(Test), [1, 2, 3, 4]) # identity


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

