# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

import numpy as np

from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from ..base import LinearRegression
from ..bayes import BayesianRidge, Ridge, ARDRegression



def test_ridge():
    """
    TODO: for this test to be robust, we should use a dataset instead
    of np.random.
    """
    alpha = 1.0

    # With more samples than features
    n_samples, n_features = 6, 5
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
    assert ridge.score (X, y) > 0.5

    # With more features than samples
    n_samples, n_features = 5, 10
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
    assert ridge.score (X, y) > .9

def test_ridge_vs_lstsq():
    """
    On alpha=0., Ridge and OLS yield the same solution.
    """

    # we need more samples than features
    n_samples, n_features = 5, 4
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    ridge = Ridge(alpha=0.)
    ols = LinearRegression()
    
    ridge.fit(X, y)
    ols.fit (X, y)
    assert np.linalg.norm (ridge.coef_ - ols.coef_) < 1e-10

    ridge.fit(X, y, fit_intercept=False)
    ols.fit (X, y, fit_intercept=False)
    assert np.linalg.norm (ridge.coef_ - ols.coef_) < 1e-10




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
    X = np.array([[1], [2], [6], [8], [10]])
    Y = np.array([1, 2, 6, 8, 10])
    clf = BayesianRidge()
    clf.fit(X, Y)
    Test = [[1], [3], [4]]
    assert(np.abs(clf.predict(Test)-[1, 3, 4]).sum()<1.e-1) # identity


def test_toy_ard_object():
    """
    Test BayesianRegression ARD classifier
    """
    X = np.array([[1], [2], [6], [8], [10]])
    Y = np.array([1, 2, 6, 8, 10])
    clf = ARDRegression()
    clf.fit(X, Y)
    Test = [[1], [3], [4]]
    assert(np.abs(clf.predict(Test)-[1, 3, 4]).sum()<1.e-1) # identity
