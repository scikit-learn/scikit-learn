# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import SkipTest
from sklearn.linear_model.bayes import BayesianRidge, ARDRegression
from sklearn.linear_model import Ridge
from sklearn import datasets


def test_bayesian_on_diabetes():
    # Test BayesianRidge on diabetes
    raise SkipTest("XFailed Test")
    diabetes = datasets.load_diabetes()
    X, y = diabetes.data, diabetes.target

    clf = BayesianRidge(compute_score=True)

    # Test with more samples than features
    clf.fit(X, y)
    # Test that scores are increasing at each iteration
    assert_array_equal(np.diff(clf.scores_) > 0, True)

    # Test with more features than samples
    X = X[:5, :]
    y = y[:5]
    clf.fit(X, y)
    # Test that scores are increasing at each iteration
    assert_array_equal(np.diff(clf.scores_) > 0, True)


def test_bayesian_ridge_parameter():
    # Test correctness of lambda_ and alpha_ parameters (GitHub issue #8224)
    X = np.array([[1, 1], [3, 4], [5, 7], [4, 1], [2, 6], [3, 10], [3, 2]])
    y = np.array([1, 2, 3, 2, 0, 4, 5]).T

    # A Ridge regression model using an alpha value equal to the ratio of
    # lambda_ and alpha_ from the Bayesian Ridge model must be identical
    br_model = BayesianRidge(compute_score=True).fit(X, y)
    rr_model = Ridge(alpha=br_model.lambda_ / br_model.alpha_).fit(X, y)
    assert_array_almost_equal(rr_model.coef_, br_model.coef_)
    assert_almost_equal(rr_model.intercept_, br_model.intercept_)


def test_toy_bayesian_ridge_object():
    # Test BayesianRidge on toy
    X = np.array([[1], [2], [6], [8], [10]])
    Y = np.array([1, 2, 6, 8, 10])
    clf = BayesianRidge(compute_score=True)
    clf.fit(X, Y)

    # Check that the model could approximately learn the identity function
    test = [[1], [3], [4]]
    assert_array_almost_equal(clf.predict(test), [1, 3, 4], 2)


def test_toy_ard_object():
    # Test BayesianRegression ARD classifier
    X = np.array([[1], [2], [3]])
    Y = np.array([1, 2, 3])
    clf = ARDRegression(compute_score=True)
    clf.fit(X, Y)

    # Check that the model could approximately learn the identity function
    test = [[1], [3], [4]]
    assert_array_almost_equal(clf.predict(test), [1, 3, 4], 2)


def test_return_std():
    # Test return_std option for both Bayesian regressors
    def f(X):
        return np.dot(X, w) + b

    def f_noise(X, noise_mult):
        return f(X) + np.random.randn(X.shape[0]) * noise_mult

    d = 5
    n_train = 50
    n_test = 10

    w = np.array([1.0, 0.0, 1.0, -1.0, 0.0])
    b = 1.0

    X = np.random.random((n_train, d))
    X_test = np.random.random((n_test, d))

    for decimal, noise_mult in enumerate([1, 0.1, 0.01]):
        y = f_noise(X, noise_mult)

        m1 = BayesianRidge()
        m1.fit(X, y)
        y_mean1, y_std1 = m1.predict(X_test, return_std=True)
        assert_array_almost_equal(y_std1, noise_mult, decimal=decimal)

        m2 = ARDRegression()
        m2.fit(X, y)
        y_mean2, y_std2 = m2.predict(X_test, return_std=True)
        assert_array_almost_equal(y_std2, noise_mult, decimal=decimal)
