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
    # Test correctness of lambda_ and alpha_ parameters (Github issue #8224)
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


def test_dtype_match():
    # Test that np.float32 input data is not cast to np.float64 when possible
    X = np.array([[1, 1], [3, 4], [5, 7], [4, 1], [2, 6], [3, 10], [3, 2]])
    y = np.array([1, 2, 3, 2, 0, 4, 5]).T
    X_32 = np.array(X).astype(np.float32)
    y_32 = np.array(y).astype(np.float32)
    X_64 = np.array(X).astype(np.float64)
    y_64 = np.array(y).astype(np.float64)

    # check type consistency
    for X, dtype in [(X_32, np.float32), (X_64, np.float64)]:
        for y in [y_32, y_64]:
            for model in [BayesianRidge(), ARDRegression()]:
                model.fit(X, y)
                y_mean, y_std = model.predict(X, return_std=True)
                assert X.dtype == dtype
                assert model.coef_.dtype == X.dtype
                assert y_mean.dtype == X.dtype
                assert y_std.dtype == X.dtype

    # check accuracy
    br_32 = BayesianRidge()
    br_32.fit(X_32, y_32)
    br_64 = BayesianRidge()
    br_64.fit(X_64, y_64)
    ard_32 = ARDRegression()
    ard_32.fit(X_32, y_32)
    ard_64 = ARDRegression()
    ard_64.fit(X_64, y_64)

    assert_array_almost_equal(br_32.coef_, br_64.coef_, decimal=5)
    assert_array_almost_equal(ard_32.coef_, ard_64.coef_, decimal=5)
