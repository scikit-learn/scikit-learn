# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from itertools import product
from math import prod

import numpy as np
import pytest
from numpy.polynomial.legendre import legroots
from scipy.stats import beta, expon, norm, uniform

from sklearn.base import BaseEstimator
from sklearn.datasets import make_friedman1
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import (
    ElasticNetCV,
    LassoCV,
    LinearRegression,
    OrthogonalMatchingPursuitCV,
)
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.polynomial_chaos import PolynomialChaosExpansion
from sklearn.utils import check_random_state
from sklearn.utils._orthogonal_polynomial import Polynomial


# Test exact coefficients
def test_exact_coefficients():
    random_state = check_random_state(123)
    distribution = uniform(loc=-1, scale=2)
    X = distribution.rvs((27, 2), random_state=random_state)
    y = (
        3.14
        + 1.74 * 0.5 * np.sqrt(5) * (3 * X[:, 0] ** 2 - 1)
        + np.sqrt(3) * X[:, 1] * 0.5 * np.sqrt(5) * (3 * X[:, 0] ** 2 - 1)
    )
    pce = PolynomialChaosExpansion(distribution, degree=3, scale_outputs=False)
    pce.fit(X, y)
    assert np.linalg.norm(y - pce.predict(X)) < 1e-12
    for multiindex, coef in zip(pce.multiindices_, pce.coef_):
        if np.all(multiindex == np.array([0, 0])):
            assert abs(coef - 3.14) < 1e-12
        elif np.all(multiindex == np.array([2, 0])):
            assert abs(coef - 1.74) < 1e-12
        elif np.all(multiindex == np.array([2, 1])):
            assert abs(coef - 1) < 1e-12
        else:
            assert abs(coef) < 1e-12


# Test 1d fitting
def test_1d_fit():
    X = np.atleast_2d([1.0, 3.0, 5.0, 6.0, 7.0, 8.0]).T
    y = (X * np.sin(X)).ravel()
    pce = PolynomialChaosExpansion(degree=5)
    pce.fit(X, y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(y - y_fit) < 1e-12


# Test fitting with distributions
def test_fit_distributions():
    # Generate data
    random_state = check_random_state(17)
    X = uniform().rvs((28, 2), random_state=random_state)
    y = np.prod((3 * X**2 + 1) / 2, axis=1)

    # Only 1 distribution is given
    pce = PolynomialChaosExpansion(uniform(), degree=6)
    pce.fit(X, y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(y - y_fit) < 1e-12

    # Multiple distributions are given
    pce = PolynomialChaosExpansion((uniform(), uniform()), degree=6)
    pce.fit(X, y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(y - y_fit) < 1e-12

    # Unmatched number of distributions throws error
    with pytest.raises(ValueError, match="number of distributions"):
        pce = PolynomialChaosExpansion((uniform(),))
        pce.fit(X, y)

    # Unmatched distribution type throws error
    with pytest.raises(ValueError, match="'dist'"):
        pce = PolynomialChaosExpansion(False)
        pce.fit(X, y)

    # Unmatched distribution types throws error
    with pytest.raises(ValueError, match="frozen distribution"):
        pce = PolynomialChaosExpansion((uniform(), uniform))
        pce.fit(X, y)


# Test fit with estimators
def test_fit_estimators():
    # Generate data
    random_state = check_random_state(17)
    X = uniform().rvs((28, 2), random_state=random_state)
    y = np.prod((3 * X**2 + 1) / 2, axis=1)

    # Passes
    pce = PolynomialChaosExpansion(
        (uniform(), uniform()), degree=12, estimator=LassoCV(fit_intercept=False)
    )
    pce.fit(X, y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(y - y_fit) / np.linalg.norm(y) < 5e-3

    # Check for fit_intercept
    with pytest.raises(ValueError, match="fit_intercept"):
        pce = PolynomialChaosExpansion(
            (uniform(), uniform()),
            degree=12,
            estimator=LassoCV(fit_intercept=True),
        )
        pce.fit(X, y)

    # Unknown estimator type raises error
    with pytest.raises(ValueError, match="fit"):
        pce = PolynomialChaosExpansion(
            (uniform(), uniform()), degree=12, estimator=False
        )
        pce.fit(X, y)


# Test estimators with noise
@pytest.mark.parametrize(
    "estimator",
    [
        LassoCV(fit_intercept=False, alphas=np.logspace(-12, 2, 25)),
        ElasticNetCV(fit_intercept=False, l1_ratio=np.logspace(-2, 0, 25)),
        OrthogonalMatchingPursuitCV(fit_intercept=False),
    ],
)
def test_estimators_with_noise(estimator):
    random_state = check_random_state(123)
    distribution = uniform(loc=-1, scale=2)
    X = distribution.rvs((27, 2), random_state=random_state)
    y = (
        3.14
        + 1.74 * 0.5 * (3 * X[:, 0] ** 2 - 1)
        + X[:, 1] * 0.5 * (3 * X[:, 0] ** 2 - 1)
    )
    y += 0.1 * random_state.randn(len(y))
    pce = PolynomialChaosExpansion(
        distribution, degree=3, scale_outputs=False, estimator=estimator
    )
    pce.fit(X, y)
    assert np.linalg.norm(y - pce.predict(X)) / np.linalg.norm(y) < 0.05


# Test input checking for fit
def test_fit_inputs():
    # Only 1 data point throws an error
    with pytest.raises(ValueError, match="more than 1 sample"):
        pce = PolynomialChaosExpansion()
        pce.fit(np.atleast_2d([1]), [1])


# Test no variation in feature
def test_no_variation():
    x1 = [1.0, 3.0, 5.0, 6.0, 7.0, 8.0]
    x2 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    X = np.vstack([x1, x2]).T
    y1 = (x1 * np.sin(x1)).ravel()
    y2 = (np.cos(x2) + 1).ravel()
    Y = np.vstack([y1, y2]).T
    pce = PolynomialChaosExpansion(degree=5)
    pce.fit(X, Y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(Y - y_fit) < 1e-12


# Test grid search
def test_grid_search():
    random_state = check_random_state(123)
    distribution = uniform(loc=-1, scale=2)
    X = distribution.rvs((27, 2), random_state=random_state)
    y = (
        3.14
        + 1.74 * 0.5 * (3 * X[:, 0] ** 2 - 1)
        + X[:, 1] * 0.5 * (3 * X[:, 0] ** 2 - 1)
    )
    param_grid = [
        {
            "degree": [0, 1, 2, 3, 4],
            "truncation": [
                "full_tensor",
                "total_degree",
                "hyperbolic_cross",
                "Zaremba_cross",
            ],
            "estimator": [
                LinearRegression(fit_intercept=False),
                LassoCV(fit_intercept=False, alphas=np.logspace(-12, 2, 25)),
            ],
        }
    ]
    pceCV = GridSearchCV(
        PolynomialChaosExpansion(distribution),
        param_grid,
        scoring="neg_root_mean_squared_error",
    )
    pceCV.fit(X, y)
    assert pceCV.best_params_["degree"] < 4


# Test grid search 1 polynomial
def test_grid_search_1_polynomial():
    random_state = check_random_state(123)
    distribution = uniform(loc=-1, scale=2)
    X = distribution.rvs((100, 5), random_state=random_state)
    y = 1 / 8 * (35 * X[:, 3] ** 4 - 30 * X[:, 3] ** 2 + 3)
    param_grid = [
        {
            "degree": range(2, 6),
        }
    ]
    estimator = LassoCV(
        fit_intercept=False, alphas=np.logspace(-12, 2, 25), max_iter=100000
    )
    pce = PolynomialChaosExpansion(
        distribution, estimator=estimator, scale_outputs=False
    )
    pceCV = GridSearchCV(pce, param_grid, cv=KFold(n_splits=5))
    pceCV.fit(X, y)
    assert pceCV.best_params_["degree"] == 4
    idx = np.argmax(pceCV.best_estimator_.coef_)
    multiindex = pceCV.best_estimator_.multiindices_[idx]
    assert np.all(multiindex == [0, 0, 0, 4, 0])
    assert np.argmax(pceCV.best_estimator_.main_sens()) == 3
    assert np.argmax(pceCV.best_estimator_.total_sens()) == 3


# Test predict inputs
def test_predict_inputs():
    random_state = check_random_state(123)
    distribution = uniform(loc=-1, scale=2)
    X = distribution.rvs((27, 2), random_state=random_state)
    y = (
        3.14
        + 1.74 * 0.5 * (3 * X[:, 0] ** 2 - 1)
        + X[:, 1] * 0.5 * (3 * X[:, 0] ** 2 - 1)
    )
    pce = PolynomialChaosExpansion(distribution, degree=3, scale_outputs=False)
    pce.fit(X, y)
    with pytest.raises(ValueError, match="expecting 2 features"):
        pce.predict(X[:, 0].reshape(-1, 1))


# Test main and total sensitivity indices
def test_main_sensitivity():
    random_state = check_random_state(17)
    distribution = uniform(0, 1)
    pce = PolynomialChaosExpansion(distribution, degree=4)
    X = distribution.rvs((60, 1), random_state=random_state)
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce.fit(X, y)
    mains = pce.main_sens()
    for j, main in enumerate(mains):
        assert abs(main - pce.joint_sens(j)) < 1e-12


# Test main sensitivity on Friedman 1 problem
def test_friedman1():
    X, y = make_friedman1(n_samples=100, n_features=10, random_state=0)
    estimator = LassoCV(
        fit_intercept=False, alphas=np.logspace(-12, 2, 25), max_iter=100000
    )
    pce = PolynomialChaosExpansion(
        uniform(-1, 2), degree=5, estimator=estimator, scale_outputs=False
    )
    pce.fit(X, y)
    assert np.sum(pce.main_sens() > 1e-14) == 5  # 5 important features


# Check mean and var of distributions
@pytest.mark.parametrize(
    "distribution, mean, var",
    [
        (uniform(loc=-3, scale=0.1), -2.95, 0.0008333333333333333),
        (norm(loc=-3, scale=0.1), -3, 0.01),
        (expon(loc=-3, scale=0.1), -2.9, 0.01),
        (beta(2, 2, loc=-3, scale=0.1), -2.95, 0.0005),
    ],
)
def test_mean_and_var(distribution, mean, var):
    random_state = check_random_state(17)
    X = distribution.rvs((2, 1), random_state=random_state)
    y = X.ravel()
    pce = PolynomialChaosExpansion(distribution, degree=1, scale_outputs=False)
    pce.fit(X, y)
    assert abs(pce.mean() - mean) < 1e-12
    assert abs(pce.var() - var) < 1e-12


# Check model with zero mean
def test_zero_mean():
    random_state = check_random_state(17)
    distribution = uniform(0, 1)
    pce = PolynomialChaosExpansion(distribution, degree=2, scale_outputs=False)
    X = distribution.rvs((60, 1), random_state=random_state)
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce.fit(X, y)
    pce.coef_ = pce.coef_[1:]
    pce.multiindices_ = pce.multiindices_[1:]
    pce.norms_ = pce.norms_[1:]
    assert pce.mean() == 0


# Check for named features
def test_pandas():
    pd = pytest.importorskip("pandas")
    distribution = uniform()
    random_state = check_random_state(17)
    X = distribution.rvs((116, 3), random_state=random_state)
    X = pd.DataFrame(data=X, columns=[f"feature{j}" for j in range(3)])
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce = PolynomialChaosExpansion(distribution, degree=6)
    pce.fit(X, y)

    assert abs(pce.joint_sens("feature0") - 25 / 91) < 1e-12
    assert abs(pce.joint_sens("feature0", "feature1") - 5 / 91) < 1e-12
    assert abs(pce.joint_sens("feature0", "feature1", "feature2") - 1 / 91) < 1e-12


# Verify input checking for joint_sens
def test_joint_sens_inputs():
    pd = pytest.importorskip("pandas")
    distribution = uniform()
    random_state = check_random_state(17)
    X = distribution.rvs((116, 3), random_state=random_state)
    X_named = pd.DataFrame(data=X, columns=[f"feature{j}" for j in range(3)])
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce = PolynomialChaosExpansion(distribution, degree=6)
    pce.fit(X_named, y)

    with pytest.raises(ValueError, match="all string"):
        pce.joint_sens("feature0", 1)

    with pytest.raises(ValueError, match="'wololo'"):
        pce.joint_sens("feature0", "wololo")

    with pytest.raises(ValueError, match="unique"):
        pce.joint_sens("feature0", "feature0")

    with pytest.raises(ValueError, match="3 features"):
        pce.joint_sens(2, 1, 3, 0)

    with pytest.raises(ValueError, match="feature names"):
        pce.fit(X, y)
        pce.joint_sens("feature0")


# Test multioutput fit
def test_multi_output_fit():
    X = np.atleast_2d([1.0, 3.0, 5.0, 6.0, 7.0, 8.0]).T
    y1 = (X * np.sin(X)).ravel()
    y2 = (X * np.cos(X)).ravel()
    Y = np.vstack([y1, y2]).T
    pce = PolynomialChaosExpansion(degree=5)
    pce.fit(X, Y)
    y_fit = pce.predict(X)
    assert np.linalg.norm(Y - y_fit) < 1e-12


# Test multioutput statistics
def test_multi_output_statistics():
    X = np.atleast_2d([1.0, 3.0, 5.0, 6.0, 7.0, 8.0]).T
    y1 = (X * np.sin(X)).ravel()
    y2 = (X * np.cos(X)).ravel()
    Y = np.vstack([y1, y2]).T
    pce = PolynomialChaosExpansion(degree=2)

    for j, y in enumerate([y1, y2]):
        assert np.sum(pce.fit(X, Y).mean()[j] - pce.fit(X, y).mean()) < 1e-12
        assert np.sum(pce.fit(X, Y).var()[j] - pce.fit(X, y).var()) < 1e-12
        assert np.sum(pce.fit(X, Y).main_sens()[j] - pce.fit(X, y).main_sens()) < 1e-12
        assert (
            np.sum(pce.fit(X, Y).total_sens()[j] - pce.fit(X, y).total_sens()) < 1e-12
        )


# Verify input checking for estimator with multioutput
def test_estimator_multioutput():
    X = np.atleast_2d([1.0, 3.0, 5.0, 6.0, 7.0, 8.0]).T
    y1 = (X * np.sin(X)).ravel()
    y2 = (X * np.cos(X)).ravel()
    Y = np.vstack([y1, y2]).T

    with pytest.raises(ValueError, match="fit_intercept=False"):
        pce = PolynomialChaosExpansion(estimator=LassoCV())
        pce.fit(X, y1)

    with pytest.raises(ValueError, match="fit_intercept=False"):
        pce = PolynomialChaosExpansion(estimator=MultiOutputRegressor(LassoCV()))
        pce.fit(X, Y)

    class DummyEstimator(BaseEstimator):
        def fit(self, X, y):
            self.coef_ = 0

    with pytest.warns(UserWarning, match="fit_intercept=False"):
        pce = PolynomialChaosExpansion(estimator=DummyEstimator())
        pce.fit(X, Y)

    class DumbDummyEstimator(BaseEstimator):
        def fit(self, X, y):
            pass

    with pytest.warns(UserWarning, match="fit_intercept=False"):
        with pytest.raises(ValueError, match="'coef_'"):
            pce = PolynomialChaosExpansion(estimator=DumbDummyEstimator())
            pce.fit(X, Y)


# Test feature selector
def test_feature_selector_valid():
    distribution = uniform(loc=-1, scale=2)
    random_state = check_random_state(17)
    X = distribution.rvs((100, 2), random_state=random_state)
    y = X[:, 0] * X[:, 1]
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    pce = PolynomialChaosExpansion(
        distribution=[distribution] * 2,
        degree=3,
        estimator=LinearRegression(fit_intercept=False),
        feature_selector=SelectFromModel(LassoCV(fit_intercept=False)),
    )
    pce.fit(X_train, y_train)

    assert pce.score(X_test, y_test) > 0.999
    assert np.linalg.norm(pce.main_sens()) < 1e-12
    assert np.linalg.norm(pce.total_sens() - np.array([1, 1])) < 1e-12
    assert np.abs(pce.joint_sens(0, 1) - 1) < 1e-12


# Verify input checking for feature selector
def test_feature_selector():
    distribution = uniform(loc=-1, scale=2)
    random_state = check_random_state(17)
    X = distribution.rvs((100, 2), random_state=random_state)
    y = X[:, 0] * X[:, 1]

    class DummyFeatureSelector(BaseEstimator):
        def fit(self, X, y):
            pass

    with pytest.raises(ValueError, match="fit and transform"):
        pce = PolynomialChaosExpansion(feature_selector=DummyFeatureSelector())
        pce.fit(X, y)


###############################################################################
# In the following unit tests we compare the results obtained with our
# implementation to the results reported in the review paper "Global
# sensitivity analysis using polynomial chaos expansions" by B. Sudret (2008).
###############################################################################


# Utility to generate sample points from Sudret (2008)
def get_samples(degree, dimension, tol=12):
    c = [0] * (degree + 2)
    c[-1] = 1
    roots = np.round(legroots(c), tol)  # Round roots themselves
    samples = list(product(*[roots for _ in range(dimension)]))
    samples.sort(
        key=lambda x: (round(np.linalg.norm(x), tol),) + tuple(np.round(x, tol))
    )  # Sort first by norm, then lexicographically
    return np.vstack(samples)


# Example 1 from Sudret (2008) - exact values
@pytest.mark.parametrize(
    "indices, SU",
    [
        ((0,), 25 / 91),
        ((1,), 25 / 91),
        ((2,), 25 / 91),
        ((0, 1), 5 / 91),
        ((0, 2), 5 / 91),
        ((1, 2), 5 / 91),
        ((0, 1, 2), 1 / 91),
    ],
)
def test_polynomial_model_exact(indices, SU):
    pce = PolynomialChaosExpansion(uniform(0, 1), degree=6)
    X = (get_samples(6, 3)[:116, :] + 1) / 2
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce.fit(X, y)

    assert abs(pce.joint_sens(*indices) - SU) < 1e-12


# Example 1 from Sudret (2008) - blind solution
# NOTE: exact values of sensitivity indices reported in the paper do not agree
# for orders 3 and 4. The values obtained with the sklearn implementation do
# satisfy the reported relative errors, so we compare to those instead. A
# possible reason could be that the sample point selection algorithm described
# in the paper does not lead to unique sample points (a lot of points have the
# same 'norm').
@pytest.mark.parametrize(
    "degree, N, relative_error",
    [
        (3, 29, (0.05, 0.12, 0.5)),
        (4, 44, (0.01, 0.04, 0.02)),
        (5, 77, (0.001, 0.001, 0.001)),
    ],
)
def test_polynomial_model(degree, N, relative_error):
    pce = PolynomialChaosExpansion(uniform(0, 1), degree=degree)
    X = (get_samples(degree, 3)[:N, :] + 1) / 2
    y = np.prod((3 * X**2 + 1) / 2, axis=1)
    pce.fit(X, y)

    for idcs in ((0,), (1,), (2,)):
        assert (pce.joint_sens(*idcs) - 25 / 91) / (25 / 91) < relative_error[0]

    for idcs in ((0, 1), (1, 2), (0, 2)):
        assert (pce.joint_sens(*idcs) - 5 / 91) / (5 / 91) < relative_error[1]

    assert (pce.joint_sens(0, 1, 2) - 1 / 91) / (1 / 91) < relative_error[2]


# Example 2 from Sudret (2008)
# NOTE: values from Table 2 in the paper
@pytest.mark.parametrize(
    "degree, N, main, total",
    [
        (
            3,
            29,
            (0.3941, 0.1330, 0.0000, 0.0000, 0.4729, 0.0000, 0.0000),
            (0.8670, 0.1330, 0.4729),
        ),
        (
            5,
            77,
            (0.3550, 0.1846, 0.0000, 0.0000, 0.4603, 0.0000, 0.0000),
            (0.8154, 0.1846, 0.4603),
        ),
        (
            7,
            157,
            (0.3114, 0.4434, 0.0000, 0.0000, 0.2452, 0.0000, 0.0000),
            (0.5566, 0.4434, 0.2452),
        ),
        (
            9,
            291,
            (0.3146, 0.4396, 0.0000, 0.0000, 0.2459, 0.0000, 0.0000),
            (0.5604, 0.4396, 0.2459),
        ),
    ],
)
def test_ishigami(degree, N, main, total):
    a = 7
    b = 0.1
    dimension = 3
    pce = PolynomialChaosExpansion(uniform(-np.pi, 2 * np.pi), degree=degree)
    X = np.pi * get_samples(degree, dimension)[:N, :]
    y = np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2 + b * X[:, 2] ** 4 * np.sin(X[:, 0])
    pce.fit(X, y)

    for j, idcs in enumerate(((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2))):
        assert abs(pce.joint_sens(*idcs) - main[j]) < 1e-4

    for j, tot in enumerate(total):
        assert abs(pce.total_sens()[0, j] - tot) < 1e-4


def test_ishigami_analytically():
    pass


# Example 3 from Sudret (2008)
# NOTE: values from Table 3 in the paper
def test_sobol_order_2():
    a = np.array([1, 2, 5, 10, 20, 50, 100, 500])
    dimension = len(a)
    degree = 2
    N = 72
    pce = PolynomialChaosExpansion(uniform(), degree=degree)
    X = (get_samples(degree, dimension)[:N, :] + 1) / 2
    y = prod((abs(4 * X_j - 2) + a_j) / (1 + a_j) for a_j, X_j in zip(a, X.T))
    pce.fit(X, y)

    exact = [0.5986, 0.3045, 0.0426, 0.0091, 0.0041, 0.0034, 0.0001, 0.0002]
    for j, sens in enumerate(exact):
        assert abs(pce.joint_sens(j) - sens) < 1e-4

    exact = [0.0140, 0.0004, 0.0005, 0.0015, 0.0019, 0.0000, 0.0000]
    for j, sens in enumerate(exact):
        assert abs(pce.joint_sens(0, j + 1) - sens) < 1e-4

    exact = [0.6170, 0.3229, 0.0473, 0.0140, 0.0122, 0.0129, 0.0055, 0.0057]
    for j, sens in enumerate(exact):
        assert abs(pce.total_sens()[0, j] - sens) < 1e-4


# Example 3 from Sudret (2008)
# NOTE: values from Table 4 in the paper
# I can exactly reproduce the values from the table for order 5, but not for
# Order 3 and order 7. See remark for Example 1 above.
@pytest.mark.parametrize(
    "degree, N, main, total, tol",
    [
        (
            3,
            73,
            (
                0.6644,
                0.2616,
                0.0581,
                0.0164,
                0.0000,
                0.0000,
                0.0000,
                0.0000,
                0.0000,
                0.0000,
            ),
            (0.6644, 0.2611, 0.0581, 0.0164),
            6e-4,
        ),
        (
            5,
            233,
            (
                0.5994,
                0.2748,
                0.0676,
                0.0198,
                0.0283,
                0.0058,
                0.0016,
                0.0021,
                0.0006,
                0.0001,
            ),
            (0.6350, 0.3057, 0.0756, 0.0220),
            1e-4,
        ),
        (
            7,
            990,  # 533
            (
                0.5999,
                0.2677,
                0.0729,
                0.0226,
                0.0191,
                0.0069,
                0.0023,
                0.0043,
                0.0013,
                0.0004,
            ),
            (0.6308, 0.2950, 0.0866, 0.0272),
            1e-1,
        ),
    ],
)
def test_sobol(degree, N, main, total, tol):
    a = np.array([1, 2, 5, 10, 20, 50, 100, 500])
    dimension = 4
    pce = PolynomialChaosExpansion(uniform(), degree=degree)
    X = (get_samples(degree, dimension)[:N, :] + 1) / 2
    X2 = np.hstack([X, 1 / 2 * np.ones_like(X)])
    y = prod((abs(4 * X_j - 2) + a_j) / (1 + a_j) for a_j, X_j in zip(a, X2.T))
    pce.fit(X, y)

    joint = (
        (0,),
        (1,),
        (2,),
        (3,),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    )
    for j, idcs in enumerate(joint):
        assert abs(pce.joint_sens(*idcs) - main[j]) < tol

    for j, tot in enumerate(total):
        assert abs(pce.total_sens()[0, j] - tot) < tol


# Test Ishigami example with Fourier basis
# NOTE: The Ishigami function is defined as
#
#    sin(x_1) + a * sin(x_2)^2 + b * x_3^4 * sin(x_1)
#
# where x_1, x_2, x_3 ~ U(-pi, pi)
# Using a Fourier basis for x_1 and x_2, and using a Legendre polynomial basis
# for x_2, this function can be fitted exactly using only 5 polynomial terms.


# Define Fourier basis
class Fourier(Polynomial):
    """A class representing Fourier basis functions (cos/sin)."""

    def _vandermonde(self, points, degree):
        points = np.asarray(points)
        n = len(points)
        V = np.zeros((n, degree + 1))

        V[:, 0] = 1.0  # Constant
        for j in range(1, degree + 1):
            k = (j + 1) // 2  # Frequency = ceil(j/2)
            if j % 2 == 1:  # Odd index → sin
                V[:, j] = np.sqrt(2) * np.sin(k * points)
            else:  # Even index → cos
                V[:, j] = np.sqrt(2) * np.cos(k * points)
        return V

    @staticmethod
    def _distribution():
        return "fourier"

    def scale_features_from_distribution(self, X, distribution):
        # Map uniform(0,1) -> [-pi, pi]
        X = super().scale_features_from_distribution(X, distribution)
        return 2 * np.pi * (X - 0.5)

    def _norm_squared(self, degree):
        # This is an orthonormal basis with respect to the uniform distribution
        # on (-pi, pi)
        return 1


# Define corresponding distribution
class FourierDistribution:
    def __init__(self, loc=-np.pi, scale=2 * np.pi):
        self._frozen = uniform(loc=loc, scale=scale)
        self.args = self._frozen.args
        self.kwds = self._frozen.kwds
        # Clone dist and rename
        self.dist = type(
            "dist",
            (),
            {"name": "fourier", "_parse_args": self._frozen.dist._parse_args},
        )()

    def __getattr__(self, name):
        return getattr(self._frozen, name)


# Test Ishigami with Fourier basis
def test_ishigami_fourier():
    # First, let's define the model
    a = 7
    b = 0.1

    def ishigami(X):
        return (
            np.sin(X[:, 0])
            + a * np.sin(X[:, 1]) ** 2
            + b * X[:, 2] ** 4 * np.sin(X[:, 0])
        )

    # Generate training data
    distribution = uniform(loc=-np.pi, scale=2 * np.pi)
    random_state = check_random_state(17)
    X = distribution.rvs((1000, 3), random_state=random_state)
    y = ishigami(X)

    # Define basis, Fourier in x2, Legendre in x3
    dist = (
        FourierDistribution(),  # Fourier in x1
        FourierDistribution(),  # Fourier in x2
        distribution,  # Legendre in x3
    )

    # Fit PCE
    pce = PolynomialChaosExpansion(
        distribution=dist,
        degree=5,
        scale_outputs=False,
    )
    pce.fit(X, y)

    # Check score
    X_test = distribution.rvs((1024, 3))
    y_test = ishigami(X_test)
    assert np.abs(pce.score(X_test, y_test) - 1) < 1e-10

    # Check coefficients
    a = 7.0
    b = 0.1
    c0 = 1.0 / 5.0
    c2 = 4.0 / 7.0
    c4 = 8.0 / 35.0
    exact = [
        a / 2.0,  # (0, 0, 0)
        (1.0 + b * np.pi**4 * c0) / np.sqrt(2),  # (1, 0, 0)
        -a / (2.0 * np.sqrt(2)),  # (0, 4, 0)
        b * np.pi**4 * c2 / (np.sqrt(2) * np.sqrt(5)),  # (1, 0, 2)
        b * np.pi**4 * c4 / (3.0 * np.sqrt(2)),  # (1, 0, 4)
    ]
    idcs = np.array([0, 1, 18, 37, 53])

    # Check expected indices
    atol = 1e-14
    assert np.allclose(pce.coef_[idcs], exact, atol=atol)

    # Check all *other* indices are small
    mask = np.ones_like(pce.coef_, dtype=bool)
    mask[idcs] = False
    assert np.all(np.abs(pce.coef_[mask]) < atol)
