"""
Testing for Theil-Sen module (sklearn.linear_model.theilsen)
"""

# Author: Florian Wilhelm <florian.wilhelm@gmail.com>
# Licence: BSD 3 clause

from __future__ import division, print_function, absolute_import

import logging
from os import devnull

import numpy as np
import numpy.testing as nptest
from numpy.linalg import norm
from scipy.optimize import minimize
from nose.tools import raises
from sklearn.linear_model import LinearRegression, TheilSen
from sklearn.linear_model.theilsen import spatial_median, modweiszfeld_step,\
    breakdown_point


def gen_toy_problem_1d():
    np.random.seed(0)
    n_samples = 100
    # Linear model y = 3*x + N(2, 0.1**2)
    x = np.random.randn(n_samples)
    w = np.array([3.])
    c = np.array([2.])
    noise = 0.1 * np.random.randn(n_samples)
    y = w * x + c + noise
    # Add some outliers
    x[42], y[42] = (-2, 4)
    x[43], y[43] = (-2.5, 8)
    x[53], y[53] = (2.5, 1)
    x[60], y[60] = (2.1, 2)
    x[72], y[72] = (1.8, -7)
    return x[:, np.newaxis], y, w, c


def gen_toy_problem_1d_no_intercept():
    np.random.seed(0)
    n_samples = 100
    # Linear model y = 3*x + N(2, 0.1**2)
    x = np.abs(np.random.randn(n_samples))
    w = np.array([3.])
    c = np.array([0.1])
    noise = 0.1 * np.random.randn(n_samples)
    y = w * x + c + noise
    # Add some outliers
    x[42], y[42] = (-2, 4)
    x[43], y[43] = (-2.5, 8)
    x[53], y[53] = (2.5, 1)
    x[60], y[60] = (2.1, 2)
    x[72], y[72] = (1.8, -7)
    return x[:, np.newaxis], y, w, c


def gen_toy_problem_2d():
    np.random.seed(0)
    n_samples=100
    # Linear model y = 5*x_1 + 10*x_2 + N(1, 0.1**2)
    X = np.random.randn(2*n_samples).reshape(n_samples, 2)
    w = np.array([5., 10.])
    c = np.array([1.])
    noise = 0.1*np.random.randn(n_samples)
    y = np.dot(X, w) + c + noise
    # Add some outliers
    n_outliers = n_samples // 10
    ix = np.random.randint(0, n_samples, n_outliers)
    y[ix] = 50*np.random.randn(n_outliers)
    return X, y, w, c


def gen_toy_problem_4d():
    np.random.seed(0)
    n_samples = 10000
    # Linear model y = 5*x_1 + 10*x_2  + 42*x_3 + 7*x_4 + N(1, 0.1**2)
    X = np.random.randn(4*n_samples).reshape(n_samples, 4)
    w = np.array([5., 10., 42., 7.])
    c = np.array([1.])
    noise = 0.1*np.random.randn(n_samples)
    y = np.dot(X, w) + c + noise
    # Add some outliers
    n_outliers = n_samples // 10
    ix = np.random.randint(0, n_samples, n_outliers)
    y[ix] = 50*np.random.randn(n_outliers)
    return X, y, w, c


def test_modweiszfeld_step_1d():
    X = np.array([1., 2., 3.]).reshape(3, 1)
    # Check startvalue is element of X and solution
    median = np.array([2.])
    new_y = modweiszfeld_step(X, median)
    nptest.assert_array_almost_equal(new_y, median)
    # Check startvalue is not the solution
    y = np.array([2.5])
    new_y = modweiszfeld_step(X, y)
    nptest.assert_array_less(median, new_y)
    nptest.assert_array_less(new_y, y)
    # Check startvalue is not the solution but element of X
    y = np.array([3.])
    new_y = modweiszfeld_step(X, y)
    nptest.assert_array_less(median, new_y)
    nptest.assert_array_less(new_y, y)
    # Check that a single vector is identity
    X = np.array([1., 2., 3.]).reshape(1, 3)
    y = X[0, ]
    new_y = modweiszfeld_step(X, y)
    nptest.assert_array_equal(y, new_y)


def test_modweiszfeld_step_2d():
    X = np.array([0., 0., 1., 1., 0., 1.]).reshape(3, 2)
    y = np.array([0.5, 0.5])
    # Check first two iterations
    new_y = modweiszfeld_step(X, y)
    nptest.assert_array_almost_equal(new_y, np.array([1 / 3, 2 / 3]))
    new_y = modweiszfeld_step(X, new_y)
    nptest.assert_array_almost_equal(new_y, np.array([0.2792408, 0.7207592]))
    # Check fix point
    y = np.array([0.21132505, 0.78867497])
    new_y = modweiszfeld_step(X, y)
    nptest.assert_array_almost_equal(new_y, y)


def test_spatial_median_1d():
    X = np.array([1., 2., 3.]).reshape(3, 1)
    true_median = np.array([2.])
    median = spatial_median(X)
    nptest.assert_array_almost_equal(median, true_median)
    # Check when maximum iteration is exceeded
    logging.basicConfig(filename=devnull)
    median = spatial_median(X, n_iter=30, tol=0.)
    nptest.assert_array_almost_equal(median, true_median)


def test_spatial_median_2d():
    X = np.array([0., 0., 1., 1., 0., 1.]).reshape(3, 2)
    median = spatial_median(X, n_iter=100, tol=1.e-6)

    def cost_func(y):
        dists = np.array([norm(x - y) for x in X])
        return np.sum(dists)

    # Check if median is solution of the Fermat-Weber location problem
    fermat_weber = minimize(cost_func, median)['x']
    nptest.assert_array_almost_equal(median, fermat_weber)


def test_theilsen_1d():
    X, y, w, c = gen_toy_problem_1d()
    # Check that Least Squares fails
    lstq = LinearRegression().fit(X, y)
    assert np.abs(lstq.coef_ - w) > 0.9
    # Check that Theil-Sen works
    theilsen = TheilSen().fit(X, y)
    nptest.assert_array_almost_equal(theilsen.coef_, w, 2)
    nptest.assert_array_almost_equal(theilsen.intercept_, c, 2)


def test_theilsen_1d_no_intercept():
    X, y, w, c = gen_toy_problem_1d_no_intercept()
    # Check that Least Squares fails
    lstq = LinearRegression(fit_intercept=False).fit(X, y)
    assert np.abs(lstq.coef_ - w - c) > 0.5
    # Check that Theil-Sen works
    theilsen = TheilSen(fit_intercept=False).fit(X, y)
    nptest.assert_array_almost_equal(theilsen.coef_, w + c, 1)
    assert theilsen.intercept_ == 0.


def test_theilsen_2d():
    X, y, w, c = gen_toy_problem_2d()
    # Check that Least Squares fails
    lstq = LinearRegression().fit(X, y)
    assert np.linalg.norm(lstq.coef_ - w) > 1.0
    # Check that Theil-Sen works
    theilsen = TheilSen().fit(X, y)
    nptest.assert_array_almost_equal(theilsen.coef_, w, 1)
    nptest.assert_array_almost_equal(theilsen.intercept_, c, 1)


def test_calc_breakdown_point():
    bp = breakdown_point(1e10, 2)
    assert np.abs(bp - (1 - 1/(np.sqrt(2)))) <= 1.e-6


@raises(AssertionError)
def test__checksubparams_too_large_subpopulation():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(n_subpopulation=100000).fit(X, y)


@raises(AssertionError)
def test__checksubparams_too_few_subsamples():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(n_subsamples=1).fit(X, y)


@raises(AssertionError)
def test__checksubparams_too_many_subsamples():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(n_subsamples=101).fit(X, y)


def test_subpopulation():
    X, y, w, c = gen_toy_problem_4d()
    theilsen = TheilSen(n_subpopulation=1000, random_state=0).fit(X, y)
    nptest.assert_array_almost_equal(theilsen.coef_, w, 1)
    nptest.assert_array_almost_equal(theilsen.intercept_, c, 1)


def test_subsamples():
    X, y, w, c = gen_toy_problem_4d()
    theilsen = TheilSen(n_subsamples=X.shape[0]).fit(X, y)
    lstq = LinearRegression().fit(X, y)
    # Check for exact the same results as Least Squares
    nptest.assert_array_almost_equal(theilsen.coef_, lstq.coef_, 9)
