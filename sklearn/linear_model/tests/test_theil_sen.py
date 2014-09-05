"""
Testing for Theil-Sen module (sklearn.linear_model.theil_sen)
"""

# Author: Florian Wilhelm <florian.wilhelm@gmail.com>
# License: BSD 3 clause

from __future__ import division, print_function, absolute_import

import logging
from os import devnull

import numpy as np
from numpy.testing import assert_array_equal, assert_array_less, \
    assert_array_almost_equal, assert_warns
from scipy.linalg import norm
from scipy.optimize import fmin_bfgs
from nose.tools import raises
from sklearn.utils import ConvergenceWarning
from sklearn.externals.joblib import cpu_count
from sklearn.externals.six.moves import xrange
from sklearn.linear_model import LinearRegression, TheilSen
from sklearn.linear_model.theil_sen import _spatial_median, _breakdown_point, \
    _modweiszfeld_step


def gen_toy_problem_1d(intercept=True):
    random_state = np.random.RandomState(0)
    if intercept:
        n_samples = 50
    else:
        n_samples = 100
    # Linear model y = 3*x + N(2, 0.1**2)
    x = random_state.normal(size=n_samples)
    w = 3.
    if intercept:
        c = 2.
    else:
        c = 0.1
    noise = 0.1 * random_state.normal(size=n_samples)
    y = w * x + c + noise
    # Add some outliers
    if intercept:
        x[42], y[42] = (-2, 4)
        x[43], y[43] = (-2.5, 8)
        x[33], y[33] = (2.5, 1)
        x[49], y[49] = (2.1, 2)
    else:
        x[42], y[42] = (-2, 4)
        x[43], y[43] = (-2.5, 8)
        x[53], y[53] = (2.5, 1)
        x[60], y[60] = (2.1, 2)
        x[72], y[72] = (1.8, -7)
    return x[:, np.newaxis], y, w, c


def gen_toy_problem_2d():
    random_state = np.random.RandomState(0)
    n_samples = 100
    # Linear model y = 5*x_1 + 10*x_2 + N(1, 0.1**2)
    X = random_state.normal(size=(n_samples, 2))
    w = np.array([5., 10.])
    c = 1.
    noise = 0.1*random_state.normal(size=n_samples)
    y = np.dot(X, w) + c + noise
    # Add some outliers
    n_outliers = n_samples // 10
    ix = random_state.randint(0, n_samples, size=n_outliers)
    y[ix] = 50*random_state.normal(size=n_outliers)
    return X, y, w, c


def gen_toy_problem_4d():
    random_state = np.random.RandomState(0)
    n_samples = 10000
    # Linear model y = 5*x_1 + 10*x_2  + 42*x_3 + 7*x_4 + N(1, 0.1**2)
    X = random_state.normal(size=(n_samples, 4))
    w = np.array([5., 10., 42., 7.])
    c = 1.
    noise = 0.1*random_state.normal(size=n_samples)
    y = np.dot(X, w) + c + noise
    # Add some outliers
    n_outliers = n_samples // 10
    ix = random_state.randint(0, n_samples, size=n_outliers)
    y[ix] = 50*random_state.normal(size=n_outliers)
    return X, y, w, c


def test_modweiszfeld_step_1d():
    X = np.array([1., 2., 3.]).reshape(3, 1)
    # Check startvalue is element of X and solution
    median = 2.
    new_y = _modweiszfeld_step(X, median)
    assert_array_almost_equal(new_y, median)
    # Check startvalue is not the solution
    y = 2.5
    new_y = _modweiszfeld_step(X, y)
    assert_array_less(median, new_y)
    assert_array_less(new_y, y)
    # Check startvalue is not the solution but element of X
    y = 3.
    new_y = _modweiszfeld_step(X, y)
    assert_array_less(median, new_y)
    assert_array_less(new_y, y)
    # Check that a single vector is identity
    X = np.array([1., 2., 3.]).reshape(1, 3)
    y = X[0, ]
    new_y = _modweiszfeld_step(X, y)
    assert_array_equal(y, new_y)


def test_modweiszfeld_step_2d():
    X = np.array([0., 0., 1., 1., 0., 1.]).reshape(3, 2)
    y = np.array([0.5, 0.5])
    # Check first two iterations
    new_y = _modweiszfeld_step(X, y)
    assert_array_almost_equal(new_y, np.array([1 / 3, 2 / 3]))
    new_y = _modweiszfeld_step(X, new_y)
    assert_array_almost_equal(new_y, np.array([0.2792408, 0.7207592]))
    # Check fix point
    y = np.array([0.21132505, 0.78867497])
    new_y = _modweiszfeld_step(X, y)
    assert_array_almost_equal(new_y, y)


def test_spatial_median_1d():
    X = np.array([1., 2., 3.]).reshape(3, 1)
    true_median = 2.
    median = _spatial_median(X)
    assert_array_almost_equal(median, true_median)
    # Check when maximum iteration is exceeded a warning is emitted
    assert_warns(ConvergenceWarning, _spatial_median, X, n_iter=30, tol=0.)


def test_spatial_median_2d():
    X = np.array([0., 0., 1., 1., 0., 1.]).reshape(3, 2)
    median = _spatial_median(X, n_iter=100, tol=1.e-6)

    def cost_func(y):
        dists = np.array([norm(x - y) for x in X])
        return np.sum(dists)

    # Check if median is solution of the Fermat-Weber location problem
    fermat_weber = fmin_bfgs(cost_func, median, disp=False)
    assert_array_almost_equal(median, fermat_weber)


def test_theil_sen_1d():
    X, y, w, c = gen_toy_problem_1d()
    # Check that Least Squares fails
    lstq = LinearRegression().fit(X, y)
    assert np.abs(lstq.coef_ - w) > 0.9
    # Check that Theil-Sen works
    theil_sen = TheilSen().fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, w, 1)
    assert_array_almost_equal(theil_sen.intercept_, c, 1)


def test_theil_sen_1d_no_intercept():
    X, y, w, c = gen_toy_problem_1d(intercept=False)
    # Check that Least Squares fails
    lstq = LinearRegression(fit_intercept=False).fit(X, y)
    assert np.abs(lstq.coef_ - w - c) > 0.5
    # Check that Theil-Sen works
    theil_sen = TheilSen(fit_intercept=False).fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, w + c, 1)
    assert theil_sen.intercept_ == 0.


def test_theil_sen_2d():
    X, y, w, c = gen_toy_problem_2d()
    # Check that Least Squares fails
    lstq = LinearRegression().fit(X, y)
    assert norm(lstq.coef_ - w) > 1.0
    # Check that Theil-Sen works
    theil_sen = TheilSen(max_subpopulation=1e3, random_state=0).fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, w, 1)
    assert_array_almost_equal(theil_sen.intercept_, c, 1)


def test_calc_breakdown_point():
    bp = _breakdown_point(1e10, 2)
    assert np.abs(bp - 1 + 1/(np.sqrt(2))) <= 1.e-6


@raises(ValueError)
def test_checksubparams_negative_subpopulation():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(max_subpopulation=-1).fit(X, y)


@raises(AssertionError)
def test_checksubparams_too_few_subsamples():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(n_subsamples=1).fit(X, y)


@raises(AssertionError)
def test_checksubparams_too_many_subsamples():
    X, y, w, c = gen_toy_problem_1d()
    TheilSen(n_subsamples=101).fit(X, y)


@raises(AssertionError)
def test_checksubparams_n_subsamples_if_less_samples_than_features():
    random_state = np.random.RandomState(0)
    n_samples, n_features = 10, 20
    X = random_state.normal(size=(n_samples, n_features))
    y = random_state.normal(size=n_samples)
    TheilSen(n_subsamples=9).fit(X, y)


def test_subpopulation():
    X, y, w, c = gen_toy_problem_4d()
    theil_sen = TheilSen(max_subpopulation=1000, random_state=0).fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, w, 1)
    assert_array_almost_equal(theil_sen.intercept_, c, 1)


def test_subsamples():
    X, y, w, c = gen_toy_problem_4d()
    theil_sen = TheilSen(n_subsamples=X.shape[0], verbose=True).fit(X, y)
    lstq = LinearRegression().fit(X, y)
    # Check for exact the same results as Least Squares
    assert_array_almost_equal(theil_sen.coef_, lstq.coef_, 9)


def test_verbosity():
    X, y, w, c = gen_toy_problem_1d()
    # Check that Theil-Sen can be verbose
    TheilSen(verbose=True).fit(X, y)
    TheilSen(verbose=True, max_subpopulation=10).fit(X, y)


def test_theil_sen_parallel():
    X, y, w, c = gen_toy_problem_2d()
    # Check that Least Squares fails
    lstq = LinearRegression().fit(X, y)
    assert norm(lstq.coef_ - w) > 1.0
    # Check that Theil-Sen works
    theil_sen = TheilSen(n_jobs=-1,
                         random_state=0,
                         max_subpopulation=2e3).fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, w, 1)
    assert_array_almost_equal(theil_sen.intercept_, c, 1)


@raises(ValueError)
def test_get_n_jobs_with_0_CPUs():
    TheilSen(n_jobs=0)._get_n_jobs()


def test_get_n_jobs():
    n_jobs = TheilSen(n_jobs=2)._get_n_jobs()
    assert n_jobs == 2
    n_jobs = TheilSen(n_jobs=-2)._get_n_jobs()
    assert n_jobs == cpu_count() - 1


def test_split_indices():
    random_state = np.random.RandomState(0)

    def gen_indices(n):
        for _ in xrange(n):
            yield random_state.randint(n, size=10)

    theil_sen = TheilSen()
    indices, starts = theil_sen._split_indices(gen_indices(10), 3)
    assert len(indices) == 3
    assert indices[0].shape[0] == 4
    assert indices[1].shape[0] == 3
    assert indices[2].shape[0] == 3
    assert_array_equal(starts, np.array([0, 4, 7]))

    # Test if more splits than elements in indices
    indices, starts = theil_sen._split_indices(gen_indices(2), 4)
    assert indices[0].shape[0] == 1
    assert indices[1].shape[0] == 1
    assert indices[2].shape[0] == 0
    assert indices[2].shape[0] == 0
    assert_array_equal(starts, np.array([0, 1, 2, 2]))


def test_less_samples_than_features():
    random_state = np.random.RandomState(0)
    n_samples, n_features = 10, 20
    X = random_state.normal(size=(n_samples, n_features))
    y = random_state.normal(size=n_samples)
    # Check that Theil-Sen falls back to Least Squares if fit_intercept=False
    theil_sen = TheilSen(fit_intercept=False).fit(X, y)
    lstq = LinearRegression(fit_intercept=False).fit(X, y)
    assert_array_almost_equal(theil_sen.coef_, lstq.coef_, 12)
    # Check fit_intercept=True case. This will not be equal to the Least
    # Squares solution since the intercept is calculated differently.
    theil_sen = TheilSen(fit_intercept=True).fit(X, y)
    y_pred = theil_sen.predict(X)
    assert_array_almost_equal(y_pred, y, 12)
