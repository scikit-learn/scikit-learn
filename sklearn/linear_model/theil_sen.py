# -*- coding: utf-8 -*-
"""
A Theil-Sen Estimator for Multiple Linear Regression Model
"""

# Author: Florian Wilhelm <florian.wilhelm@gmail.com>
#
# License: BSD 3 clause

from __future__ import division, print_function, absolute_import

import warnings
from itertools import combinations

import numpy as np
from scipy import linalg
from scipy.special import binom
from scipy.linalg.lapack import get_lapack_funcs

from .base import LinearModel
from ..base import RegressorMixin
from ..utils import check_array, check_random_state, ConvergenceWarning
from ..utils import check_consistent_length, _get_n_jobs
from ..utils.random import choice
from ..externals.joblib import Parallel, delayed
from ..externals.six.moves import xrange as range

_EPSILON = np.finfo(np.double).eps


def _modified_weiszfeld_step(X, x_old):
    """Modified Weiszfeld step.

    This function defines one iteration step in order to approximate the
    spatial median (L1 median). It is a form of an iteratively re-weighted
    least squares method.

    Parameters
    ----------
    X : array, shape = [n_samples, n_features]
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    x_old : array, shape = [n_features]
        Current start vector.

    Returns
    -------
    x_new : array, shape = [n_features]
        New iteration step.

    References
    ----------
    - On Computation of Spatial Median for Robust Data Mining, 2005
      T. Kärkkäinen and S. Äyrämö
      http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
    """
    diff = X - x_old
    diff_norm = np.sqrt(np.sum(diff ** 2, axis=1))
    mask = diff_norm >= _EPSILON
    # x_old equals one of our samples
    is_x_old_in_X = int(mask.sum() < X.shape[0])

    diff = diff[mask]
    diff_norm = diff_norm[mask][:, np.newaxis]
    quotient_norm = linalg.norm(np.sum(diff / diff_norm, axis=0))

    if quotient_norm > _EPSILON:  # to avoid division by zero
        new_direction = (np.sum(X[mask, :] / diff_norm, axis=0)
                         / np.sum(1 / diff_norm, axis=0))
    else:
        new_direction = 1.
        quotient_norm = 1.

    return (max(0., 1. - is_x_old_in_X / quotient_norm) * new_direction
            + min(1., is_x_old_in_X / quotient_norm) * x_old)


def _spatial_median(X, max_iter=300, tol=1.e-3):
    """Spatial median (L1 median).

    The spatial median is member of a class of so-called M-estimators which
    are defined by an optimization problem. Given a number of p points in an
    n-dimensional space, the point x minimizing the sum of all distances to the
    p other points is called spatial median.

    Parameters
    ----------
    X : array, shape = [n_samples, n_features]
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    max_iter : int, optional
        Maximum number of iterations.  Default is 300.

    tol : float, optional
        Stop the algorithm if spatial_median has converged. Default is 1.e-3.

    Returns
    -------
    spatial_median : array, shape = [n_features]
        Spatial median.

    n_iter: int
        Number of iterations needed.

    References
    ----------
    - On Computation of Spatial Median for Robust Data Mining, 2005
      T. Kärkkäinen and S. Äyrämö
      http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf
    """
    if X.shape[1] == 1:
        return 1, np.median(X.ravel())

    tol **= 2  # We are computing the tol on the squared norm
    spatial_median_old = np.mean(X, axis=0)

    for n_iter in range(max_iter):
        spatial_median = _modified_weiszfeld_step(X, spatial_median_old)
        if np.sum((spatial_median_old - spatial_median) ** 2) < tol:
            break
        else:
            spatial_median_old = spatial_median
    else:
        warnings.warn("Maximum number of iterations {max_iter} reached in "
                      "spatial median for TheilSen regressor."
                      "".format(max_iter=max_iter), ConvergenceWarning)

    return n_iter, spatial_median


def _breakdown_point(n_samples, n_subsamples):
    """Approximation of the breakdown point.

    Parameters
    ----------
    n_samples : int
        Number of samples.

    n_subsamples : int
        Number of subsamples to consider.

    Returns
    -------
    breakdown_point : float
        Approximation of breakdown point.
    """
    return 1 - (0.5 ** (1 / n_subsamples) * (n_samples - n_subsamples + 1) +
                n_subsamples - 1) / n_samples


def _lstsq(X, y, indices, fit_intercept):
    """Least Squares Estimator for TheilSenRegressor class.

    This function calculates the least squares method on a subset of rows of X
    and y defined by the indices array. Optionally, an intercept column is
    added if intercept is set to true.

    Parameters
    ----------
    X : array, shape = [n_samples, n_features]
        Design matrix, where n_samples is the number of samples and
        n_features is the number of features.

    y : array, shape = [n_samples]
        Target vector, where n_samples is the number of samples.

    indices : array, shape = [n_subpopulation, n_subsamples]
        Indices of all subsamples with respect to the chosen subpopulation.

    fit_intercept : bool
        Fit intercept or not.

    Returns
    -------
    weights : array, shape = [n_subpopulation, n_features + intercept]
        Solution matrix of n_subpopulation solved least square problems.
    """
    fit_intercept = int(fit_intercept)
    n_features = X.shape[1] + fit_intercept
    n_subsamples = indices.shape[1]
    weights = np.empty((indices.shape[0], n_features))
    X_subpopulation = np.ones((n_subsamples, n_features))
    # gelss need to pad y_subpopulation to be of the max dim of X_subpopulation
    y_subpopulation = np.zeros((max(n_subsamples, n_features)))
    lstsq, = get_lapack_funcs(('gelss',), (X_subpopulation, y_subpopulation))

    for index, subset in enumerate(indices):
        X_subpopulation[:, fit_intercept:] = X[subset, :]
        y_subpopulation[:n_subsamples] = y[subset]
        weights[index] = lstsq(X_subpopulation,
                               y_subpopulation)[1][:n_features]

    return weights


class TheilSenRegressor(LinearModel, RegressorMixin):
    """Theil-Sen Estimator: robust multivariate regression model.

    The algorithm calculates least square solutions on subsets with size
    n_subsamples of the samples in X. Any value of n_subsamples between the
    number of features and samples leads to an estimator with a compromise
    between robustness and efficiency. Since the number of least square
    solutions is "n_samples choose n_subsamples", it can be extremely large
    and can therefore be limited with max_subpopulation. If this limit is
    reached, the subsets are chosen randomly. In a final step, the spatial
    median (or L1 median) is calculated of all least square solutions.

    Parameters
    ----------
    fit_intercept : boolean, optional, default True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    max_subpopulation : int, optional, default 1e4
        Instead of computing with a set of cardinality 'n choose k', where n is
        the number of samples and k is the number of subsamples (at least
        number of features), consider only a stochastic subpopulation of a
        given maximal size if 'n choose k' is larger than max_subpopulation.
        For other than small problem sizes this parameter will determine
        memory usage and runtime if n_subsamples is not changed.

    n_subsamples : int, optional, default None
        Number of samples to calculate the parameters. This is at least the
        number of features (plus 1 if fit_intercept=True) and the number of
        samples as a maximum. A lower number leads to a higher breakdown
        point and a low efficiency while a high number leads to a low
        breakdown point and a high efficiency. If None, take the
        minimum number of subsamples leading to maximal robustness.
        If n_subsamples is set to n_samples, Theil-Sen is identical to least
        squares.

    max_iter : int, optional, default 300
        Maximum number of iterations for the calculation of spatial median.

    tol : float, optional, default 1.e-3
        Tolerance when calculating spatial median.

    random_state : RandomState or an int seed, optional, default None
        A random number generator instance to define the state of the
        random permutations generator.

    n_jobs : integer, optional, default 1
        Number of CPUs to use during the cross validation. If ``-1``, use
        all the CPUs.

    verbose : boolean, optional, default False
        Verbose mode when fitting the model.

    Attributes
    ----------
    `coef_` : array, shape = (n_features)
        Coefficients of the regression model (median of distribution).

    `intercept_` : float
        Estimated intercept of regression model.

    `breakdown_` : float
        Approximated breakdown point.

    `n_iter_` : int
        Number of iterations needed for the spatial median.

    n_subpopulation_ : int
        Number of combinations taken into account from 'n choose k', where n is
        the number of samples and k is the number of subsamples.

    References
    ----------
    - Theil-Sen Estimators in a Multiple Linear Regression Model, 2009
      Xin Dang, Hanxiang Peng, Xueqin Wang and Heping Zhang
      http://www.math.iupui.edu/~hpeng/MTSE_0908.pdf
    """

    def __init__(self, fit_intercept=True, copy_X=True,
                 max_subpopulation=1e4, n_subsamples=None, max_iter=300,
                 tol=1.e-3, random_state=None, n_jobs=1, verbose=False):
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.max_subpopulation = int(max_subpopulation)
        self.n_subsamples = n_subsamples
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

    def _check_subparams(self, n_samples, n_features):
        n_subsamples = self.n_subsamples

        if self.fit_intercept:
            n_dim = n_features + 1
        else:
            n_dim = n_features

        if n_subsamples is not None:
            if n_subsamples > n_samples:
                raise ValueError("Invalid parameter since n_subsamples > "
                                 "n_samples ({0} > {1}).".format(n_subsamples,
                                                                 n_samples))
            if n_samples >= n_features:
                if n_dim > n_subsamples:
                    plus_1 = "+1" if self.fit_intercept else ""
                    raise ValueError("Invalid parameter since n_features{0} "
                                     "> n_subsamples ({1} > {2})."
                                     "".format(plus_1, n_dim, n_samples))
            else:  # if n_samples < n_features
                if n_subsamples != n_samples:
                    raise ValueError("Invalid parameter since n_subsamples != "
                                     "n_samples ({0} != {1}) while n_samples "
                                     "< n_features.".format(n_subsamples,
                                                            n_samples))
        else:
            n_subsamples = min(n_dim, n_samples)

        if self.max_subpopulation <= 0:
            raise ValueError("Subpopulation must be strictly positive "
                             "({0} <= 0).".format(self.max_subpopulation))

        all_combinations = max(1, np.rint(binom(n_samples, n_subsamples)))
        n_subpopulation = int(min(self.max_subpopulation, all_combinations))

        return n_subsamples, n_subpopulation

    def fit(self, X, y):
        """Fit linear model.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]
            Training data
        y : numpy array of shape [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        random_state = check_random_state(self.random_state)
        X = check_array(X)
        y = check_array(y, ensure_2d=False)
        check_consistent_length(X, y)
        n_samples, n_features = X.shape
        n_subsamples, self.n_subpopulation_ = self._check_subparams(n_samples,
                                                                    n_features)
        self.breakdown_ = _breakdown_point(n_samples, n_subsamples)

        if self.verbose:
            print("Breakdown point: {0}".format(self.breakdown_))
            print("Number of samples: {0}".format(n_samples))
            tol_outliers = int(self.breakdown_ * n_samples)
            print("Tolerable outliers: {0}".format(tol_outliers))
            print("Number of subpopulations: {0}".format(
                self.n_subpopulation_))

        # Determine indices of subpopulation
        if np.rint(binom(n_samples, n_subsamples)) <= self.max_subpopulation:
            indices = list(combinations(range(n_samples), n_subsamples))
        else:
            indices = [choice(n_samples,
                              size=n_subsamples,
                              replace=False,
                              random_state=random_state)
                       for _ in range(self.n_subpopulation_)]

        n_jobs = _get_n_jobs(self.n_jobs)
        index_list = np.array_split(indices, n_jobs)
        weights = Parallel(n_jobs=n_jobs,
                           verbose=self.verbose)(
            delayed(_lstsq)(X, y, index_list[job], self.fit_intercept)
            for job in range(n_jobs))
        weights = np.vstack(weights)
        self.n_iter_, coefs = _spatial_median(weights,
                                              max_iter=self.max_iter,
                                              tol=self.tol)

        if self.fit_intercept:
            self.intercept_ = coefs[0]
            self.coef_ = coefs[1:]
        else:
            self.intercept_ = 0.
            self.coef_ = coefs

        return self
