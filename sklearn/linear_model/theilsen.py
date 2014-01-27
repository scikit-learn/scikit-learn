"""
A Theil-Sen Estimator for Multiple Linear Regression Model
"""

# Author: Florian Wilhelm <florian.wilhelm@gmail.com>
#
# License: BSD 3 clause

from __future__ import division, print_function, absolute_import

import logging
from itertools import combinations

import numpy as np
from numpy.linalg import norm, lstsq
from scipy.special import binom
from .base import LinearModel
from ..base import RegressorMixin
from ..utils import check_arrays, check_random_state

_logger = logging.getLogger(__name__)


def modweiszfeld_step(X, y):
    """Modified Weiszfeld step.

    Parameters
    ----------
    X : array, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    y : array, shape = [n_features]
        Current start vector.

    Returns
    -------
    new_y : array, shape = [n_features]
        New iteration step.
    """
    X, y = check_arrays(X.T, y, sparse_format='dense', dtype=np.float)
    eta = 0.
    res = np.zeros(y.shape)
    T_nom = np.zeros(y.shape)
    T_denom = 0.
    for x in X.T:
        diff = x - y
        normdiff = norm(diff)
        if normdiff < 1.e-6:
            eta = 1.
            continue
        res += diff / normdiff
        T_denom += 1. / normdiff
        T_nom += x / normdiff
    if T_denom != 0.:
        T = T_nom / T_denom
    else:
        T = T_nom
    r = norm(res)
    if r < 1.e-6:
        r = 1.
    return max(0., 1. - eta/r)*T + min(1., eta/r)*y


def spatial_median(X, n_iter=300, tol=1.e-3):
    """Spatial median (L1 median).

    Parameters
    ----------
    X : array, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    n_iter : int, optional
        Maximum number of iterations.  Default is 300.

    tol : float, optional
        Stop the algorithm if spmed has converged. Default is 1.e-3.

    Returns
    -------
    spmed : array, shape = [n_features]
        Spatial median.
    """
    X = check_arrays(X, sparse_format='dense', dtype=np.float)[0]
    spmed_old = np.mean(X, axis=0)
    for _ in xrange(n_iter):
        spmed = modweiszfeld_step(X, spmed_old)
        if norm(spmed_old - spmed) < tol:
            return spmed
        else:
            spmed_old = spmed
    _logger.warn("Maximum number of iterations reached.")
    return spmed


def breakdown_point(n_samples, n_subsamples):
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
    return 1 - 0.5**(1/n_subsamples)*(n_samples - n_subsamples + 1)/n_samples


class TheilSen(LinearModel, RegressorMixin):
    """Theil-Sen Estimator for a multiple linear regression model.

    Parameters
    ----------
    fit_intercept : boolean, optional, default True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    n_subpopulation : int, optional, default None
        Instead of computing with a set of cardinality 'n choose k', where n is
        the number of samples and k is min_subsamples (at least number of
        features), consider only a stochastic subpopulation of cardinality
        n_subpopulation. If None, no stochastic sampling of subpopulation is
        done.

    min_subsamples : int, optional, default None
        Number of samples to calculate the parameters. This is at least the
        number of features (plus 1 if fit_intercept=True) and the number of
        samples as a maximum. A lower number leads to a higher breakdown
        point and a low efficiency while a high number leads to a low
        breakdown point and a high efficiency. If None, take the
        minimum number of subsamples leading to maximal robustness.

    n_iter : int, optional, default 300
        Maximum number of iterations for the calculation of spatial median.

    tol : float, optional, default 1.e-3
        Tolerance when calculating spatial median.

    random_state : RandomState or an int seed, optional, default None
        A random number generator instance to define the state of the
        random permutations generator.

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

    `random_state_` : RandomState
        The current random state.
    """

    def __init__(self, fit_intercept=True, copy_X=True,
                 n_subpopulation=None, n_subsamples=None, n_iter=300,
                 tol=1.e-3, random_state=None, verbose=False):
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.n_subpopulation = n_subpopulation
        self.n_subsamples = n_subsamples
        self.n_iter = n_iter
        self.tol = tol
        self.random_state = random_state
        self.verbose = verbose

    def _check_subparams(self, n_samples, n_features):
        if self.fit_intercept:
            n_dim = n_features + 1
            fst = 1
        else:
            n_dim = n_features
            fst = 0
        n_subpopulation = self.n_subpopulation
        n_subsamples = self.n_subsamples
        if n_subsamples is not None:
            assert n_dim <= n_subsamples <= n_samples
        else:
            n_subsamples = n_dim
        if n_subpopulation is not None:
            assert n_subpopulation <= binom(n_samples, n_subsamples)
        else:
            n_subpopulation = int(binom(n_samples, n_subsamples))

        return fst, n_dim, n_subsamples, n_subpopulation

    def _subpop_iter(self, n_samples, n_dim):
        for s in xrange(self.n_subpopulation):
            yield self.random_state_.randint(0, n_samples, n_dim)

    def _get_indices(self, n_samples, n_dim):
        if self.n_subpopulation is None:
            return combinations(xrange(n_samples), n_dim)
        else:
            return self._subpop_iter(n_samples, n_dim)

    def fit(self, X, y):
        self.random_state_ = check_random_state(self.random_state)
        X, y = check_arrays(X, y, sparse_format='dense', dtype=np.float)
        n_samples, n_features = X.shape
        fst, n_dim, n_ss, n_sp = self._check_subparams(n_samples, n_features)
        self.breakdown_ = breakdown_point(n_samples, n_dim)

        indices = self._get_indices(n_samples, n_ss)
        weights = np.empty((n_sp, n_dim))
        for i, ix in enumerate(indices):
            X_sub = np.ones((n_ss, n_dim))
            X_sub[:, fst:] = X[list(ix), :]
            weights[i, :] = lstsq(X_sub, y[list(ix)])[0]

        coefs = spatial_median(weights, n_iter=self.n_iter, tol=self.tol)

        if self.fit_intercept:
            self.intercept_ = coefs[0]
            self.coef_ = coefs[1:]
        else:
            self.intercept_ = 0.
            self.coef_ = coefs

        return self
