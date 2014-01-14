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
from numpy.linalg import norm, solve
from scipy.special import binom
from .base import LinearModel
from ..base import RegressorMixin
from ..utils import check_arrays

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
    T = T_nom / T_denom
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


class TheilSen(LinearModel, RegressorMixin):
    """Theil-Sen Estimator for a multiple linear regression model.

    Parameters
    ----------
    n_iter : int, optional, default 300
        Maximum number of iterations for the calculation of spatial median.

    tol : float, optional, default 1.e-3
        Tolerance when calculating spatial median.

    n_samples : int, optional, default 1e6
        Number of samples to calculate the parameters.

    fit_intercept : boolean, optional, default True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.


    Attributes
    ----------
    `coef_` : array, shape = (n_features)
        Coefficients of the regression model (median of distribution).

    `intercept_` : float
       Estimated intercept of regression model.
    """

    def __init__(self, n_iter=300, n_samples=1e6, tol=1.e-3,
                 fit_intercept=True, copy_X=True):
        self.n_iter = n_iter
        self.n_samples = int(n_samples)
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X

    def fit(self, X, y):
        X, y = check_arrays(X, y, sparse_format='dense', dtype=np.float)
        n_samples, n_features = X.shape
        if self.fit_intercept:
            n_dim = n_features + 1
            fst = 1
        else:
            n_dim = n_features
            fst = 0
        indices = combinations(xrange(n_samples), n_dim)
        weights = np.empty((int(binom(n_samples, n_dim)), n_dim))
        for i, ix in enumerate(indices):
            X_sub = np.ones((n_dim, n_dim))
            X_sub[:, fst:] = X[list(ix), :]
            weights[i, :] = solve(X_sub, y[list(ix)])

        coefs = spatial_median(weights, n_iter=self.n_iter, tol=self.tol)

        if self.fit_intercept:
            self.intercept_ = coefs[0]
            self.coef_ = coefs[1:]
        else:
            self.intercept_ = 0.
            self.coef_ = coefs

        return self
