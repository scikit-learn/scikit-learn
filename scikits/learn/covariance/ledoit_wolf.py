"""
Basic covariance estimators
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np
from scipy import linalg

from ..base import BaseEstimator
from ..utils.extmath import fast_logdet as exact_logdet

################################################################################
# Covariance estimator

class Covariance(BaseEstimator):
    """Basic covariance estimator

    Parameters
    ----------
    store_covariance : bool
        Specify if the estimated covariance is stored

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix
        (stored only is store_covariance is True)

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated precision matrix

    """
    def __init__(self, store_covariance=True):
        self.store_covariance = True


    def _set_precision(self, covariance):
        if self.store_covariance:
            self.covariance_ = covariance
        else:
            self.covariance_ = None
        self.precision_ = linalg.inv(covariance)


    def fit(self, X, **params):
        self._set_params(**params)
        n_samples = X.shape[0]
        covariance_ = np.dot(X.T, X) / n_samples
        self._set_precision(covariance_)
        return self


    def score(self, X_test):
        n_samples = X_test.shape[0]
        test_cov = np.dot(X_test.T, X_test) / n_samples
        return self.log_likelihood(test_cov)


    def log_likelihood(self, test_cov):
        return -np.sum(test_cov*self.precision_) + \
                                            exact_logdet(self.precision_)


################################################################################
# ShrunkCovariance estimator

def shrunk_covariance(X, shrinkage=0.1, data_is_cov=False):
    """ Calculate a covariance matrix shrunk on the diagonal

    Returns
    -------
    regularised_cov: 2D ndarray
        Regularized covariance

    Notes
    -----
    The regularized covariance is given by

        (1 - shrinkage)*cov
                + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features

    """
    n_samples, n_features = X.shape
    if data_is_cov:
        emp_cov = X
    else:
        emp_cov = np.dot(X.T, X) / n_samples
    mu = np.trace(emp_cov) / n_features
    shrunk_cov = (1.-shrinkage)*emp_cov + shrinkage*mu*np.eye(n_features)
    return shrunk_cov


class ShrunkCovariance(Covariance):
    """Covariance estimator with shrinkage

    Parameters
    ----------
    store_covariance : bool
        Specify if the estimated covariance is stored

    shrinkage : float
        Shrinkage (in [0, 1])

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix
        (stored only is store_covariance is True)

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated precision matrix

    Notes
    -----
    The regularized covariance is given by

        (1 - shrinkage)*cov
                + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features

    """
    def __init__(self, store_covariance=True, shrinkage=None):
        self.store_covariance = True
        self.shrinkage = shrinkage


    def fit(self, X, **params):
        self._set_params(**params)
        covariance_ = shrunk_covariance(X, self.shrinkage,
                                                    data_is_cov=False)
        self._set_precision(covariance_)
        return self


################################################################################
# Ledoit-Wolf estimator

def ledoit_wolf(X, return_shrinkage=False):
    """ Estimates the shrunk Ledoit-Wolf covariance matrix.

        Parameters
        ----------
        X: 2D ndarray, shape (n, p)
            The data matrix, with p features and n samples.
        return_shrinkage: boolean, optional
            If return_shrinkage is True, the regularisation_factor is
            returned.

        Returns
        -------
        regularised_cov: 2D ndarray
            Regularized covariance
        shrinkage: float
            Regularisation factor

        Notes
        -----
        The regularised covariance is::

            (1 - shrinkage)*cov
                    + shrinkage * mu * np.identity(n_features)

        where mu = trace(cov) / n_features
    """
    n_samples, n_features = X.shape
    if n_features == 1:
        if return_shrinkage:
            return np.atleast_2d(X.std()), 0
        return np.atleast_2d(X.std())
    cov = np.dot(X.T, X) / n_samples
    i = np.identity(n_features)
    mu = np.trace(cov) / n_features
    delta = ((cov - mu*i)**2).sum() / n_features
    X2 = X**2
    beta_ = 1./(n_features*n_samples) * np.sum(
                            np.dot(X2.T, X2)/n_samples - cov**2
                )

    beta = min(beta_, delta)
    alpha = delta - beta
    if not return_shrinkage:
        return beta/delta * mu * i + alpha/delta * cov
    else:
        return beta/delta * mu * i + alpha/delta * cov, beta/delta


class LedoitWolf(Covariance):
    """
    LedoitWolf Estimator

    Parameters
    ----------
    store_covariance : bool
        Specify if the estimated covariance is stored

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix
        (stored only is store_covariance is True)

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated precision matrix

    `shrinkage_` : float
        Scalar used to regularize the precision matrix estimation

    Notes
    -----
    The regularised covariance is::

        (1 - shrinkage)*cov
                + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features

    Reference :
    "A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices",
    Ledoit and Wolf, Journal of Multivariate Analysis, Volume 88, Issue 2,
    February 2004, pages 365-411.
    """

    def __init__(self, store_covariance=True):
        super(LedoitWolf, self).__init__(store_covariance)


    def fit(self, X):
        covariance_, shrinkage = ledoit_wolf(X, return_shrinkage=True)
        self.shrinkage_ = shrinkage
        self._set_precision(covariance_)
        return self


