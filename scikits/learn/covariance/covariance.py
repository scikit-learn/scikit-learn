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

