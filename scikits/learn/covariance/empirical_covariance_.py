"""
Maximum likelihood covariance estimator.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np
from scipy import linalg

from ..base import BaseEstimator
from ..utils.extmath import fast_logdet as exact_logdet


def log_likelihood(emp_cov, precision):
    """Computes the negative log_likelihood of the data

    Params
    ------
    emp_cov: 2D ndarray (n_features, n_features)
      Maximum Likelihood Estimator of covariance
    precision: 2D ndarray (n_features, n_features)
      The precision matrix of the covariance model to be tested

    """
    return -np.sum(emp_cov*precision) + exact_logdet(precision)


class EmpiricalCovariance(BaseEstimator):
    """Maximum likelihood covariance estimator

    Parameters
    ----------
    store_precision : bool
        Specifies if the estimated precision is stored

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    """
    def __init__(self, store_precision=True):
        self.store_precision = store_precision
    
    def _set_estimates(self, covariance):
        """Saves the covariance and precision estimates
        
        Storage is done accordingly to `self.store_precision`.
        Precision stored only if invertible.

        Params
        ------
        covariance: 2D ndarray, shape (n_features, n_features)
          Estimated covariance matrix to be stored, and from which the precision
          is computed.

        """
        covariance = np.atleast_2d(covariance)
        # set covariance
        self.covariance_ = covariance
        # set precision
        if self.store_precision:
            self.precision_ = linalg.pinv(covariance)
        else:
            self.precision_ = None
    
    def fit(self, X, assume_centered=False, **params):
        self._set_params(**params)
        if assume_centered:
            covariance = np.dot(X.T, X) / X.shape[0]
        else:
            covariance = np.cov(X.T, bias=1)
        self._set_estimates(covariance)
        
        return self

    def score(self, X_test, assume_centered=False):
        # compute empirical covariance of the test set
        if assume_centered:
            test_cov = np.dot(X_test.T, X_test) / X_test.shape[0]
        else:
            test_cov = np.cov(X_test.T, bias=1)
        # compute log likelihood
        if self.store_precision:
            res = log_likelihood(test_cov, self.precision_)
        else:
            res = log_likelihood(test_cov, linalg.pinv(self.covariance_))
        
        return res
        
    def mse(self, real_cov):
        diff = real_cov - self.covariance_
        
        return np.sum(diff**2)
