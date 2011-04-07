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


def empirical_covariance(X):
    """Computes the maximum likelihood estimator of the covariance matrix	

    Parameters
    ----------
    X: 2D ndarray, shape (n_samples, n_features)
      Data from which to compute the covariance estimate	

    Returns
    -------	
    cov: 2D ndarray, shape (n_features, n_features)
      Maximum Likelihood Estimator of covariance

    """
    X = np.asanyarray(X, dtype=np.float64, order='C')
    if X.ndim == 1:
        X = np.atleast_2d(X).T

    return np.dot(X.T, X) / X.shape[0]

class EmpiricalCovariance(BaseEstimator):
    """Maximum likelihood covariance estimator

    Parameters
    ----------
    store_precision : bool
        Specify if the estimated precision is stored

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix
        (stored only if store_covariance is True)

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated pseudo inverse matrix

    """
    def __init__(self, store_precision=True):
        self.store_precision = store_precision
    
    def _set_estimates(self, covariance):
        """Saves the covariance and precision estimates
        
        Storage is done accordingly to `self.store_covariance`
        and `self.store_precision`.
        Precision stored only if invertible.

        Params
        ------
        covariance: 2D ndarray, shape (n_features, n_features)
          Estimated covariance matrix to be stored, and from which the precision
          is computed.

        """
        # set covariance
        self.covariance_ = covariance
        # set precision
        if self.store_precision:
            self.precision_ = linalg.pinv(covariance)
        else:
            self.precision_ = None
    
    def fit(self, X, **params):
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64, order='C')
        covariance = empirical_covariance(X)
        self._set_estimates(covariance)
        
        return self

    def score(self, X_test):
        X_test = np.asanyarray(X_test, dtype=np.float64, order='C')
        test_cov = empirical_covariance(X_test)
        if self.store_precision:
            res = log_likelihood(test_cov, self.precision_)
        else:
            res = log_likelihood(test_cov, linalg.pinv(self.covariance_))
        
        return res
        
    def mse(self, real_cov):
        diff = real_cov - self.covariance_
        
        return np.sum(diff**2)
