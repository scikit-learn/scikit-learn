"""
Covariance estimators

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
    """Covariance estimator

    Parameters
    ----------
    store_covariance : bool
        Specify if the estimated covariance is stored
    store_precision : bool
        Specify if the estimated precision is stored

    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix
        (stored only if store_covariance is True)

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated precision matrix

    """
    def __init__(self, store_covariance=True, store_precision=True):
        self.store_covariance = store_covariance
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
        if self.store_covariance:
            self.covariance_ = covariance
        else:
            self.covariance_ = None
        # set precision
        if self.store_precision:
            try:
                self.precision_ = linalg.inv(covariance)
            except:
                # /!\ add warning message or relevant exception handling here
                self.precision_ = None
        else:
            self.precision_ = None
    
    def fit(self, X, **params):
        self._set_params(**params)
        n_samples = X.shape[0]
        covariance_ = np.dot(X.T, X) / n_samples
        self._set_estimates(covariance_)
        
        return self


    def score(self, X_test):
        n_samples = X_test.shape[0]
        test_cov = np.dot(X_test.T, X_test) / n_samples
        return self.log_likelihood(test_cov)


    def log_likelihood(self, test_cov):
        return -np.sum(test_cov*self.precision_) + \
                                            exact_logdet(self.precision_)

