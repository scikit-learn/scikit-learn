"""
Covariance estimators using shrinkage.

Shrinkage corresponds to regularising `cov` using a convex combination:
shrunk_cov = (1-shrinkage)*cov + shrinkage*structured_estimate.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np

from .covariance import Covariance

###############################################################################
# ShrunkCovariance estimator

def shrunk_covariance(emp_cov, shrinkage=None):
    """Calculates a covariance matrix shrunk on the diagonal
    
    Params
    ------
    emp_cov: 2D ndarray, shape (n_features, n_features)
      Covariance matrix to be shrunk
    shrinkage: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Returns
    -------
    shrunk_cov: 2D ndarray
      Shrunk covariance
    
    Notes
    -----
    The regularized (shrunk) covariance is given by
    
    (1 - shrinkage)*cov
      + shrinkage*mu*np.identity(n_features)
    
    where mu = trace(cov) / n_features
    
    """
    emp_cov = np.atleast_2d(emp_cov)
    n_features = emp_cov.shape[0]

    if shrinkage is None:
        shrinkage = np.trace(emp_cov) / n_features
    
    mu = np.trace(emp_cov) / n_features
    shrunk_cov = (1.-shrinkage)*emp_cov + shrinkage*mu*np.eye(n_features)
    
    return shrunk_cov


class ShrunkCovariance(Covariance):
    """Covariance estimator with shrinkage
    
    Parameters
    ----------
    shrinkage: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Attributes
    ----------
    `shrinkage_`: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Notes
    -----
    The regularized covariance is given by
    
    (1 - shrinkage)*cov
      + shrinkage*mu*np.identity(n_features)
    
    where mu = trace(cov) / n_features
    
    """
    def __init__(self, store_covariance=True, store_precision=True,
                 shrinkage=None):
        super(ShrunkCovariance, self).__init__(store_covariance,
                                               store_precision)
        self.shrinkage_ = shrinkage


    def fit(self, X, **params):
        self._set_params(**params)
        emp_cov = Covariance(store_precision=False)
        emp_cov.fit(X)
        covariance = shrunk_covariance(emp_cov.covariance_, self.shrinkage_)
        self._set_estimates(covariance)
        
        return self
