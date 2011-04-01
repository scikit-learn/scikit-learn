"""
OAS covariance estimator.

OAS is a particular form of shrinkage, where the shrinkage
coefficient is computed using the OSA formula as given in
"Shrinkage Algorithms for MMSE Covariance Estimation"
Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np

from .covariance import Covariance
from .shrunk_covariance import ShrunkCovariance

################################################################################
# OAS estimator

def oas(X):
    """Estimates the shrunk OAS covariance matrix.

    Parameters
    ----------
    X: 2D ndarray, shape (n_samples, n_features)
      Data from which to compute the covariance estimate
    
    Returns
    -------
    shrunk_cov: 2D ndarray, shape (n_features, n_features)
      Shrunk covariance
    shrinkage: float
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Notes
    -----
    The regularised (shrunk) covariance is:
    
    (1 - shrinkage)*cov
      + shrinkage * mu * np.identity(n_features)
    
    where mu = trace(cov) / n_features
    
    """
    X = np.atleast_2d(X)
    n_samples, n_features = X.shape
    
    emp_cov = Covariance(store_precision=False).fit(X).covariance_
    structured_estimate = np.eye(n_features)
    mu = np.trace(emp_cov) / n_features
    
    num = (((-1.)/n_features) * np.trace(emp_cov)**2) + np.sum(emp_cov**2)**2
    den = ((n_samples-1.)/n_features) \
        * (np.sum(emp_cov**2)**2 - ((np.trace(emp_cov)**2)/float(n_features)))

    shrinkage = min(num/den, 1.)
    shrunk_cov = shrinkage*mu*structured_estimate + (1.-shrinkage)*emp_cov
    
    return shrunk_cov, shrinkage


class OAS(ShrunkCovariance):
    """
    OAS Estimator

    Notes
    -----
    The regularised covariance is::

        (1 - shrinkage)*cov
                + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features
    and shinkage is given by the OAS formula (see Reference)

    Reference
    ---------
    "Shrinkage Algorithms for MMSE Covariance Estimation"
    Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
    
    """

    def __init__(self, store_covariance=True, store_precision=True):
        super(OAS, self).__init__(store_covariance, store_precision)


    def fit(self, X):
        covariance, shrinkage = oas(X)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)
        
        return self


