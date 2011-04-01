"""
Ledoit-Wolf covariance estimator.

Ledoit-Wolf is a particular form of shrinkage, where the shrinkage
coefficient is computed using O.Ledoit and M.Wolf's formula as
described in
"A Well-Conditioned Estimator for Large-Dimensional
Covariance Matrices",
Ledoit and Wolf, Journal of Multivariate
Analysis, Volume 88, Issue 2, February 2004, pages 365-411.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np

from .covariance import Covariance
from .shrunk_covariance import ShrunkCovariance

################################################################################
# Ledoit-Wolf estimator

def ledoit_wolf(X):
    """Estimates the shrunk Ledoit-Wolf covariance matrix.

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
    
    if n_features == 1:
        return np.atleast_2d(X.std()), 0
    
    emp_cov = Covariance(store_precision=False).fit(X).covariance_
    structured_estimate = np.eye(n_features)
    mu = np.trace(emp_cov) / n_features
    delta = ((emp_cov - mu*structured_estimate)**2).sum() / n_features
    X2 = X**2
    beta_ = 1./(n_features*n_samples) \
        * np.sum(np.dot(X2.T, X2)/n_samples - emp_cov**2)

    beta = min(beta_, delta)
    shrinkage = beta/delta
    shrunk_cov = shrinkage*mu*structured_estimate + (1.-shrinkage)*emp_cov
    
    return shrunk_cov, shrinkage


class LedoitWolf(ShrunkCovariance):
    """
    LedoitWolf Estimator

    Notes
    -----
    The regularised covariance is::

        (1 - shrinkage)*cov
                + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features
    and shinkage is given by the Ledoit and Wolf formula (see Reference)

    Reference
    ---------
    "A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices",
    Ledoit and Wolf, Journal of Multivariate Analysis, Volume 88, Issue 2,
    February 2004, pages 365-411.
    
    """

    def __init__(self, store_covariance=True, store_precision=True):
        super(LedoitWolf, self).__init__(store_covariance, store_precision)


    def fit(self, X):
        covariance, shrinkage = ledoit_wolf(X)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)
        
        return self


