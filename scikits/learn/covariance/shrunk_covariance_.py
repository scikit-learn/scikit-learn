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

from .base_covariance_ import BaseCovariance

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


class ShrunkCovariance(BaseCovariance):
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
        emp_cov = BaseCovariance(store_precision=False)
        emp_cov.fit(X)
        covariance = shrunk_covariance(emp_cov.covariance_, self.shrinkage_)
        self._set_estimates(covariance)
        
        return self

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
    
    emp_cov = BaseCovariance(store_precision=False).fit(X).covariance_
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
    """LedoitWolf Estimator
    
    Ledoit-Wolf is a particular form of shrinkage, where the shrinkage
    coefficient is computed using O.Ledoit and M.Wolf's formula as
    described in "A Well-Conditioned Estimator for Large-Dimensional
    Covariance Matrices", Ledoit and Wolf, Journal of Multivariate
    Analysis, Volume 88, Issue 2, February 2004, pages 365-411.

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
    
    emp_cov = BaseCovariance(store_precision=False).fit(X).covariance_
    structured_estimator = np.eye(n_features)
    mu = np.trace(emp_cov) / n_features
    
    # formula from Chen et al.'s implementation
    num = np.sum(emp_cov**2) + (np.trace(emp_cov)**2)
    den = (n_samples+1.) \
        * (np.sum(emp_cov**2) - ((np.trace(emp_cov)**2)/n_features))
    
    shrinkage = min(num/den, 1.)
    shrunk_cov = (1.-shrinkage)*emp_cov + shrinkage*mu*structured_estimator
    
    return shrunk_cov, shrinkage


class OAS(ShrunkCovariance):
    """
    OAS Estimator

    OAS is a particular form of shrinkage described in
    "Shrinkage Algorithms for MMSE Covariance Estimation"
    Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
    
    The formula used here does not correspond to the one given in the
    article. It has been taken from the matlab programm available from the
    authors webpage (https://tbayes.eecs.umich.edu/yilun/covestimation).
    
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
