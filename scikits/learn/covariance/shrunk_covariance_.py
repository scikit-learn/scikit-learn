"""
Covariance estimators using shrinkage.

Shrinkage corresponds to regularising `cov` using a convex combination:
shrunk_cov = (1-shrinkage)*cov + shrinkage*structured_estimate.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

from __future__ import division
import numpy as np

from .empirical_covariance_ import EmpiricalCovariance

###############################################################################
# ShrunkCovariance estimator

def shrunk_covariance(emp_cov, shrinkage=0.1):
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
    
    mu = np.trace(emp_cov) / n_features
    shrunk_cov = (1.-shrinkage)*emp_cov
    shrunk_cov.flat[::n_features + 1] += shrinkage*mu

    return shrunk_cov

class ShrunkCovariance(EmpiricalCovariance):
    """Covariance estimator with shrinkage
    
    Parameters
    ----------
    store_precision : bool
      Specify if the estimated precision is stored
      
    shrinkage: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)
        
    `shrinkage`: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
    Notes
    -----
    The regularized covariance is given by
    
    (1 - shrinkage)*cov
      + shrinkage*mu*np.identity(n_features)
    
    where mu = trace(cov) / n_features
    
    """
    def __init__(self, store_precision=True, shrinkage=0.1):
        self.store_precision = store_precision
        self.shrinkage = shrinkage


    def fit(self, X, assume_centered=False, **params):
        self._set_params(**params)
        if assume_centered:
            empirical_cov = np.dot(X.T, X) / X.shape[0]
        else:
            empirical_cov = np.cov(X.T, bias=1)
        covariance = shrunk_covariance(empirical_cov, self.shrinkage)
        self._set_estimates(covariance)
        
        return self

################################################################################
# Ledoit-Wolf estimator

def ledoit_wolf(X, assume_centered=False):
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
    X = np.asanyarray(X)
    # for only one feature, the result is the same whatever the shrinkage
    if X.ndim == 1:
        if not assume_centered:
            X = X - X.mean()
            print np.atleast_2d((X**2).mean())
        return np.atleast_2d((X**2).mean()), 0.
    n_samples, n_features = X.shape
    
    # optionaly center data
    if not assume_centered:
        emp_cov = np.cov(X.T, bias=1)
        X = X - X.mean(0)
    else:
        emp_cov = np.dot(X.T, X) / n_samples
    
    mu = np.trace(emp_cov) / n_features
    delta_ = emp_cov.copy()
    delta_.flat[::n_features + 1] -= mu
    delta = (delta_**2).sum() / n_features
    X2 = X**2
    beta_ = 1./(n_features*n_samples) \
        * np.sum(np.dot(X2.T, X2)/n_samples - emp_cov**2)
    
    beta = min(beta_, delta)
    shrinkage = beta/delta
    shrunk_cov = (1.-shrinkage)*emp_cov
    shrunk_cov.flat[::n_features + 1] += shrinkage*mu
    
    return shrunk_cov, shrinkage


class LedoitWolf(EmpiricalCovariance):
    """LedoitWolf Estimator
    
    Ledoit-Wolf is a particular form of shrinkage, where the shrinkage
    coefficient is computed using O.Ledoit and M.Wolf's formula as
    described in "A Well-Conditioned Estimator for Large-Dimensional
    Covariance Matrices", Ledoit and Wolf, Journal of Multivariate
    Analysis, Volume 88, Issue 2, February 2004, pages 365-411.
    
    Parameters
    ----------
    store_precision : bool
        Specify if the estimated precision is stored
    
    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    `shrinkage_`: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.

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
    def fit(self, X, assume_centered=False):
        covariance, shrinkage = ledoit_wolf(X, assume_centered=assume_centered)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)
        
        return self

################################################################################
# OAS estimator

def oas(X, assume_centered=False):
    """Estimates the covariance matrix with the Oracle Approximating Shrinkage.

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
    X = np.asanyarray(X)
    # for only one feature, the result is the same whatever the shrinkage
    if X.ndim == 1:
        if not assume_centered:
            X = X - X.mean()
        return np.atleast_2d((X**2).mean()), 0.
    n_samples, n_features = X.shape
    
    if not assume_centered:
        emp_cov = np.cov(X.T, bias=1)
    else:
        emp_cov = np.dot(X.T, X) / n_samples
    mu = np.trace(emp_cov) / n_features
    
    # formula from Chen et al.'s **implementation**
    alpha = np.mean(emp_cov**2)
    num = alpha + mu**2
    den = (n_samples+1.) * (alpha - (mu**2)/n_features)
    
    shrinkage = min(num/den, 1.)
    shrunk_cov = (1.-shrinkage)*emp_cov
    shrunk_cov.flat[::n_features+1] += shrinkage*mu
    
    return shrunk_cov, shrinkage


class OAS(EmpiricalCovariance):
    """
    Oracle Approximating Shrinkage Estimator

    OAS is a particular form of shrinkage described in
    "Shrinkage Algorithms for MMSE Covariance Estimation"
    Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
    
    The formula used here does not correspond to the one given in the
    article. It has been taken from the matlab programm available from the
    authors webpage (https://tbayes.eecs.umich.edu/yilun/covestimation).
    
    Parameters
    ----------
    store_precision : bool
        Specify if the estimated precision is stored
    
    Attributes
    ----------
    `covariance_` : 2D ndarray, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : 2D ndarray, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)
        
    `shrinkage_`: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.
    
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
    def fit(self, X, assume_centered=False):
        covariance, shrinkage = oas(X, assume_centered=assume_centered)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)
        
        return self
