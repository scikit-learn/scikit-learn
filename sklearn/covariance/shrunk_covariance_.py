"""
Covariance estimators using shrinkage.

Shrinkage corresponds to regularising `cov` using a convex combination:
shrunk_cov = (1-shrinkage)*cov + shrinkage*structured_estimate.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD Style.

# avoid division truncation
from __future__ import division
import warnings
import numpy as np

from .empirical_covariance_ import empirical_covariance, EmpiricalCovariance
from ..utils import array2d


###############################################################################
# ShrunkCovariance estimator

def shrunk_covariance(emp_cov, shrinkage=0.1):
    """Calculates a covariance matrix shrunk on the diagonal

    Parameters
    ----------
    emp_cov: array-like, shape (n_features, n_features)
      Covariance matrix to be shrunk

    shrinkage: float, 0 <= shrinkage <= 1
      coefficient in the convex combination used for the computation
      of the shrunk estimate.

    Returns
    -------
    shrunk_cov: array-like
      shrunk covariance

    Notes
    -----
    The regularized (shrunk) covariance is given by

    (1 - shrinkage)*cov
      + shrinkage*mu*np.identity(n_features)

    where mu = trace(cov) / n_features

    """
    emp_cov = array2d(emp_cov)
    n_features = emp_cov.shape[0]

    mu = np.trace(emp_cov) / n_features
    shrunk_cov = (1. - shrinkage) * emp_cov
    shrunk_cov.flat[::n_features + 1] += shrinkage * mu

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
    `covariance_` : array-like, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : array-like, shape (n_features, n_features)
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

    def fit(self, X, assume_centered=False):
        """ Fits the shrunk covariance model
        according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
          Training data, where n_samples is the number of samples
          and n_features is the number of features.

        assume_centered: Boolean
          If True, data are not centered before computation.
          Useful to work with data whose mean is significantly equal to
          zero but is not exactly zero.
          If False, data are centered before computation.

        Returns
        -------
        self : object
            Returns self.

        """
        empirical_cov = empirical_covariance(
            X, assume_centered=assume_centered)
        covariance = shrunk_covariance(empirical_cov, self.shrinkage)
        self._set_estimates(covariance)

        return self


###############################################################################
# Ledoit-Wolf estimator

def ledoit_wolf_shrinkage(X, assume_centered=False, block_size=1000):
    """Estimates the shrunk Ledoit-Wolf covariance matrix.

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data from which to compute the Ledoit-Wolf shrunk covariance shrinkage

    assume_centered: Boolean
      If True, data are not centered before computation.
      Useful to work with data whose mean is significantly equal to
      zero but is not exactly zero.
      If False, data are centered before computation.

    block_size: int,
      Size of the blocks into which the covariance matrix will be split.

    Returns
    -------
    shrinkage: float
      Coefficient in the convex combination used for the computation
      of the shrunk estimate.

    Notes
    -----
    The regularized (shrunk) covariance is:

    (1 - shrinkage)*cov
      + shrinkage * mu * np.identity(n_features)

    where mu = trace(cov) / n_features

    """
    X = np.asarray(X)
    # for only one feature, the result is the same whatever the shrinkage
    if len(X.shape) == 2 and X.shape[1] == 1:
        return 0.
    if X.ndim == 1:
        X = np.reshape(X, (1, -1))
        warnings.warn("Only one sample available. " \
                          "You may want to reshape your data array")
        n_samples = 1
        n_features = X.size
    else:
        n_samples, n_features = X.shape

    # optionaly center data
    if not assume_centered:
        X = X - X.mean(0)

    # number of blocks to split the covariance matrix into
    n_splits = int(n_features / block_size)
    emp_cov_trace = np.sum(X ** 2, 0) / n_samples
    mu = np.sum(emp_cov_trace) / n_features
    X2 = X ** 2
    beta_ = 0.  # sum of the coefficients of <X2.T, X2>
    delta_ = 0.  # sum of the *squared* coefficients of <X.T, X>
    # starting block computation
    for i in xrange(n_splits):
        for j in xrange(n_splits):
            rows = np.arange(block_size * i, block_size * (i + 1))
            cols = np.arange(block_size * j, block_size * (j + 1))
            beta_ += np.sum(np.dot(X2.T[rows], X2[:, cols]))
            delta_ += np.sum(np.dot(X.T[rows], X[:, cols]) ** 2)
        rows = np.arange(block_size * i, block_size * (i + 1))
        beta_ += np.sum(np.dot(X2.T[rows], X2[:, block_size * n_splits:]))
        delta_ += np.sum(
            np.dot(X.T[rows], X[:, block_size * n_splits:]) ** 2)
    for j in xrange(n_splits):
        cols = np.arange(block_size * j, block_size * (j + 1))
        beta_ += np.sum(np.dot(X2.T[block_size * n_splits:], X2[:, cols]))
        delta_ += np.sum(
            np.dot(X.T[block_size * n_splits:], X[:, cols]) ** 2)
    delta_ += np.sum(np.dot(X.T[block_size * n_splits:],
                            X[:, block_size * n_splits:]) ** 2)
    delta_ /= n_samples ** 2
    beta_ += np.sum(np.dot(
            X2.T[block_size * n_splits:], X2[:, block_size * n_splits:]))
    # use delta_ to compute beta
    beta = 1. / (n_features * n_samples) * (beta_ / n_samples - delta_)
    # delta is the sum of the squared coefficients of (<X.T,X> - mu*Id) / p
    delta = delta_ - 2. * mu * emp_cov_trace.sum() + n_features * mu ** 2
    delta /= n_features
    # get final beta as the min between beta and delta
    beta = min(beta, delta)
    # finally get shrinkage
    shrinkage = beta / delta
    # only return the shrunk covariance if it is not to big

    return shrinkage


def ledoit_wolf(X, assume_centered=False, block_size=1000):
    """Estimates the shrunk Ledoit-Wolf covariance matrix.

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data from which to compute the covariance estimate

    assume_centered: Boolean
      If True, data are not centered before computation.
      Usefull to work with data whose mean is significantly equal to
      zero but is not exactly zero.
      If False, data are centered before computation.

    block_size: int,
      Size of the blocks into which the covariance matrix will be split.
      If n_features > `block_size`, an error will be raised since the
      shrunk covariance matrix will be considered as too large regarding
      the available memory.

    Returns
    -------
    shrunk_cov: array-like, shape (n_features, n_features)
      Shrunk covariance.

    shrinkage: float
      Coefficient in the convex combination used for the computation
      of the shrunk estimate.

    Notes
    -----
    The regularized (shrunk) covariance is:

    (1 - shrinkage)*cov
      + shrinkage * mu * np.identity(n_features)

    where mu = trace(cov) / n_features

    """
    X = np.asarray(X)
    # for only one feature, the result is the same whatever the shrinkage
    if len(X.shape) == 2 and X.shape[1] == 1:
        if not assume_centered:
            X = X - X.mean()
        return np.atleast_2d((X ** 2).mean()), 0.
    if X.ndim == 1:
        X = np.reshape(X, (1, -1))
        warnings.warn("Only one sample available. " \
                          "You may want to reshape your data array")
        n_samples = 1
        n_features = X.size
    else:
        n_samples, n_features = X.shape

    if n_features > block_size:
        raise MemoryError("LW: n_features is too large, " +
                          "try increasing block_size")

    # get Ledoit-Wolf shrinkage
    shrinkage = ledoit_wolf_shrinkage(
        X, assume_centered=assume_centered, block_size=block_size)
    emp_cov = empirical_covariance(X, assume_centered=assume_centered)
    mu = np.sum(np.trace(emp_cov)) / n_features
    shrunk_cov = (1. - shrinkage) * emp_cov
    shrunk_cov.flat[::n_features + 1] += shrinkage * mu

    return shrunk_cov, shrinkage


class LedoitWolf(EmpiricalCovariance):
    """LedoitWolf Estimator

    Ledoit-Wolf is a particular form of shrinkage, where the shrinkage
    coefficient is computed using O. Ledoit and M. Wolf's formula as
    described in "A Well-Conditioned Estimator for Large-Dimensional
    Covariance Matrices", Ledoit and Wolf, Journal of Multivariate
    Analysis, Volume 88, Issue 2, February 2004, pages 365-411.

    Parameters
    ----------
    store_precision : bool
        Specify if the estimated precision is stored
    block_size: int,
      Size of the blocks into which the covariance matrix will be split
      during its Ledoit-Wolf estimation.
      If n_features > `block_size`, an error will be raised since the
      shrunk covariance matrix will be considered as too large regarding
      the available memory.

    Attributes
    ----------
    `covariance_` : array-like, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : array-like, shape (n_features, n_features)
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
    and shinkage is given by the Ledoit and Wolf formula (see References)

    References
    ----------
    "A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices",
    Ledoit and Wolf, Journal of Multivariate Analysis, Volume 88, Issue 2,
    February 2004, pages 365-411.

    """
    def __init__(self, store_precision=True, block_size=1000):
        EmpiricalCovariance.__init__(self, store_precision=store_precision)
        self.block_size = block_size

    def fit(self, X, assume_centered=False):
        """ Fits the Ledoit-Wolf shrunk covariance model
        according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
          Training data, where n_samples is the number of samples
          and n_features is the number of features.

        assume_centered: Boolean
          If True, data are not centered before computation.
          Useful to work with data whose mean is significantly equal to
          zero but is not exactly zero.
          If False, data are centered before computation.

        Returns
        -------
        self : object
            Returns self.

        """
        # only return the shrunk covariance if it is not too big
        covariance, shrinkage = ledoit_wolf(
            X, assume_centered=assume_centered, block_size=self.block_size)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)

        return self


###############################################################################
# OAS estimator

def oas(X, assume_centered=False):
    """Estimate covariance with the Oracle Approximating Shrinkage algorithm.

    Parameters
    ----------
    X: array-like, shape (n_samples, n_features)
      Data from which to compute the covariance estimate

    assume_centered: boolean
      If True, data are not centered before computation.
      Useful to work with data whose mean is significantly equal to
      zero but is not exactly zero.
      If False, data are centered before computation.

    Returns
    -------
    shrunk_cov: array-like, shape (n_features, n_features)
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
    X = np.asarray(X)
    # for only one feature, the result is the same whatever the shrinkage
    if len(X.shape) == 2 and X.shape[1] == 1:
        if not assume_centered:
            X = X - X.mean()
        return np.atleast_2d((X ** 2).mean()), 0.
    if X.ndim == 1:
        X = np.reshape(X, (1, -1))
        warnings.warn("Only one sample available. " \
                          "You may want to reshape your data array")
        n_samples = 1
        n_features = X.size
    else:
        n_samples, n_features = X.shape

    emp_cov = empirical_covariance(X, assume_centered=assume_centered)
    mu = np.trace(emp_cov) / n_features

    # formula from Chen et al.'s **implementation**
    alpha = np.mean(emp_cov ** 2)
    num = alpha + mu ** 2
    den = (n_samples + 1.) * (alpha - (mu ** 2) / n_features)

    shrinkage = min(num / den, 1.)
    shrunk_cov = (1. - shrinkage) * emp_cov
    shrunk_cov.flat[::n_features + 1] += shrinkage * mu

    return shrunk_cov, shrinkage


class OAS(EmpiricalCovariance):
    """
    Oracle Approximating Shrinkage Estimator

    OAS is a particular form of shrinkage described in
    "Shrinkage Algorithms for MMSE Covariance Estimation"
    Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.

    The formula used here does not correspond to the one given in the
    article. It has been taken from the Matlab program available from the
    authors' webpage (https://tbayes.eecs.umich.edu/yilun/covestimation).

    Parameters
    ----------
    store_precision : bool
        Specify if the estimated precision is stored

    Attributes
    ----------
    `covariance_` : array-like, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : array-like, shape (n_features, n_features)
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
    and shinkage is given by the OAS formula (see References)

    References
    ----------
    "Shrinkage Algorithms for MMSE Covariance Estimation"
    Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.

    """
    def fit(self, X, assume_centered=False):
        """ Fits the Oracle Approximating Shrinkage covariance model
        according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
          Training data, where n_samples is the number of samples
          and n_features is the number of features.

        assume_centered: boolean
          If True, data are not centered before computation.
          Useful to work with data whose mean is significantly equal to
          zero but is not exactly zero.
          If False, data are centered before computation.

        Returns
        -------
        self : object
            Returns self.

        """
        covariance, shrinkage = oas(X, assume_centered=assume_centered)
        self.shrinkage_ = shrinkage
        self._set_estimates(covariance)

        return self
