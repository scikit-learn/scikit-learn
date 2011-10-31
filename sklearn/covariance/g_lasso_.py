"""GLasso: sparse inverse covariance estimation with an l1-penalized
estimator.
"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

import numpy as np
from scipy import linalg

from .empirical_covariance_ import empirical_covariance, \
                EmpiricalCovariance, log_likelihood
from ..utils.extmath import fast_logdet

from ..linear_model import lars_path

################################################################################
# Helper functions to compute the objective and dual objective functions
# of the l1-penalized estimator
def _objective(mle, precision_, alpha):
    cost = (-log_likelihood(mle, precision_)
            + alpha*np.abs(precision_).sum())
    return cost


def _dual_objective(mle, covariance_, alpha):
    dual_var = covariance_ - mle
    dual_var.flat[::dual_var.shape[0]+1] = 0
    dual_var /= np.maximum(np.abs(dual_var)/alpha, 1)
    B = mle + dual_var
    # It might be necessary to enforce B to be symetric here, if it is
    # not garantied by the estimator. The glasso estimator garanties
    # this.
    #B = B + B.T
    #B *= .5
    return fast_logdet(B) + covariance_.shape[0]


def _dual_gap(emp_cov, precision_, alpha):
    """ Expression of the dual gap given in Duchi "Projected Subgradient
    Methods for Learning Sparse Gaussians"
    """
    gap = np.sum(emp_cov * precision_)
    gap -= precision_.shape[0]
    gap += alpha*np.abs(precision_).sum()
           #             - np.abs(np.diag(precision_)).sum())
    return gap


################################################################################
# The g-lasso algorithm

def g_lasso(X, alpha, tol=1e-4, maxiter=100,
                          eps=10*np.finfo(np.float).eps,
            ):
    # XXX: need to be able to give an initial guess
    _, n_features = X.shape
    covariance_ = empirical_covariance(X)
    if alpha == 0:
        return covariance_
    mle = covariance_.copy()
    covariance_.flat[::n_features + 1] += alpha
    indices = np.arange(n_features)
    precision_ = linalg.inv(covariance_)
    for i in xrange(maxiter):
        for idx in xrange(n_features):
            sub_covariance = covariance_[indices != idx].T[indices != idx]
            row = mle[idx, indices != idx]
            _, _, coefs = lars_path(sub_covariance, row,
                                    Xy=row, Gram=sub_covariance,
                                    alpha_min=alpha/(n_features-1),
                                    copy_Gram=True,
                                    method='lars')
            coefs = coefs[:, -1]
            precision_[idx, idx] = 1./(covariance_[idx, idx] -
                        np.dot(covariance_[indices != idx, idx], coefs))
            precision_[indices != idx, idx] = -precision_[idx, idx]*coefs
            precision_[idx, indices != idx] = -precision_[idx, idx]*coefs
            # We need to regularize the coefs in case the system is
            # ill-conditioned
            #coefs[coefs != 0] += eps
            # update the precision matrix
            coefs = np.dot(sub_covariance, coefs)
            #print 'MMMMM', coefs.max()
            covariance_[idx, indices != idx] = coefs
            covariance_[indices != idx, idx] = coefs
        cost = _objective(mle, precision_, alpha)
        d_gap = _dual_gap(mle, precision_, alpha)
        print cost, d_gap
        if np.abs(d_gap) < tol:
            break
        #_dual_objective(mle, covariance_, alpha), \
    return covariance_, precision_


class GLasso(EmpiricalCovariance):
    """GLasso: sparse inverse covariance estimation with an l1-penalized
    estimator.


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


    """

    def __init__(self, alpha=.01, tol=1e-4):
        self.alpha = alpha
        self.tol = tol
        # XXX: verbosity

    def fit(self, X, y=None):
        self.covariance_ = g_lasso(X, alpha=self.alpha, tol=self.tol)
        return self
