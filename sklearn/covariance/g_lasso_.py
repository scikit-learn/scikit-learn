"""GLasso: sparse inverse covariance estimation with an l1-penalized
estimator.
"""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD Style
# Copyright: INRIA
import warnings

import numpy as np
from scipy import linalg

from .empirical_covariance_ import empirical_covariance, \
                EmpiricalCovariance, log_likelihood

from ..linear_model import lars_path
from ..linear_model import cd_fast
from ..cross_validation import check_cv
from ..externals.joblib import Parallel, delayed


################################################################################
# Helper functions to compute the objective and dual objective functions
# of the l1-penalized estimator
def _objective(mle, precision_, alpha):
    cost = (-log_likelihood(mle, precision_)
            + alpha*(np.abs(precision_).sum()
                        - np.abs(np.diag(precision_)).sum()))
    return cost


def _dual_gap(emp_cov, precision_, alpha):
    """ Expression of the dual gap given in Duchi "Projected Subgradient
    Methods for Learning Sparse Gaussians"
    """
    gap = np.sum(emp_cov * precision_)
    gap -= precision_.shape[0]
    gap += alpha*(np.abs(precision_).sum()
                        - np.abs(np.diag(precision_)).sum())
    return gap


################################################################################
# The g-lasso algorithm

def g_lasso(X, alpha, cov_init=None, mode='cd', tol=1e-4,
            max_iter=100, verbose=False, return_costs=False):
    """ l1-penalized covariance estimator

    Parameters
    ----------
    X: 2D ndarray, shape (n_samples, n_features)
        Data from which to compute the covariance estimate

    alpha: positive float
        The regularization parameter: the higher alpha, the more
        regularization, the sparser the inverse covariance
    cov_init: 2D array (n_features, n_features), optional
        The initial guess for the covariance
    mode: {'cd', 'lars'}
        The Lasso solver to use: coordinate descent or LARS. Use LARS for
        very sparse underlying graphs, where p > n. Elsewhere prefer cd
        which is more numerically stable.
    tol: positive float, optional
        The tolerance to declare convergence: if the dual gap goes below
        this value, iterations are stopped
    max_iter: integer, optional
        The maximum number of iterations
    verbose: boolean, optional
        If verbose is True, the objective function and dual gap are
        printed at each iteration
    return_costs: boolean, optional
        If return_costs is True, the objective function and dual gap
        at each iteration are returned

    Returns
    -------
    covariance_: 2D ndarray, shape (n_features, n_features)
        The estimated covariance matrix
    precision_: 2D ndarray, shape (n_features, n_features)
        The estimated (sparse) precision matrix
    costs: list of (objective, dual_gap) pairs
        The list of values of the objective function and the dual gap at
        each iteration. Returned only if return_costs is True

    """
    _, n_features = X.shape
    mle = empirical_covariance(X)
    if alpha == 0:
        return mle, linalg.inv(mle)
    if cov_init is None:
        covariance_ = mle.copy()
    else:
        covariance_ = cov_init.copy()
        covariance_.flat[::n_features + 1] = mle.flat[::n_features + 1]
    indices = np.arange(n_features)
    precision_ = linalg.inv(covariance_)
    costs = list()
    for i in xrange(max_iter):
        for idx in xrange(n_features):
            sub_covariance = covariance_[indices != idx].T[indices != idx]
            row = mle[idx, indices != idx]
            if mode == 'cd':
                # Use coordinate descent
                coefs = -precision_[indices != idx, idx]/precision_[idx, idx]
                coefs, _, _ = cd_fast.enet_coordinate_descent_gram(coefs,
                                            alpha, 0, sub_covariance,
                                            row, row, max_iter, tol)
            else:
                # Use LARS
                _, _, coefs = lars_path(sub_covariance, row,
                                        Xy=row, Gram=sub_covariance,
                                        alpha_min=alpha/(n_features-1),
                                        copy_Gram=True,
                                        method='lars')
                coefs = coefs[:, -1]
            # Update the precision matrix
            precision_[idx, idx] = 1./(covariance_[idx, idx] -
                        np.dot(covariance_[indices != idx, idx], coefs))
            precision_[indices != idx, idx] = -precision_[idx, idx]*coefs
            precision_[idx, indices != idx] = -precision_[idx, idx]*coefs
            coefs = np.dot(sub_covariance, coefs)
            covariance_[idx, indices != idx] = coefs
            covariance_[indices != idx, idx] = coefs
        d_gap = _dual_gap(mle, precision_, alpha)
        if verbose or return_costs:
            cost = _objective(mle, precision_, alpha)
            if verbose:
                print '[g_lasso] Iteration % 3i, cost %.4f, dual gap %.3e' % (
                                                    i, cost, d_gap)
            if return_costs:
                costs.append((cost, d_gap))
        if np.abs(d_gap) < tol:
            break
    else:
        warnings.warn('g_lasso: did not converge after %i iteration:'
                        'dual gap: %.3e' % (max_iter, d_gap))
    if return_costs:
        return covariance_, precision_, costs
    return covariance_, precision_


class GLasso(EmpiricalCovariance):
    """GLasso: sparse inverse covariance estimation with an l1-penalized
    estimator.

    Attributes
    ----------
    `covariance_` : array-like, shape (n_features, n_features)
        Estimated covariance matrix

    `precision_` : array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.

    """

    def __init__(self, alpha=.01, mode='cd', tol=1e-4,
                 max_iter=100, verbose=False):
        """ l1-penalized covariance estimator

        Parameters
        ----------
        alpha: positive float, optional
            The regularization parameter: the higher alpha, the more
            regularization, the sparser the inverse covariance
        cov_init: 2D array (n_features, n_features), optional
            The initial guess for the covariance
        mode: {'cd', 'lars'}
            The Lasso solver to use: coordinate descent or LARS. Use LARS for
            very sparse underlying graphs, where p > n. Elsewhere prefer cd
            which is more numerically stable.
        tol: positive float, optional
            The tolerance to declare convergence: if the dual gap goes below
            this value, iterations are stopped
        max_iter: integer, optional
            The maximum number of iterations
        verbose: boolean, optional
            If verbose is True, the objective function and dual gap are
            plotted at each iteration

    """
        self.alpha = alpha
        self.mode = mode
        self.tol = tol
        self.max_iter = max_iter
        self.verbose = verbose

    def fit(self, X, y=None):
        self.covariance_, self.precision_ = g_lasso(X,
                                        alpha=self.alpha, mode=self.mode,
                                        tol=self.tol, max_iter=self.max_iter,
                                        verbose=self.verbose,
                                        )
        return self


################################################################################
# Cross-validation with GLasso
def g_lasso_path(X, alphas, X_test=False, mode='cd', tol=1e-4,
                 max_iter=100, verbose=False):
    """ l1-penalized covariance estimator

    Parameters
    ----------
    X: 2D ndarray, shape (n_samples, n_features)
        Data from which to compute the covariance estimate
    alphas: list of positive floats
        The list of regularization parameters, decreasing order
    X_test: 2D array, shape (n_test_samples, n_features), optional
        Optional test matrix to measure generalisation error
    mode: {'cd', 'lars'}
        The Lasso solver to use: coordinate descent or LARS. Use LARS for
        very sparse underlying graphs, where p > n. Elsewhere prefer cd
        which is more numerically stable.
    tol: positive float, optional
        The tolerance to declare convergence: if the dual gap goes below
        this value, iterations are stopped
    max_iter: integer, optional
        The maximum number of iterations
    verbose: boolean, optional
        If verbose is True, the objective function and dual gap are
        printed at each iteration

    Returns
    -------
    covariances_: List of 2D ndarray, shape (n_features, n_features)
        The estimated covariance matrices
    precisions_: List of 2D ndarray, shape (n_features, n_features)
        The estimated (sparse) precision matrices
    scores_: List of float
        The generalisation error (log-likelihood) on the test data.
        Returned only if test data is passed.
    """
    covariance_ = empirical_covariance(X)
    covariances_ = list()
    precisions_ = list()
    scores_ = list()
    if X_test is not None:
        test_emp_cov = empirical_covariance(X_test)
    for alpha in alphas:
        covariance_, precision_ = g_lasso(X, alpha=alpha,
                                    cov_init=covariance_, mode=mode, tol=tol,
                                    max_iter=max_iter, verbose=verbose)
        covariances_.append(covariance_)
        precisions_.append(precision_)
        if X_test is not None:
            scores_.append(log_likelihood(test_emp_cov, precision_))
    if X_test is not None:
        return covariances_, precisions_, scores_
    return covariances_, precisions_


DEFAULT_ALPHAS = np.r_[0, np.logspace(-2.5, .5, 20)][::-1]


class GLassoCV(GLasso):

    def __init__(self, alphas=DEFAULT_ALPHAS, cv=None, tol=1e-4,
                 max_iter=100, mode='cd', n_jobs=1, verbose=False):
        """ l1-penalized covariance estimator

        Parameters
        ----------
        alpha: positive float, optional
            The regularization parameter: the higher alpha, the more
            regularization, the sparser the inverse covariance
        cov_init: 2D array (n_features, n_features), optional
            The initial guess for the covariance
        mode: {'cd', 'lars'}
            The Lasso solver to use: coordinate descent or LARS. Use LARS for
            very sparse underlying graphs, where p > n. Elsewhere prefer cd
            which is more numerically stable.
        tol: positive float, optional
            The tolerance to declare convergence: if the dual gap goes below
            this value, iterations are stopped
        max_iter: integer, optional
            The maximum number of iterations
        cv : crossvalidation generator, optional
            see sklearn.cross_validation module. If None is passed, default to
            a 5-fold strategy
        verbose: boolean, optional
            If verbose is True, the objective function and dual gap are
            plotted at each iteration

        """
        self.alphas = alphas
        self.mode = mode
        self.tol = tol
        self.max_iter = max_iter
        self.verbose = verbose
        self.cv = cv
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        X = np.asarray(X)

        # init cross-validation generator
        cv = check_cv(self.cv, X, y, classifier=False)

        path = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                        delayed(g_lasso_path)(X[train], alphas=self.alphas,
                                    X_test=X[test], mode=self.mode,
                                    tol=self.tol, max_iter=self.max_iter,
                                    verbose=self.verbose)
                        for train, test in cv)

        path = zip(*path)
        best_score = -np.inf
        all_scores = np.array(path[2]).T
        for alpha, scores in zip(self.alphas, all_scores):
            this_score = np.mean(scores)
            if this_score >= best_score:
                best_score = this_score
                best_alpha = alpha
        # XXX: need to store the scores in a name consistent with the
        # rest of the CV estimators
        self.all_scores = all_scores
        self.alpha_ = best_alpha
        self.covariance_, self.precision_ = g_lasso(X, alpha=best_alpha,
                        mode=self.mode, tol=self.tol, max_iter=self.max_iter,
                        verbose=self.verbose)
        return self
