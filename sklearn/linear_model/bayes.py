"""
Various bayesian regression
"""
from __future__ import print_function

# Authors: V. Michel, F. Pedregosa, A. Gramfort
# License: BSD 3 clause

from math import log
import numpy as np
from scipy import linalg

from .base import LinearModel
from ..base import RegressorMixin
from ..utils.extmath import fast_logdet, pinvh
from ..utils import check_X_y


###############################################################################
# BayesianRidge regression

class BayesianRidge(LinearModel, RegressorMixin):
    """Bayesian ridge regression

    Fit a Bayesian ridge model and optimize the regularization parameters
    lambda (precision of the weights) and alpha (precision of the noise).

    Parameters
    ----------
    n_iter : int, optional
        Maximum number of iterations.  Default is 300.

    tol : float, optional
        Stop the algorithm if w has converged. Default is 1.e-3.

    alpha_1 : float, optional
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the alpha parameter. Default is 1.e-6

    alpha_2 : float, optional
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the alpha parameter.
        Default is 1.e-6.

    lambda_1 : float, optional
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the lambda parameter. Default is 1.e-6.

    lambda_2 : float, optional
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the lambda parameter.
        Default is 1.e-6

    compute_score : boolean, optional
        If True, compute the objective function at each step of the model.
        Default is False

    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).
        Default is True.

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    verbose : boolean, optional, default False
        Verbose mode when fitting the model.


    Attributes
    ----------
    coef_ : array, shape = (n_features)
        Coefficients of the regression model (mean of distribution)

    alpha_ : float
       estimated precision of the noise.

    lambda_ : array, shape = (n_features)
       estimated precisions of the weights.

    scores_ : float
        if computed, value of the objective function (to be maximized)

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.BayesianRidge()
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    ... # doctest: +NORMALIZE_WHITESPACE
    BayesianRidge(alpha_1=1e-06, alpha_2=1e-06, compute_score=False,
            copy_X=True, fit_intercept=True, lambda_1=1e-06, lambda_2=1e-06,
            n_iter=300, normalize=False, tol=0.001, verbose=False)
    >>> clf.predict([[1, 1]])
    array([ 1.])

    Notes
    -----
    See examples/linear_model/plot_bayesian_ridge.py for an example.
    """

    def __init__(self, n_iter=300, tol=1.e-3, alpha_1=1.e-6, alpha_2=1.e-6,
                 lambda_1=1.e-6, lambda_2=1.e-6, compute_score=False,
                 fit_intercept=True, normalize=False, copy_X=True,
                 verbose=False):
        self.n_iter = n_iter
        self.tol = tol
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.compute_score = compute_score
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X
        self.verbose = verbose

    def fit(self, X, y):
        """Fit the model

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(X, y, dtype=np.float)
        X, y, X_mean, y_mean, X_std = self._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X)
        n_samples, n_features = X.shape

        ### Initialization of the values of the parameters
        alpha_ = 1. / np.var(y)
        lambda_ = 1.

        verbose = self.verbose
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        alpha_1 = self.alpha_1
        alpha_2 = self.alpha_2

        self.scores_ = list()
        coef_old_ = None

        XT_y = np.dot(X.T, y)
        U, S, Vh = linalg.svd(X, full_matrices=False)
        eigen_vals_ = S ** 2

        ### Convergence loop of the bayesian ridge regression
        for iter_ in range(self.n_iter):

            ### Compute mu and sigma
            # sigma_ = lambda_ / alpha_ * np.eye(n_features) + np.dot(X.T, X)
            # coef_ = sigma_^-1 * XT * y
            if n_samples > n_features:
                coef_ = np.dot(Vh.T,
                               Vh / (eigen_vals_ + lambda_ / alpha_)[:, None])
                coef_ = np.dot(coef_, XT_y)
                if self.compute_score:
                    logdet_sigma_ = - np.sum(
                        np.log(lambda_ + alpha_ * eigen_vals_))
            else:
                coef_ = np.dot(X.T, np.dot(
                    U / (eigen_vals_ + lambda_ / alpha_)[None, :], U.T))
                coef_ = np.dot(coef_, y)
                if self.compute_score:
                    logdet_sigma_ = lambda_ * np.ones(n_features)
                    logdet_sigma_[:n_samples] += alpha_ * eigen_vals_
                    logdet_sigma_ = - np.sum(np.log(logdet_sigma_))

            ### Update alpha and lambda
            rmse_ = np.sum((y - np.dot(X, coef_)) ** 2)
            gamma_ = (np.sum((alpha_ * eigen_vals_)
                      / (lambda_ + alpha_ * eigen_vals_)))
            lambda_ = ((gamma_ + 2 * lambda_1)
                       / (np.sum(coef_ ** 2) + 2 * lambda_2))
            alpha_ = ((n_samples - gamma_ + 2 * alpha_1)
                      / (rmse_ + 2 * alpha_2))

            ### Compute the objective function
            if self.compute_score:
                s = lambda_1 * log(lambda_) - lambda_2 * lambda_
                s += alpha_1 * log(alpha_) - alpha_2 * alpha_
                s += 0.5 * (n_features * log(lambda_)
                            + n_samples * log(alpha_)
                            - alpha_ * rmse_
                            - (lambda_ * np.sum(coef_ ** 2))
                            - logdet_sigma_
                            - n_samples * log(2 * np.pi))
                self.scores_.append(s)

            ### Check for convergence
            if iter_ != 0 and np.sum(np.abs(coef_old_ - coef_)) < self.tol:
                if verbose:
                    print("Convergence after ", str(iter_), " iterations")
                break
            coef_old_ = np.copy(coef_)

        self.alpha_ = alpha_
        self.lambda_ = lambda_
        self.coef_ = coef_

        self._set_intercept(X_mean, y_mean, X_std)
        return self


###############################################################################
# ARD (Automatic Relevance Determination) regression


class ARDRegression(LinearModel, RegressorMixin):
    """Bayesian ARD regression.

    Fit the weights of a regression model, using an ARD prior. The weights of
    the regression model are assumed to be in Gaussian distributions.
    Also estimate the parameters lambda (precisions of the distributions of the
    weights) and alpha (precision of the distribution of the noise).
    The estimation is done by an iterative procedures (Evidence Maximization)

    Parameters
    ----------
    n_iter : int, optional
        Maximum number of iterations. Default is 300

    tol : float, optional
        Stop the algorithm if w has converged. Default is 1.e-3.

    alpha_1 : float, optional
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the alpha parameter. Default is 1.e-6.

    alpha_2 : float, optional
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the alpha parameter. Default is 1.e-6.

    lambda_1 : float, optional
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the lambda parameter. Default is 1.e-6.

    lambda_2 : float, optional
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the lambda parameter. Default is 1.e-6.

    compute_score : boolean, optional
        If True, compute the objective function at each step of the model.
        Default is False.

    threshold_lambda : float, optional
        threshold for removing (pruning) weights with high precision from
        the computation. Default is 1.e+4.

    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).
        Default is True.

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    copy_X : boolean, optional, default True.
        If True, X will be copied; else, it may be overwritten.

    verbose : boolean, optional, default False
        Verbose mode when fitting the model.

    Attributes
    ----------
    coef_ : array, shape = (n_features)
        Coefficients of the regression model (mean of distribution)

    alpha_ : float
       estimated precision of the noise.

    lambda_ : array, shape = (n_features)
       estimated precisions of the weights.

    sigma_ : array, shape = (n_features, n_features)
        estimated variance-covariance matrix of the weights

    scores_ : float
        if computed, value of the objective function (to be maximized)

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.ARDRegression()
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    ... # doctest: +NORMALIZE_WHITESPACE
    ARDRegression(alpha_1=1e-06, alpha_2=1e-06, compute_score=False,
            copy_X=True, fit_intercept=True, lambda_1=1e-06, lambda_2=1e-06,
            n_iter=300, normalize=False, threshold_lambda=10000.0, tol=0.001,
            verbose=False)
    >>> clf.predict([[1, 1]])
    array([ 1.])

    Notes
    --------
    See examples/linear_model/plot_ard.py for an example.
    """

    def __init__(self, n_iter=300, tol=1.e-3, alpha_1=1.e-6, alpha_2=1.e-6,
                 lambda_1=1.e-6, lambda_2=1.e-6, compute_score=False,
                 threshold_lambda=1.e+4, fit_intercept=True, normalize=False,
                 copy_X=True, verbose=False):
        self.n_iter = n_iter
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.compute_score = compute_score
        self.threshold_lambda = threshold_lambda
        self.copy_X = copy_X
        self.verbose = verbose

    def fit(self, X, y):
        """Fit the ARDRegression model according to the given training data
        and parameters.

        Iterative procedure to maximize the evidence

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers)

        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(X, y, dtype=np.float)

        n_samples, n_features = X.shape
        coef_ = np.zeros(n_features)

        X, y, X_mean, y_mean, X_std = self._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X)

        ### Launch the convergence loop
        keep_lambda = np.ones(n_features, dtype=bool)

        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        alpha_1 = self.alpha_1
        alpha_2 = self.alpha_2
        verbose = self.verbose

        ### Initialization of the values of the parameters
        alpha_ = 1. / np.var(y)
        lambda_ = np.ones(n_features)

        self.scores_ = list()
        coef_old_ = None

        ### Iterative procedure of ARDRegression
        for iter_ in range(self.n_iter):
            ### Compute mu and sigma (using Woodbury matrix identity)
            sigma_ = pinvh(np.eye(n_samples) / alpha_ +
                           np.dot(X[:, keep_lambda] *
                           np.reshape(1. / lambda_[keep_lambda], [1, -1]),
                           X[:, keep_lambda].T))
            sigma_ = np.dot(sigma_, X[:, keep_lambda]
                            * np.reshape(1. / lambda_[keep_lambda], [1, -1]))
            sigma_ = - np.dot(np.reshape(1. / lambda_[keep_lambda], [-1, 1])
                              * X[:, keep_lambda].T, sigma_)
            sigma_.flat[::(sigma_.shape[1] + 1)] += 1. / lambda_[keep_lambda]
            coef_[keep_lambda] = alpha_ * np.dot(
                sigma_, np.dot(X[:, keep_lambda].T, y))

            ### Update alpha and lambda
            rmse_ = np.sum((y - np.dot(X, coef_)) ** 2)
            gamma_ = 1. - lambda_[keep_lambda] * np.diag(sigma_)
            lambda_[keep_lambda] = ((gamma_ + 2. * lambda_1)
                                    / ((coef_[keep_lambda]) ** 2
                                       + 2. * lambda_2))
            alpha_ = ((n_samples - gamma_.sum() + 2. * alpha_1)
                      / (rmse_ + 2. * alpha_2))

            ### Prune the weights with a precision over a threshold
            keep_lambda = lambda_ < self.threshold_lambda
            coef_[~keep_lambda] = 0

            ### Compute the objective function
            if self.compute_score:
                s = (lambda_1 * np.log(lambda_) - lambda_2 * lambda_).sum()
                s += alpha_1 * log(alpha_) - alpha_2 * alpha_
                s += 0.5 * (fast_logdet(sigma_) + n_samples * log(alpha_)
                                                + np.sum(np.log(lambda_)))
                s -= 0.5 * (alpha_ * rmse_ + (lambda_ * coef_ ** 2).sum())
                self.scores_.append(s)

            ### Check for convergence
            if iter_ > 0 and np.sum(np.abs(coef_old_ - coef_)) < self.tol:
                if verbose:
                    print("Converged after %s iterations" % iter_)
                break
            coef_old_ = np.copy(coef_)

        self.coef_ = coef_
        self.alpha_ = alpha_
        self.sigma_ = sigma_
        self.lambda_ = lambda_
        self._set_intercept(X_mean, y_mean, X_std)
        return self
