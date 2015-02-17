"""Gaussian processes regression. """

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve
from scipy.optimize import fmin_l_bfgs_b

from sklearn.base import BaseEstimator


class GaussianProcessRegression(BaseEstimator):
    """ Gaussian process regression (GPR).

    The implementation is based on Algorithm 2.1 of ``Gaussian Processes
    for Machine Learning'' (GPML) by Rasmussen and Williams.

    In addition to standard sklearn estimators, GaussianProcessRegression
       * allows prediction without prior fitting (based on the GP prior)
       * provides an additional method sample(X), which evaluates samples drawn
         from the GPR (prior or posterior) at given inputs
       * exposes a method log_marginal_likelihood(theta), which can be used
         externally for other ways of selecting hyperparamters, e.g., via
         Markov chain Monte Carlo.

    Parameters
    ----------
    kernel : Kernel object
        The kernel specifying the covariance function of the GP.

    y_err : float, optional (default: 1e-10)
        Value added to the diagonal of the kernel matrix during fitting.
        Larger values correspond to increased noise level in the observations
        and reduce potential numerical issue during fitting.

    Attributes
    ----------
    X_fit_:

    y_fit_:

    theta_:

    L_:

    alpha_:
    """

    def __init__(self, kernel, y_err=1e-10):
        self.kernel = kernel
        self.y_err = y_err

    def fit(self, X, y):
        self.X_fit_ = np.asarray(X)
        self.y_fit_ = np.asarray(y)

        if self.kernel.has_bounds:
            # Choose hyperparameters based on maximizing the log-marginal
            # likelihood
            def obj_func(theta):
                lml, grad = self.log_marginal_likelihood(theta,
                                                         eval_gradient=True)
                return -lml, -grad
            self.theta_, lml, _ = fmin_l_bfgs_b(obj_func, self.kernel.params,
                                                bounds=self.kernel.bounds)
            self.kernel.params = self.theta_
        else:
            self.theta_ = self.kernel.params

        # Precompute quantities required for predictions which are independent
        # of actual query points
        K = self.kernel.auto_correlation(self.X_fit_)
        K[np.diag_indices_from(K)] += self.y_err
        self.L_ = cholesky(K, lower=True)  # Line 2
        self.alpha_ = cho_solve((self.L_, True), self.y_fit_)  # Line 3

        return self

    def predict(self, X, return_cov=False):
        X = np.asarray(X)

        if not hasattr(self, "X_fit_"):  # Unfitted; predict based on GP prior
            y_mean = np.zeros(X.shape[0])
            if return_cov:
                y_cov = self.kernel.auto_correlation(X)
                return y_mean, y_cov
            else:
                return y_mean
        else:  # Predict based on GP posterior
            K_trans = self.kernel.cross_correlation(X, self.X_fit_)
            y_mean = K_trans.dot(self.alpha_)  # Line 4 (y_mean = f_star)
            if return_cov:
                v = cho_solve((self.L_, True), K_trans.T)  # Line 5
                y_cov = \
                    self.kernel.auto_correlation(X) - K_trans.dot(v)  # Line 6
                return y_mean, y_cov
            else:
                return y_mean

    def sample(self, X, n_samples=1):
        y_mean, y_cov = self.predict(X, return_cov=True)
        y_samples = \
            np.random.multivariate_normal(y_mean, y_cov, n_samples).T
        return y_samples

    def log_marginal_likelihood(self, theta, eval_gradient=False):
        import copy  # XXX: Avoid deepcopy
        kernel = copy.deepcopy(self.kernel)
        kernel.params = theta

        if eval_gradient:
            K, K_gradient = \
                kernel.auto_correlation(self.X_fit_, eval_gradient=True)
        else:
            K = kernel.auto_correlation(self.X_fit_)

        K[np.diag_indices_from(K)] += self.y_err
        L = cholesky(K, lower=True)  # Line 2
        alpha = cho_solve((L, True), self.y_fit_)  # Line 3

        # Compute log-likelihood (compare line 7)
        log_likelihood = -0.5*self.y_fit_.dot(alpha)
        log_likelihood -= np.log(np.diag(L)).sum()
        log_likelihood -= K.shape[0] / 2 * np.log(2 * np.pi)

        if eval_gradient:  # compare Equation 5.9 from GPML
            tmp = np.outer(alpha, alpha)
            tmp -= cho_solve((L, True), np.eye(K.shape[0]))
            gradient = 0.5 * np.trace(tmp.dot(K_gradient))
            return log_likelihood, gradient
        else:
            return log_likelihood
