"""Gaussian processes regression. """

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve, solve_triangular
from scipy.optimize import fmin_l_bfgs_b

from sklearn.base import BaseEstimator, clone
from sklearn.gaussian_process.kernels import RBF
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_X_y, check_array


class GaussianProcessRegressor(BaseEstimator):
    """ Gaussian process regression (GPR).

    The implementation is based on Algorithm 2.1 of ``Gaussian Processes
    for Machine Learning'' (GPML) by Rasmussen and Williams.

    In addition to standard sklearn estimators, GaussianProcessRegressor
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

    optimizer : string, optional (default: "fmin_l_bfgs_b")
        A string specifying the optimization algorithm used for optimizing the
        kernel's parameters. Default uses 'fmin_l_bfgs_b' algorithm from
        scipy.optimize. If None, the kernel's paramters are kept fixed.
        Available optimizers are::

            'fmin_l_bfgs_b'

    Attributes
    ----------
    X_fit_:

    y_fit_:

    theta_:

    L_:

    alpha_:
    """

    def __init__(self, kernel=None, y_err=1e-10, optimizer="fmin_l_bfgs_b"):
        self.kernel = kernel
        self.y_err = y_err
        self.optimizer = optimizer

    def fit(self, X, y):
        if self.kernel is None:  # Use an RBF kernel as default
            self.kernel_ = RBF()
        else:
            self.kernel_ = clone(self.kernel)

        X, y = check_X_y(X, y)

        self.X_fit_ = X
        self.y_fit_ = y

        if self.optimizer == "fmin_l_bfgs_b":
            # Choose hyperparameters based on maximizing the log-marginal
            # likelihood using fmin_l_bfgs_b
            def obj_func(theta):
                lml, grad = self.log_marginal_likelihood(theta,
                                                         eval_gradient=True)
                return -lml, -grad
            self.theta_, _, _ = fmin_l_bfgs_b(obj_func, self.kernel_.theta,
                                              bounds=self.kernel_.bounds)
            self.kernel_.theta = self.theta_
        elif self.optimizer is None:
            self.theta_ = self.kernel_.theta
        else:
            raise ValueError("Unknown optimizer %s." % self.optimizer)

        # Precompute quantities required for predictions which are independent
        # of actual query points
        K = self.kernel_(self.X_fit_)
        K[np.diag_indices_from(K)] += self.y_err
        self.L_ = cholesky(K, lower=True)  # Line 2
        self.alpha_ = cho_solve((self.L_, True), self.y_fit_)  # Line 3

        return self

    def predict(self, X, return_std=False, return_cov=False):
        if return_std and return_cov:
            raise RuntimeError(
                "Not returning standard deviation of predictions when "
                "returning full covariance.")

        X = check_array(X)

        if not hasattr(self, "X_fit_"):  # Unfitted; predict based on GP prior
            y_mean = np.zeros(X.shape[0])
            if return_cov:
                y_cov = self.kernel(X)
                return y_mean, y_cov
            elif return_std:
                y_var = np.apply_along_axis(self.kernel, 1, X)[:, 0]
                return y_mean, np.sqrt(y_var)
            else:
                return y_mean
        else:  # Predict based on GP posterior
            K_trans = self.kernel_(X, self.X_fit_)
            y_mean = K_trans.dot(self.alpha_)  # Line 4 (y_mean = f_star)
            if return_cov:
                v = cho_solve((self.L_, True), K_trans.T)  # Line 5
                y_cov = self.kernel_(X) - K_trans.dot(v)  # Line 6
                return y_mean, y_cov
            elif return_std:
                # compute inverse K_inv of K based on its cholesky
                # decomposition L and its inverse L_inv
                L_inv = solve_triangular(self.L_.T, np.eye(self.L_.shape[0]))
                K_inv = L_inv.dot(L_inv.T)
                # Compute variance of predictive distribution
                y_var = np.apply_along_axis(self.kernel_, 1, X)[:, 0]
                y_var -= np.sum(K_trans.T[:, np.newaxis] * K_trans.T
                                * K_inv[:, :, np.newaxis],
                                axis=(0, 1))
                return y_mean, np.sqrt(y_var)
            else:
                return y_mean

    def sample_y(self, X, n_samples=1, random_state=0):
        rng = check_random_state(random_state)

        y_mean, y_cov = self.predict(X, return_cov=True)
        y_samples = rng.multivariate_normal(y_mean, y_cov, n_samples).T
        return y_samples

    def log_marginal_likelihood(self, theta, eval_gradient=False):
        kernel = self.kernel_.clone_with_theta(theta)

        if eval_gradient:
            K, K_gradient = kernel(self.X_fit_, eval_gradient=True)
        else:
            K = kernel(self.X_fit_)

        K[np.diag_indices_from(K)] += self.y_err
        try:
            L = cholesky(K, lower=True)  # Line 2
        except np.linalg.LinAlgError:
            return (-np.inf, np.zeros_like(theta))\
                 if eval_gradient else -np.inf

        alpha = cho_solve((L, True), self.y_fit_)  # Line 3

        # Compute log-likelihood (compare line 7)
        log_likelihood = -0.5*self.y_fit_.dot(alpha)
        log_likelihood -= np.log(np.diag(L)).sum()
        log_likelihood -= K.shape[0] / 2 * np.log(2 * np.pi)

        if eval_gradient:  # compare Equation 5.9 from GPML
            tmp = np.outer(alpha, alpha)
            tmp -= cho_solve((L, True), np.eye(K.shape[0]))
            # Compute "0.5 * trace(tmp.dot(K_gradient))" without constructing
            # the full matrix tmp.dot(K_gradient) since only its diagonal is
            # required
            gradient = 0.5 * np.einsum("ij,ijk->k", tmp, K_gradient)
            return log_likelihood, gradient
        else:
            return log_likelihood
