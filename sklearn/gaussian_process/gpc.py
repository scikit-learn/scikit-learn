"""Gaussian processes classification based on the Laplace approximation. """

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import erf

from sklearn.base import BaseEstimator
from sklearn.gaussian_process.kernels import RBF
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array


# Values required for approximating the logistic sigmoid by
# error functions. coefs are obtained via:
# x = np.array([0, 0.6, 2, 3.5, 4.5, np.inf])
# b = logistic(x)
# A = (erf(np.dot(x, self.lambdas)) + 1) / 2
# coefs = lstsq(A, b)[0]
LAMBDAS = np.array([0.41, 0.4, 0.37, 0.44, 0.39])[:, np.newaxis]
COEFS = np.array([-1854.8214151, 3516.89893646, 221.29346712,
                  128.12323805, -2010.49422654])[:, np.newaxis]


class GaussianProcessClassification(BaseEstimator):
    """ Gaussian process classification (GPC).

    The implementation is based on Algorithm 3.1, 3.2, and 5.1 of
    ``Gaussian Processes for Machine Learning'' (GPML) by Rasmussen and
    Williams.

    Internally, the Laplace approximation is used for approximating the
    non-Gaussian posterior by a Gaussian.

    Currently, the implementation is restricted to
      * using the logistic link function
      * binary classification
    """

    def __init__(self, kernel=RBF(), jitter=0.0):
        self.kernel = kernel
        self.jitter = jitter

    def fit(self, X, y):
        X, y = check_X_y(X, y)

        # XXX: Assert that y is binary and labels are {0, 1}
        self.X_fit_ = X
        self.y_fit_ = y

        if self.kernel.has_bounds:
            # Choose hyperparameters based on maximizing the log-marginal
            # likelihood
            def obj_func(theta):
                lml, grad = self.log_marginal_likelihood(theta,
                                                         eval_gradient=True)
                return -lml, -grad
            self.theta_, _, _ = fmin_l_bfgs_b(obj_func, self.kernel.params,
                                              bounds=self.kernel.bounds)
            self.kernel.params = self.theta_
        else:
            self.theta_ = self.kernel.params

        # Precompute quantities required for predictions which are independent
        # of actual query points
        self.K_ = self.kernel(self.X_fit_)
        self.K_[np.diag_indices_from(self.K_)] += self.jitter

        self.f_, _, (self.pi, self.W_sr, self.L, _, _) = \
            self._posterior_mode(self.K_, return_temporaries=True)

        return self

    def predict(self, X):
        check_is_fitted(self, ["X_fit_", "y_fit_", "K_", "f_"])
        X = check_array(X)

        # As discussed on Section 3.4.2 of GPML, for making hard binary
        # decisions, it is enough to compute the MAP of the posterior and
        # pass it through the link function
        K_star = \
            self.kernel(self.X_fit_, X)  # K_star =k(x_star)
        f_star = K_star.T.dot(self.y_fit_ - self.pi)  # Line 4 (Algorithm 3.2)

        return f_star > 0

    def predict_proba(self, X):
        check_is_fitted(self, ["X_fit_", "y_fit_", "K_", "f_"])
        X = check_array(X)

        # Based on Algorithm 3.2 of GPML
        K_star = self.kernel(self.X_fit_, X)  # K_star =k(x_star)
        f_star = K_star.T.dot(self.y_fit_ - self.pi)  # Line 4
        v = solve(self.L, self.W_sr.dot(K_star))  # Line 5
        var_f_star = self.kernel(X) - v.T.dot(v)  # Line 6

        # Line 7:
        # Approximate \int log(z) * N(z | f_star, var_f_star)
        # Approximation is due to Williams & Barber, "Bayesian Classification
        # with Gaussian Processes", Appendix A:
        # Approximate the logistic sigmoid by a linear combination of
        # 5 error functions
        # See
        # blitiri.blogspot.de/2012/11/gaussian-integral-of-error-function.html
        # for information on how this integral can be computed
        alpha = 1 / (2 * np.diag(var_f_star))
        gamma = LAMBDAS * f_star
        integrals = np.sqrt(np.pi / alpha) \
            * erf(gamma * np.sqrt(alpha / (alpha + LAMBDAS**2))) \
            / (2 * np.sqrt(np.diag(var_f_star) * 2 * np.pi))
        pi_star = (COEFS * integrals).sum(axis=0) + .5 * COEFS.sum()

        return pi_star

    def log_marginal_likelihood(self, theta, eval_gradient=False):
        import copy  # XXX
        kernel = copy.deepcopy(self.kernel)
        kernel.params = theta

        if eval_gradient:
            K, K_gradient = kernel(self.X_fit_, eval_gradient=True)
        else:
            K = kernel(self.X_fit_)

        K[np.diag_indices_from(K)] += self.jitter

        # Return temporaries
        f, Z, (pi, W_sr, L, b, a) = \
            self._posterior_mode(K, return_temporaries=True)

        if not eval_gradient:
            return Z

        # Compute gradient based on Algorithm 5.1 of GPML
        d_Z = np.empty(theta.shape[0])
        R = W_sr.dot(cho_solve((L, True), W_sr))  # Line 7
        C = solve(L, W_sr.dot(K))  # Line 8
        # Line 9:
        s_2 = -0.5*(np.diag(K) - np.diag(C.T.dot(C))) \
            * (pi * (1 - pi) * (1 - 2*pi))  # third derivative
        for j in range(d_Z.shape[0]):
            C = K_gradient[..., j]   # Line 11
            s_1 = .5 * a.T.dot(C).dot(a) - .5 * np.trace(R.dot(C))  # Line 12

            b = C.dot(self.y_fit_ - pi)  # Line 13
            s_3 = b - K.dot(R).dot(b)  # Line 14

            d_Z[j] = s_1 + s_2.T.dot(s_3)  # Line 15

        return Z, d_Z

    def _posterior_mode(self, K, return_temporaries=False):
        # Based on Algorithm 3.1 of GPML
        f = np.zeros_like(self.y_fit_, dtype=np.float64)
        log_marginal_likelihood = -np.inf
        while True:
            # Line 4
            pi = 1 / (1 + np.exp(-f))
            W = pi * (1 - pi)
            # Line 5
            W_sr = np.diag(np.sqrt(W))  # XXX: avoid creating square matrix?
            B = np.eye(W.shape[0]) + W_sr.dot(K).dot(W_sr)
            L = cholesky(B, lower=True)
            # Line 6
            b = W * f + (self.y_fit_ - pi)
            # Line 7
            a = b - W_sr.dot(cho_solve((L, True), W_sr.dot(K).dot(b)))
            # Line 8
            f = K.dot(a)

            # Line 10: Compute log marginal likelihood in loop and use as
            #          convergence criterion
            lml = -0.5*a.T.dot(f) \
                - np.log(1 + np.exp(-(self.y_fit_*2 - 1)*f)).sum() \
                - np.log(np.diag(L)).sum()
            # Check if we have converged (log marginal likelihood does
            # not decrease)
            # XXX: more complex convergence criterion
            if lml - log_marginal_likelihood < 1e-10:
                break
            log_marginal_likelihood = lml

        if return_temporaries:
            return f, log_marginal_likelihood, (pi, W_sr, L, b, a)
        else:
            return f, log_marginal_likelihood
