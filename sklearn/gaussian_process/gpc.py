"""Gaussian processes classification."""

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import warnings
from operator import itemgetter

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import erf

from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.gaussian_process.kernels import RBF
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array
from sklearn.utils import check_random_state
from sklearn.preprocessing import LabelEncoder

# Values required for approximating the logistic sigmoid by
# error functions. coefs are obtained via:
# x = np.array([0, 0.6, 2, 3.5, 4.5, np.inf])
# b = logistic(x)
# A = (erf(np.dot(x, self.lambdas)) + 1) / 2
# coefs = lstsq(A, b)[0]
LAMBDAS = np.array([0.41, 0.4, 0.37, 0.44, 0.39])[:, np.newaxis]
COEFS = np.array([-1854.8214151, 3516.89893646, 221.29346712,
                  128.12323805, -2010.49422654])[:, np.newaxis]


class GaussianProcessClassifier(BaseEstimator, ClassifierMixin):
    """Gaussian process classification (GPC) based on Laplace approximation.

    The implementation is based on Algorithm 3.1, 3.2, and 5.1 of
    ``Gaussian Processes for Machine Learning'' (GPML) by Rasmussen and
    Williams.

    Internally, the Laplace approximation is used for approximating the
    non-Gaussian posterior by a Gaussian.

    Currently, the implementation is restricted to
      * using the logistic link function
      * and binary classification

    Parameters
    ----------
    kernel : kernel object
        The kernel specifying the covariance function of the GP. If None is
        passed, the kernel "1.0 * RBF(1.0)" is used as default. Note that
        the kernel's hyperparameters are optimized during fitting.

    jitter : float, optional (default: 0.0)
        Value added to the diagonal of the kernel matrix during fitting.
        Larger values correspond to increased noise level in the observations
        and reduce potential numerical issue during fitting.

    optimizer : string, optional (default: "fmin_l_bfgs_b")
        A string specifying the optimization algorithm used for optimizing the
        kernel's parameters. Default uses 'fmin_l_bfgs_b' algorithm from
        scipy.optimize. If None, the kernel's parameters are kept fixed.
        Available optimizers are::

            'fmin_l_bfgs_b'

    n_restarts_optimizer: int, optional (default: 1)
        The number of restarts of the optimizer for finding the kernel's
        parameters which maximize the log-marginal likelihood. The first run
        of the optimizer is performed from the kernel's initial parameters,
        the remaining ones (if any) from thetas sampled log-uniform randomly
        from the space of allowed theta-values. If greater than 1, all bounds
        must be finite.

    warm_start : bool, optional (default: False)
        If warm-starts are enabled, the solution of the last Newton iteration
        on the Laplace approximation of the posterior mode is used as
        initialization for the next call of _posterior_mode(). This can speed
        up convergence when _posterior_mode is called several times on similar
        problems as in hyperparameter optimization.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Attributes
    ----------
    X_fit_ : array-like, shape = (n_samples, n_features)
        Feature values in training data (also required for prediction)

    y_fit_: array-like, shape = (n_samples,)
        Target values in training data (also required for prediction)

    classes_ : array-like, shape = (n_classes,)
        Unique class labels.

    kernel_: kernel object
        The kernel used for prediction. The structure of the kernel is the
        same as the one passed as parameter but with optimized hyperparameters

    L_: array-like, shape = (n_samples, n_samples)
        Lower-triangular Cholesky decomposition of the kernel in X_fit_

    pi_: array-like, shape = (n_samples,)
        The probabilities of the positive class for the training points X_fit_

    W_sr_: array-like, shape = (n_samples,)
        Square root of W, the Hessian of log-likelihood of the latent function
        values for the observed labels. Since W is diagonal, only the diagonal
        of sqrt(W) is stored.
    """

    def __init__(self, kernel=None, jitter=0.0, optimizer="fmin_l_bfgs_b",
                 n_restarts_optimizer=1, warm_start=False, random_state=None):
        self.kernel = kernel
        self.jitter = jitter
        self.optimizer = optimizer
        self.n_restarts_optimizer = n_restarts_optimizer
        self.warm_start = warm_start
        self.rng = check_random_state(random_state)

    def fit(self, X, y):
        """Fit Gaussian process regression model

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Training data

        y : array-like, shape = (n_samples,)
            Target values, must be binary

        Returns
        -------
        self : returns an instance of self.
        """
        if self.kernel is None:  # Use an RBF kernel as default
            self.kernel_ = 1.0 * RBF(1.0)
        else:
            self.kernel_ = clone(self.kernel)

        X, y = check_X_y(X, y)

        self.X_fit_ = X

        # Encode class labels and check that it is a binary classification
        # problem
        label_encoder = LabelEncoder()
        self.y_fit_ = label_encoder.fit_transform(y)
        self.classes_ = label_encoder.classes_
        if self.classes_.size > 2:
            raise ValueError("GaussianProcessClassifier supports only binary "
                             "classification. y contains classes %s"
                             % self.classes_)
        elif self.classes_.size == 1:
            warnings.warn("Only one class label (%s) occurrs in training set."
                          % self.classes_)
            self.classes_ = np.array([self.classes_[0], self.classes_[0]])

        if self.optimizer in ["fmin_l_bfgs_b"]:
            # Choose hyperparameters based on maximizing the log-marginal
            # likelihood (potentially starting from several initial values)
            def obj_func(theta):
                lml, grad = self.log_marginal_likelihood(theta,
                                                         eval_gradient=True)
                return -lml, -grad

            # First optimize starting from theta specified in kernel
            optima = [(self._constrained_optimization(obj_func,
                                                      self.kernel_.theta,
                                                      self.kernel_.bounds))]

            # Additional runs are performed from log-uniform chosen initial
            # theta
            if self.n_restarts_optimizer > 1:
                if not np.isfinite(self.kernel_.bounds).all():
                    raise ValueError(
                        "Multiple optimizer restarts (n_restarts_optimizer>1) "
                        "requires that all bounds are finite.")
                bounds = self.kernel_.bounds
                for iteration in range(1, self.n_restarts_optimizer):
                    theta_initial = np.exp(self.rng.uniform(bounds[:, 0],
                                                            bounds[:, 1]))
                    optima.append(
                        self._constrained_optimization(obj_func, theta_initial,
                                                       bounds))
            # Select result from run with minimal (negative) log-marginal
            # likelihood
            self.kernel_.theta = optima[np.argmin(map(itemgetter(1), optima))][0]
        elif self.optimizer is None:
            pass
        else:
            raise ValueError("Unknown optimizer %s." % self.optimizer)

        # Precompute quantities required for predictions which are independent
        # of actual query points
        K = self.kernel_(self.X_fit_)
        K[np.diag_indices_from(K)] += self.jitter

        _, (self.pi_, self.W_sr_, self.L_, _, _) = \
            self._posterior_mode(K, return_temporaries=True)

        return self

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)

        Returns
        -------
        C : array, shape = (n_samples,)
            Predicted target values for X, values are from classes_
        """
        check_is_fitted(self, ["X_fit_", "y_fit_", "pi_", "W_sr_", "L_"])
        X = check_array(X)

        # As discussed on Section 3.4.2 of GPML, for making hard binary
        # decisions, it is enough to compute the MAP of the posterior and
        # pass it through the link function
        K_star = self.kernel_(self.X_fit_, X)  # K_star =k(x_star)
        f_star = K_star.T.dot(self.y_fit_ - self.pi_)  # Line 4 (Algorithm 3.2)

        return np.where(f_star > 0, self.classes_[1], self.classes_[0])

    def predict_proba(self, X):
        """Return probability estimates for the test vector X.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)

        Returns
        -------
        C : array-like, shape = (n_samples, n_classes)
            Returns the probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.
        """
        check_is_fitted(self, ["X_fit_", "y_fit_", "pi_", "W_sr_", "L_"])
        X = check_array(X)

        # Based on Algorithm 3.2 of GPML
        K_star = self.kernel_(self.X_fit_, X)  # K_star =k(x_star)
        f_star = K_star.T.dot(self.y_fit_ - self.pi_)  # Line 4
        v = solve(self.L_, self.W_sr_[:, np.newaxis] * K_star)  # Line 5
        # Line 6 (compute np.diag(v.T.dot(v)) via einsum)
        var_f_star = self.kernel_.diag(X) - np.einsum("ij,ij->j", v, v)

        # Line 7:
        # Approximate \int log(z) * N(z | f_star, var_f_star)
        # Approximation is due to Williams & Barber, "Bayesian Classification
        # with Gaussian Processes", Appendix A: Approximate the logistic
        # sigmoid by a linear combination of 5 error functions.
        # For information on how this integral can be computed see
        # blitiri.blogspot.de/2012/11/gaussian-integral-of-error-function.html
        alpha = 1 / (2 * var_f_star)
        gamma = LAMBDAS * f_star
        integrals = np.sqrt(np.pi / alpha) \
            * erf(gamma * np.sqrt(alpha / (alpha + LAMBDAS**2))) \
            / (2 * np.sqrt(var_f_star * 2 * np.pi))
        pi_star = (COEFS * integrals).sum(axis=0) + .5 * COEFS.sum()

        return np.vstack((1 - pi_star, pi_star)).T

    def log_marginal_likelihood(self, theta, eval_gradient=False):
        """Returns log-marginal likelihood of theta for training data.

        Parameters
        ----------
        theta : array-like, shape = (n_kernel_params,)
            Kernel hyperparameters for which the log-marginal likelihood is
            evaluated

        eval_gradient : bool, default: False
            If True, the gradient of the log-marginal likelihood with respect
            to the kernel hyperparameters at position theta is returned
            additionally.

        Returns
        -------
        log_likelihood : float
            Log-marginal likelihood of theta for training data.

        log_likelihood_gradient : array, shape = (n_kernel_params,), optional
            Gradient of the log-marginal likelihood with respect to the kernel
            hyperparameters at position theta.
            Only returned when eval_gradient is True.
        """
        kernel = self.kernel_.clone_with_theta(theta)

        if eval_gradient:
            K, K_gradient = kernel(self.X_fit_, eval_gradient=True)
        else:
            K = kernel(self.X_fit_)

        K[np.diag_indices_from(K)] += self.jitter

        # Compute log-marginal-likelihood Z and also store some temporaries
        # which can be reused for computing Z's gradient
        Z, (pi, W_sr, L, b, a) = \
            self._posterior_mode(K, return_temporaries=True)

        if not eval_gradient:
            return Z

        # Compute gradient based on Algorithm 5.1 of GPML
        d_Z = np.empty(theta.shape[0])
        # XXX: Get rid of the np.diag() in the next line
        R = W_sr[:, np.newaxis] * cho_solve((L, True), np.diag(W_sr))  # Line 7
        C = solve(L, W_sr[:, np.newaxis] * K)  # Line 8
        # Line 9:
        s_2 = -0.5*(np.diag(K) - np.diag(C.T.dot(C))) \
            * (pi * (1 - pi) * (1 - 2*pi))  # third derivative
        for j in range(d_Z.shape[0]):
            C = K_gradient[:, :, j]   # Line 11
            s_1 = .5 * a.T.dot(C).dot(a) - .5 * np.trace(R.dot(C))  # Line 12

            b = C.dot(self.y_fit_ - pi)  # Line 13
            s_3 = b - K.dot(R).dot(b)  # Line 14

            d_Z[j] = s_1 + s_2.T.dot(s_3)  # Line 15

        return Z, d_Z

    def _posterior_mode(self, K, return_temporaries=False):
        """Mode-finding for binary Laplace GPC and fixed kernel.

        This approximates the posterior of the latent function values for given
        inputs and target observations with a Gaussian approximation and uses
        Newton's iteration to find the mode of this approximation.
        """
        # Based on Algorithm 3.1 of GPML

        # If warm_start are enabled, we reuse the last solution for the
        # posterior mode as initialization; otherwise, we initialize with 0
        if self.warm_start and hasattr(self, "f_cached") \
           and self.f_cached.shape == self.y_fit_.shape:
            f = self.f_cached
        else:
            f = np.zeros_like(self.y_fit_, dtype=np.float64)

        # Use Newton's iteration method to find mode of Laplace approximation
        log_marginal_likelihood = -np.inf
        while True:
            # Line 4
            pi = 1 / (1 + np.exp(-f))
            W = pi * (1 - pi)
            # Line 5
            W_sr = np.sqrt(W)
            W_sr_K = W_sr[:, np.newaxis] * K
            B = np.eye(W.shape[0]) + W_sr_K * W_sr
            L = cholesky(B, lower=True)
            # Line 6
            b = W * f + (self.y_fit_ - pi)
            # Line 7
            a = b - W_sr * cho_solve((L, True), W_sr_K.dot(b))
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

        self.f_cached = f  # Remember solution for later warm-starts
        if return_temporaries:
            return log_marginal_likelihood, (pi, W_sr, L, b, a)
        else:
            return log_marginal_likelihood

    def _constrained_optimization(self, obj_func, initial_theta, bounds):
        if self.optimizer in ["fmin_l_bfgs_b"]:
            theta_opt, func_min, _ = \
                fmin_l_bfgs_b(obj_func, initial_theta, bounds=bounds)
        return theta_opt, func_min
