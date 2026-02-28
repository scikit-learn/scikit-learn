"""Gaussian processes regression."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Integral, Real

import numpy as np
import scipy.optimize
from scipy.linalg import cho_solve, cholesky, solve_triangular

from sklearn.base import (
    BaseEstimator,
    MultiOutputMixin,
    RegressorMixin,
    _fit_context,
    clone,
)
from sklearn.gaussian_process.kernels import RBF, Kernel, _is_sequence_like
from sklearn.gaussian_process.kernels import ConstantKernel as C
from sklearn.preprocessing._data import _handle_zeros_in_scale
from sklearn.utils import check_random_state
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.optimize import _check_optimize_result
from sklearn.utils.validation import validate_data

GPR_CHOLESKY_LOWER = True


class GaussianProcessRegressor(MultiOutputMixin, RegressorMixin, BaseEstimator):
    """Gaussian process regression (GPR).

    The implementation is based on Algorithm 2.1 of [RW2006]_.

    In addition to standard scikit-learn estimator API,
    :class:`GaussianProcessRegressor`:

    * allows prediction without prior fitting, based on the Gaussian process (GP) prior,
    * provides an additional method :meth:`sample_y`,
      which evaluates samples drawn from the GP (prior or posterior) at given inputs,
    * exposes a method :meth:`log_marginal_likelihood`,
      which can be used externally for other ways of selecting hyperparameters,
      e.g., via Markov chain Monte Carlo.

    To learn the difference
    between a point-estimate approach vs. a more Bayesian modelling approach,
    refer to the example entitled
    :ref:`sphx_glr_auto_examples_gaussian_process_plot_compare_gpr_krr.py`.

    Read more in the :ref:`User Guide <gaussian_process>`.

    .. versionadded:: 0.18

    Parameters
    ----------
    kernel : kernel instance, default=None
        The kernel specifying the covariance function of the GP.
        If `None` is passed,
        :class:`~sklearn.gaussian_process.kernels.RBF` with fixed parameters is used.
        Note that the kernel hyperparameters are optimized during fitting
        unless the bounds are marked as "fixed".

    alpha : float or ndarray of shape (n_samples,), default=1e-10
        Value added to the diagonal of the kernel matrix during fitting.
        This can prevent a potential numerical issue during fitting,
        by ensuring that the calculated values form a positive definite matrix.
        It can also be interpreted as
        the variance of a Gaussian measurement noise on the training target values.
        Note that this is different
        from using a :class:`~sklearn.gaussian_process.kernels.WhiteKernel`.
        If an array is passed,
        it must have the same number of entries as the number of training samples
        and is used as datapoint-dependent noise level.
        Allowing to specify the noise level directly as a parameter is mainly
        for convenience and for consistency with :class:`~sklearn.linear_model.Ridge`.
        For an example
        illustrating how the `alpha` parameter controls the noise variance in GPR,
        see :ref:`sphx_glr_auto_examples_gaussian_process_plot_gpr_noisy_targets.py`.

    optimizer : "fmin_l_bfgs_b", callable or None, default="fmin_l_bfgs_b"
        Used for optimizing the kernel hyperparameters.
        Can either be one of the internally supported optimizers specified by a string,
        or an externally defined optimizer passed as a callable.
        If a callable is passed,
        it must have the signature:

            def optimizer(obj_func, initial_theta, bounds):
                # * 'obj_func': the objective function to be minimized,
                #   which takes the hyperparameters `theta` as an argument
                #   and an optional argument `eval_gradient`,
                #   which determines if the gradient is returned
                #   additionally to the function value,
                # * 'initial_theta': the initial value for theta,
                #   which can be used by local optimizers,
                # * 'bounds': the bounds on the values of theta.
                ....
                # Returned are the best found hyperparameters theta
                # and the corresponding value of the target function.
                return theta_opt, func_min

        Per default,
        the L-BFGS-B algorithm from `scipy.optimize.minimize` is used.
        If `None` is passed,
        the kernel hyperparameters are kept fixed.
        Available internal optimizers are: `{'fmin_l_bfgs_b'}`.

    n_restarts_optimizer : int, default=0
        The number of restarts of the optimizer
        for finding the kernel hyperparameters maximizing the log-marginal likelihood,
        in addition to its first run from the initial hyperparameters.
        These repetititions starts from new initial hyperparameters
        uniformly sampled over the bounded hyperparameter space.

    normalize_y : bool, default=False
        Whether to normalize the target values `y`
        by removing the mean and scaling to unit-variance.
        This is recommended for cases where zero-mean, unit-variance priors are used.
        Note that,
        in this implementation,
        the normalisation is reversed before the GP predictions are reported.

        .. versionchanged:: 0.23

    copy_X_train : bool, default=True
        If `True`,
        a persistent copy of the training data is stored in the object.
        Otherwise,
        just a reference to the training data is stored,
        which might cause predictions to change if the data is modified externally.

    n_targets : int, default=None
        The number of dimensions of the target values.
        Used to decide the number of outputs when sampling from the prior distributions
        (i.e. calling :meth:`sample_y` before :meth:`fit`).
        This argument is ignored once :meth:`fit` has been called.

        .. versionadded:: 1.3

    random_state : int, RandomState instance or None, default=None
        Determines random number generation used to initialize the centers.
        Pass an integer for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    X_train_ : array-like of shape (n_samples, n_features) or list of object
        Feature vectors or other representations of training data
        (also required for prediction).

    y_train_ : array-like of shape (n_samples,) or (n_samples, n_targets)
        Target values in training data (also required for prediction);
        these target values are normalized when ``normalize_y`` is ``True``.

    kernel_ : kernel instance
        The kernel used for prediction.
        The structure of the kernel is the same as the one passed as argument
        but with optimized hyperparameters.

    L_ : array-like of shape (n_samples, n_samples)
        Lower-triangular Cholesky decomposition of the kernel in ``X_train_``.

    alpha_ : array-like of shape (n_samples,)
        Dual coefficients of training data points in kernel space.

    log_marginal_likelihood_value_ : float
        The log-marginal-likelihood of ``self.kernel_.theta``.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    GaussianProcessClassifier : Gaussian process classification (GPC)
        based on Laplace approximation.

    References
    ----------
    .. [RW2006] `Carl E. Rasmussen and Christopher K.I. Williams,
       "Gaussian Processes for Machine Learning",
       MIT Press 2006 <https://www.gaussianprocess.org/gpml/chapters/RW.pdf>`_

    Examples
    --------
    >>> from sklearn.datasets import make_friedman2
    >>> from sklearn.gaussian_process import GaussianProcessRegressor
    >>> from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
    >>> X, y = make_friedman2(n_samples=500, noise=0, random_state=0)
    >>> kernel = DotProduct() + WhiteKernel()
    >>> gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(X, y)
    >>> gpr.score(X, y)
    0.3680...
    >>> gpr.predict(X[:2,:], return_std=True)
    (array([653.0, 592.1]), array([316.6, 316.6]))
    """

    _parameter_constraints: dict = {
        "kernel": [None, Kernel],
        "alpha": [Interval(Real, 0, None, closed="left"), np.ndarray],
        "optimizer": [StrOptions({"fmin_l_bfgs_b"}), callable, None],
        "n_restarts_optimizer": [Interval(Integral, 0, None, closed="left")],
        "normalize_y": ["boolean"],
        "copy_X_train": ["boolean"],
        "n_targets": [Interval(Integral, 1, None, closed="left"), None],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        kernel=None,
        *,
        alpha=1e-10,
        optimizer="fmin_l_bfgs_b",
        n_restarts_optimizer=0,
        normalize_y=False,
        copy_X_train=True,
        n_targets=None,
        random_state=None,
    ):
        self.kernel = kernel
        self.alpha = alpha
        self.optimizer = optimizer
        self.n_restarts_optimizer = n_restarts_optimizer
        self.normalize_y = normalize_y
        self.copy_X_train = copy_X_train
        self.n_targets = n_targets
        self.random_state = random_state

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y):
        """Fit GPR model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Feature vectors or other representations of training data.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values.

        Returns
        -------
        self : object
            GaussianProcessRegressor class instance.
        """
        self.kernel_ = self.__get_kernel(True)

        self._rng = check_random_state(self.random_state)

        dtype, ensure_2d = self.__get_dtype_ensure_2d()
        X, y = validate_data(
            self,
            X,
            y,
            multi_output=True,
            y_numeric=True,
            ensure_2d=ensure_2d,
            dtype=dtype,
        )

        n_targets_seen = y.shape[1] if y.ndim > 1 else 1
        if self.n_targets is not None and n_targets_seen != self.n_targets:
            raise ValueError(
                "The number of targets seen in `y` is different from the parameter "
                f"`n_targets`. Got {n_targets_seen} != {self.n_targets}."
            )

        # Normalize target value
        if self.normalize_y:
            self._y_train_mean = np.mean(y, axis=0)
            self._y_train_std = _handle_zeros_in_scale(np.std(y, axis=0), copy=False)

            # Remove mean and make unit variance
            y = (y - self._y_train_mean) / self._y_train_std

        else:
            shape_y_stats = (y.shape[1],) if y.ndim == 2 else 1
            self._y_train_mean = np.zeros(shape=shape_y_stats)
            self._y_train_std = np.ones(shape=shape_y_stats)

        alpha_length = len(self.alpha) if _is_sequence_like(self.alpha) else 0
        n_samples = len(y)
        if alpha_length != n_samples:
            if alpha_length == 1:
                self.alpha = self.alpha[0]
            elif alpha_length != 0:
                raise ValueError(
                    "alpha must be a scalar or an array with same number of "
                    f"entries as y. ({alpha_length} != {n_samples})"
                )

        self.X_train_ = np.copy(X) if self.copy_X_train else X
        self.y_train_ = np.copy(y) if self.copy_X_train and not self.normalize_y else y

        if self.optimizer is not None and self.kernel_.n_dims > 0:
            # Find hyperparameters maximizing the log-marginal likelihood (LML):
            bounds = self.kernel_.bounds
            optima = [self.__maximize_lml(self.kernel_.theta, bounds)]
            if self.n_restarts_optimizer > 0:
                # Repeat LML maximization from log-uniform chosen initial theta:
                if not np.isfinite(self.kernel_.bounds).all():
                    raise ValueError(
                        "Multiple optimizer restarts (n_restarts_optimizer>0) "
                        "requires that all bounds are finite."
                    )
                for _ in range(self.n_restarts_optimizer):
                    initial_theta = self._rng.uniform(bounds[:, 0], bounds[:, 1])
                    optima.append(self.__maximize_lml(initial_theta, bounds))

            # Select hyperparameters maximizing the LML across repetitions:
            theta, neg_lml = optima[min(enumerate(optima), key=lambda x: x[1][1])[0]]
            self.kernel_.theta = theta
            self.kernel_._check_bounds_params()
            self.log_marginal_likelihood_value_ = -neg_lml
        else:
            self.log_marginal_likelihood_value_ = self.log_marginal_likelihood(
                self.kernel_.theta, clone_kernel=False
            )

        # Precompute quantities required for predictions
        # which are independent of actual query points
        # Alg. 2.1, page 19, line 2 -> L = cholesky(K + sigma^2 I)
        K = self.kernel_(self.X_train_)
        K[np.diag_indices_from(K)] += self.alpha
        try:
            self.L_ = cholesky(K, lower=GPR_CHOLESKY_LOWER, check_finite=False)
        except np.linalg.LinAlgError as exc:
            exc.args = (
                (
                    f"The kernel, {self.kernel_}, "
                    "is not returning a positive definite matrix. "
                    "Try gradually increasing the 'alpha' parameter "
                    "of your GaussianProcessRegressor estimator."
                ),
            ) + exc.args
            raise
        # Alg 2.1, page 19, line 3 -> alpha = L^T \ (L \ y)
        self.alpha_ = cho_solve(
            (self.L_, GPR_CHOLESKY_LOWER),
            self.y_train_,
            check_finite=False,
        )
        return self

    def predict(self, X, return_std=False, return_cov=False):
        """Predict using the GPR model.

        We can also predict based on an unfitted model by using the GP prior.
        In addition to the mean of the predictive distribution,
        optionally also returns its standard deviation (`return_std=True`)
        or covariance (`return_cov=True`).
        Note that at most one of the two can be requested.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Query points where the GP is evaluated.

        return_std : bool, default=False
            If `True`,
            the standard-deviation at the query points is returned along with the mean.

        return_cov : bool, default=False
            If `True`,
            the covariance between the query points is returned along with the mean.

        Returns
        -------
        y_mean : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Mean at query points.

        y_std : ndarray of shape (n_samples,) or (n_samples, n_targets), optional
            Standard deviation at query points.
            Only returned when `return_std` is True.

        y_cov : ndarray of shape (n_samples, n_samples) or \
                (n_samples, n_samples, n_targets), optional
            Covariance between the query points.
            Only returned when `return_cov` is True.
        """
        if return_std and return_cov:
            raise RuntimeError(
                "At most one of return_std or return_cov can be requested."
            )

        dtype, ensure_2d = self.__get_dtype_ensure_2d()
        X = validate_data(self, X, ensure_2d=ensure_2d, dtype=dtype, reset=False)

        if not hasattr(self, "X_train_"):
            # GPR is unfitted; predict based on GP prior
            kernel = self.__get_kernel(False)
            n_targets = self.n_targets if self.n_targets is not None else 1
            y_mean = np.zeros(shape=(X.shape[0], n_targets)).squeeze()
            if return_cov:
                # Return prior mean and prior covariance
                y_cov = kernel(X)
                y_cov = self.__repeat_by_target(y_cov, n_targets)
                return y_mean, y_cov
            elif return_std:
                # Return prior mean and prior standard deviation
                y_var = kernel.diag(X)
                y_std = self.__repeat_by_target(np.sqrt(y_var), n_targets)
                return y_mean, y_std
            else:
                # Return prior mean
                return y_mean
        else:
            # GPR is fitted; predict based on GP posterior

            # Alg 2.1, page 19, line 4 -> f*_bar = K(X_test, X_train) . alpha
            K_test_train = self.kernel_(X, self.X_train_)
            y_mean = K_test_train @ self.alpha_

            # Undo normalisation
            y_mean = self._y_train_std * y_mean + self._y_train_mean

            y_mean = self.__remove_target_axis_if_scalar(y_mean, 1)
            if not return_cov and not return_std:
                # Return posterior mean
                return y_mean

            # Alg 2.1, page 19, line 5 -> V = L \ K(X_test, X_train)^T
            V = solve_triangular(
                self.L_, K_test_train.T, lower=GPR_CHOLESKY_LOWER, check_finite=False
            )

            if return_cov:
                # Return posterior mean and posterior covariance

                # Alg 2.1, page 19, line 6 -> K(X_test, X_test) - V^T . V
                y_cov = self.kernel_(X) - V.T @ V

                # Undo normalisation
                y_cov = np.outer(y_cov, self._y_train_std**2).reshape(*y_cov.shape, -1)

                y_cov = self.__remove_target_axis_if_scalar(y_cov, 2)
                return y_mean, y_cov
            else:
                # Return posterior mean and posterior standard deviation

                # Posterior variance = diag(A), with A = K(X_test, X_test) - V^T . V
                # Efficiency: compute only the diagonal elements of A.
                y_var = self.kernel_.diag(X).copy() - np.einsum("ij,ji->i", V.T, V)
                y_var = self.__replace_negative_variances_by_zero(y_var)

                # Undo normalisation
                y_var = np.outer(y_var, self._y_train_std**2).reshape(*y_var.shape, -1)

                y_var = self.__remove_target_axis_if_scalar(y_var, 1)
                return y_mean, np.sqrt(y_var)

    @staticmethod
    def __replace_negative_variances_by_zero(y_var):
        y_var_negative = y_var < 0
        if np.any(y_var_negative):
            # There are negative variances due to numerical issues.
            warnings.warn(
                "Predicted variances smaller than 0. Setting those variances to 0."
            )
            y_var[y_var_negative] = 0.0

        return y_var

    def __get_dtype_ensure_2d(self):
        if self.kernel is None or self.kernel.requires_vector_input:
            return "numeric", True

        return None, False

    def sample_y(self, X, n_samples=1, random_state=0):
        """Draw GP samples and evaluate them at query points.

        Parameters
        ----------
        X : array-like of shape (n_samples_X, n_features) or list of object
            Query points.

        n_samples : int, default=1
            Number of GP samples.

        random_state : int, RandomState instance or None, default=0
            Determines random number generation to randomly draw GP samples.
            Pass an integer for reproducible results across multiple function calls.
            See :term:`Glossary <random_state>`.

        Returns
        -------
        y_samples : ndarray of shape (n_samples_X, n_samples), or \
            (n_samples_X, n_targets, n_samples)
            GP samples evaluated at query points.
        """
        sample = check_random_state(random_state).multivariate_normal
        y_mean, y_cov = self.predict(X, return_cov=True)
        if y_mean.ndim == 1:
            return sample(y_mean, y_cov, n_samples).T

        return np.hstack(
            [
                sample(y_mean[:, i], y_cov[..., i], n_samples).T[:, np.newaxis]
                for i in range(y_mean.shape[1])
            ]
        )

    def log_marginal_likelihood(
        self, theta=None, eval_gradient=False, clone_kernel=True
    ):
        """Return log-marginal likelihood of kernel hyperparameters for training data.

        Parameters
        ----------
        theta : array-like of shape (n_kernel_params,) default=None
            Kernel hyperparameters value
            at which the log-marginal likelihood is evaluated.
            If ``None``,
            the value ``self.kernel_.theta`` is considered.

        eval_gradient : bool, default=False
            If ``True``,
            the gradient of the log-marginal likelihood
            with respect to the kernel hyperparameters is evaluated at ``theta``,
            which must not be ``None``.

        clone_kernel : bool, default=True
            If ``True``,
            the kernel attribute is copied.
            If ``False``,
            the kernel attribute is modified,
            but may result in a performance improvement.

        Returns
        -------
        log_likelihood : float
            Log-marginal likelihood of theta for training data.

        log_likelihood_gradient : ndarray of shape (n_kernel_params,), optional
            Gradient of the log-marginal likelihood
            with respect to the kernel hyperparameters at ``theta``.
            Only returned when ``eval_gradient`` is ``True``.
        """
        if theta is None:
            if eval_gradient:
                raise ValueError("Gradient can only be evaluated for theta!=None")
            return self.log_marginal_likelihood_value_

        if clone_kernel:
            kernel = self.kernel_.clone_with_theta(theta)
        else:
            kernel = self.kernel_
            kernel.theta = theta

        if eval_gradient:
            K, K_gradient = kernel(self.X_train_, eval_gradient=True)
        else:
            K = kernel(self.X_train_)

        # Alg. 2.1, page 19, line 2 -> L = cholesky(K + sigma^2 I)
        K[np.diag_indices_from(K)] += self.alpha
        try:
            L = cholesky(K, lower=GPR_CHOLESKY_LOWER, check_finite=False)
        except np.linalg.LinAlgError:
            return (-np.inf, np.zeros_like(theta)) if eval_gradient else -np.inf

        # Target values as a (n_samples, n_targets)-shaped array
        y_train = self.y_train_
        if y_train.ndim == 1:
            y_train = y_train[:, np.newaxis]

        # Alg 2.1, page 19, line 3 -> alpha = L^T \ (L \ y)
        alpha = cho_solve((L, GPR_CHOLESKY_LOWER), y_train, check_finite=False)

        # Alg 2.1, page 19, line 7 (log-likelihood for single target)
        # -> -0.5 . y^T . alpha - sum(log(diag(L))) - n_samples / 2 log(2*pi)
        # We apply this equation to each target:
        n_samples = len(y_train)
        log_likelihoods = (
            -0.5 * np.einsum("ik,ik->k", y_train, alpha)
            - np.log(np.diag(L)).sum()
            - n_samples / 2 * np.log(2 * np.pi)
        )
        # The final log-likelihood is the sum of these log-likelihoods:
        log_likelihood = log_likelihoods.sum(axis=-1)
        if not eval_gradient:
            return log_likelihood

        # Eq. 5.9, p. 114, and footnote 5 in p. 114
        # (log-likelihood gradient for single target)
        # -> 0.5 * trace((alpha . alpha^T - K^-1) . K_gradient)
        # We apply this equation to each target.
        # First,
        # we compute alpha . alpha^T - K^-1 shaped as (n_samples, n_samples, n_targets):
        K_inv = cho_solve(
            (L, GPR_CHOLESKY_LOWER), np.eye(n_samples), check_finite=False
        )
        inner_term = np.einsum("ik,jk->ijk", alpha, alpha) - K_inv[..., np.newaxis]
        # Then, the log-likelihoods shaped as (n_features, n_targets):
        log_likelihood_gradients = 0.5 * np.einsum(
            "ijl,jik->kl", inner_term, K_gradient
        )
        # The final log-likelihood gradient is the sum of these gradients:
        log_likelihood_gradient = log_likelihood_gradients.sum(axis=-1)
        return log_likelihood, log_likelihood_gradient

    def __maximize_lml(self, initial_theta, bounds):
        if self.optimizer == "fmin_l_bfgs_b":
            result = scipy.optimize.minimize(
                self.__evaluate_neg_lml,
                initial_theta,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
            )
            _check_optimize_result("lbfgs", result)
            return result.x, result.fun
        elif callable(self.optimizer):
            return self.optimizer(self.__evaluate_neg_lml, initial_theta, bounds=bounds)
        else:
            raise ValueError(f"Unknown optimizer {self.optimizer}.")

    def __evaluate_neg_lml(self, theta, eval_gradient=True):
        result = self.log_marginal_likelihood(
            theta=theta, eval_gradient=eval_gradient, clone_kernel=False
        )
        return (-result[0], -result[1]) if isinstance(result, tuple) else -result

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.requires_fit = False
        return tags

    def __get_kernel(self, clone_):
        if self.kernel is None:
            return C(constant_value_bounds="fixed") * RBF(length_scale_bounds="fixed")

        return clone(self.kernel) if clone_ else self.kernel

    @staticmethod
    def __remove_target_axis_if_scalar(array_, axis):
        if array_.ndim > axis and array_.shape[axis] == 1:
            return np.squeeze(array_, axis=axis)

        return array_

    @staticmethod
    def __repeat_by_target(array_, n_targets):
        if n_targets > 1:
            return np.expand_dims(array_, -1).repeat(n_targets, axis=-1)

        return array_
