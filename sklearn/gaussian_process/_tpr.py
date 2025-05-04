"""T processes regression."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Real

import numpy as np
from scipy.linalg import cho_solve
from scipy.special import gammaln
from scipy.stats import multivariate_t

from ..base import _fit_context
from ..utils import check_random_state
from ..utils._param_validation import Interval
from ._gpr import GPR_CHOLESKY_LOWER, GaussianProcessRegressor


class TProcessRegressor(GaussianProcessRegressor):
    """T Process Regressor (TPR).

    This implementation is primarily based of [SW2014]_. However, due to their unusual
    parametrization of the Student T distribution [TW2018]_ is also referenced.
    Lastly the kernel parameter optimization is largely based off [RW2006]_ as was
    done in GPs.

    In addition to the :class:`GaussianProcessRegressor`, :class:`TProcessRegressor`:
        * allows variance of predictions to be predicted
        * allows variance of predictions to scaled using Bayesian Methods.

    Parameters
    ----------
    kernel : kernel instance, default=None
        The kernel specifying the covariance function of the TP. Unlike a GP, the
        covariance function is used to construct the shape parameter of the
        inverse-gamma prior over the multivariate normal's covariance matrix.

    degrees_of_freedom_init : float, default=2
        The starting degrees of freedom of the multivariate T distribution. As
        observations are added, the degrees of freedom will increase. This
        value orriginates as the degrees of freedom of the inverse gamma prior.

    alpha : float or ndarray of shape (n_samples,), default=1e-10
        Value added to the diagonal of the kernel matrix during fitting.
        This can prevent a potential numerical issue during fitting, by
        ensuring that the calculated values form a positive definite matrix.
        It can also be interpreted as the variance of additional Gaussian
        measurement noise on the training observations. Note however this is
        not independent noise, it is noise added to kernel's prior. Further note,
        that this is different from using a `WhiteKernel`. If an array is passed,
        it must have the same number of entries as the data used for fitting and is
        used as datapoint-dependent noise level. Allowing to specify the
        noise level directly as a parameter is mainly for convenience and
        for consistency with :class:`~sklearn.linear_model.Ridge`.

    optimizer : "fmin_l_bfgs_b", callable or None, default="fmin_l_bfgs_b"
        Can either be one of the internally supported optimizers for optimizing
        the kernel's parameters, specified by a string, or an externally
        defined optimizer passed as a callable. If a callable is passed, it
        must have the signature::

            def optimizer(obj_func, initial_theta, bounds):
                # * 'obj_func': the objective function to be minimized, which
                #   takes the hyperparameters theta as a parameter and an
                #   optional flag eval_gradient, which determines if the
                #   gradient is returned additionally to the function value
                # * 'initial_theta': the initial value for theta, which can be
                #   used by local optimizers
                # * 'bounds': the bounds on the values of theta
                ....
                # Returned are the best found hyperparameters theta and
                # the corresponding value of the target function.
                return theta_opt, func_min

        Per default, the L-BFGS-B algorithm from `scipy.optimize.minimize`
        is used. If None is passed, the kernel's parameters are kept fixed.
        Available internal optimizers are: `{'fmin_l_bfgs_b'}`.

    n_restarts_optimizer : int, default=0
        The number of restarts of the optimizer for finding the kernel's
        parameters which maximize the log-marginal likelihood. The first run
        of the optimizer is performed from the kernel's initial parameters,
        the remaining ones (if any) from thetas sampled log-uniform randomly
        from the space of allowed theta-values. If greater than 0, all bounds
        must be finite. Note that `n_restarts_optimizer == 0` implies that one
        run is performed.

    normalize_y : bool, default=False
        Whether or not to normalize the target values `y` by removing the mean
        and scaling to unit-variance. This is recommended for cases where
        zero-mean, unit-variance priors are used. Note that, in this
        implementation, the normalisation is reversed before the GP predictions
        are reported.

    copy_X_train : bool, default=True
        If True, a persistent copy of the training data is stored in the
        object. Otherwise, just a reference to the training data is stored,
        which might cause predictions to change if the data is modified
        externally.

    n_targets : int, default=None
        The number of dimensions of the target values. Used to decide the number
        of outputs when sampling from the prior distributions (i.e. calling
        :meth:`sample_y` before :meth:`fit`). This parameter is ignored once
        :meth:`fit` has been called.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation used to initialize the centers.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    X_train_ : array-like of shape (n_samples, n_features) or list of object
        Feature vectors or other representations of training data (also
        required for prediction).

    y_train_ : array-like of shape (n_samples,) or (n_samples, n_targets)
        Target values in training data (also required for prediction).

    kernel_ : kernel instance
        The kernel used for prediction. The structure of the kernel is the
        same as the one passed as parameter but with optimized hyperparameters.

    L_ : array-like of shape (n_samples, n_samples)
        Lower-triangular Cholesky decomposition of the kernel in ``X_train_``.

    alpha_ : array-like of shape (n_samples,)
        Dual coefficients of training data points in kernel space.

    log_marginal_likelihood_value_ : float
        The log-marginal-likelihood of ``self.kernel_.theta``.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    GaussianProcessRegressor : Gaussian process regression (GPR).

    References
    ----------
    .. [SW2014] `Amar Shah, Andrew Gordon Wilson, Zoubin Ghahramani,
       "Student-t Processes as Alternatives to Gaussian Processes",
       arxiv > stat > arXiv:1402.4306 2014 <https://arxiv.org/abs/1402.4306>`_

    .. [TW2018] `Brendan D. Tracey, David H. Wolpert,
       "Upgrading from Gaussian Processes to Student's-T Processes"
        arxiv > stat > arXiv:1801.06147 2018 <https://arxiv.org/abs/1801.06147>`_

    .. [RW2006] `Carl E. Rasmussen and Christopher K.I. Williams,
       "Gaussian Processes for Machine Learning",
       MIT Press 2006 <https://www.gaussianprocess.org/gpml/chapters/RW.pdf>`_

    .. [S2024] `Stevens Conrad
       "Towards Analytic Informative Path Planning",
       The University of Sydney Library 2024 <https://hdl.handle.net/2123/33299>`_

    Examples
    --------
    >>> from sklearn.datasets import make_friedman2
    >>> from sklearn.gaussian_process import TProcessRegressor
    >>> from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
    >>> X, y = make_friedman2(n_samples=500, noise=0, random_state=0)
    >>> kernel = DotProduct() + WhiteKernel()
    >>> tpr = TProcessRegressor(kernel=kernel,
    ...         random_state=0).fit(X, y)
    >>> tpr.score(X, y)
    0.8690...
    >>> tpr.predict(X[:2,:], return_std=True)
    (array([754.5..., 526.2...]), array([147.8..., 148.0...]))
    """

    _parameter_constraints: dict = {
        **GaussianProcessRegressor._parameter_constraints,
        "degrees_of_freedom_init": [Interval(Real, 3, None, closed="left"), np.ndarray],
    }

    def __init__(
        self,
        kernel=None,
        degrees_of_freedom_init=3,
        *,
        alpha=1e-10,
        optimizer="fmin_l_bfgs_b",
        n_restarts_optimizer=0,
        normalize_y=False,
        copy_X_train=True,
        n_targets=None,
        random_state=None,
    ):
        super().__init__(
            kernel,
            alpha=alpha,
            optimizer=optimizer,
            n_restarts_optimizer=n_restarts_optimizer,
            normalize_y=normalize_y,
            copy_X_train=copy_X_train,
            n_targets=n_targets,
            random_state=random_state,
        )
        self.degrees_of_freedom_init = degrees_of_freedom_init

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y):
        """Fit T process regression model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Feature vectors or other representations of training data.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values.

        Returns
        -------
        self : object
            TProcessRegressor class instance.
        """
        super().fit(X, y)
        return self

    def predict(
        self,
        X,
        return_std=False,
        return_cov=False,
        return_t_shape=False,
        return_t_shape_matrix=False,
    ):
        """Predict using the T process regression model.

        The joint predictive distribution at query points is a multivariate
        Student-t distribution. In addition to the mean of the predictive
        distribution, optionally also returns its standard deviation
        (`return_std=True`), covariance (`return_cov=True`) or shape parameter
        (`return_t_shape=True`). Note that at most one of the three can be
        requested.

        We can also predict based on an unfitted model by using the TP prior.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Query points where the GP is evaluated.

        return_std : bool, default=False
            If True, the standard-deviation of the predictive distribution at
            the query points is returned along with the mean.

        return_cov : bool, default=False
            If True, the covariance of the joint predictive distribution at
            the query points is returned along with the mean.

        return_t_shape : bool, default=False
            If True, the shape parameter of the predictive t distribution at
            the query points is returned along with the mean.

        return_t_shape_matrix : bool, default=False
            If True, the shape parameter of the joint predictive t distribution
            the query points is returned along with the mean.

        Returns
        -------
        y_mean : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Mean of predictive distribution at query points.

        y_std : ndarray of shape (n_samples,) or (n_samples, n_targets), optional
            Standard deviation of predictive distribution at query points.
            Only returned when `return_std` is True.

        y_cov : ndarray of shape (n_samples, n_samples) or \
                (n_samples, n_samples, n_targets), optional
            Covariance of joint predictive distribution at query points.
            Only returned when `return_cov` is True.

        y_t_shape : ndarray of shape (n_samples, n_samples) or \
                (n_samples, n_samples, n_targets), optional
            Shape of joint predictive t distribution at query points.
            Only returned when `return_cov` is True.
        """
        if [return_std, return_cov, return_t_shape, return_t_shape_matrix].count(
            True
        ) > 1:
            raise RuntimeError(
                "At most one of return_std, return_cov, return_t_shape or "
                + "return_t_shape_matrix can be requested."
            )

        n = getattr(self, "_n", 0)
        degrees_of_freedom = self.degrees_of_freedom_init + n

        ### Spread may be either std or cov ###
        if not any([return_std, return_cov, return_t_shape, return_t_shape_matrix]):
            return super().predict(X, return_std, return_cov)
        elif return_cov or return_t_shape_matrix:  # Return a matrix:
            y_mean, y_spread = super().predict(X, return_cov=True)
        else:  # Return spread metric at points
            y_mean, y_spread = super().predict(X, return_std=True)

        ### Adjust depending on desired posterior ###
        if n > 0 and (return_t_shape or return_t_shape_matrix):
            scaling_factor = (
                self._mdis + self.degrees_of_freedom_init - 2
            ) / degrees_of_freedom
        elif n > 0 and (return_std or return_cov):
            scaling_factor = (self._mdis + self.degrees_of_freedom_init - 2) / (
                degrees_of_freedom - 2
            )
        elif n == 0 and (return_t_shape or return_t_shape_matrix):
            scaling_factor = (degrees_of_freedom - 2) / degrees_of_freedom
        else:
            scaling_factor = 1

        if return_std:
            y_spread = y_spread * np.sqrt(scaling_factor)
        else:
            y_spread = y_spread * scaling_factor

        return y_mean, y_spread

    def sample_y(self, X, n_samples=1, random_state=0):
        """Draw samples from T-process and evaluate at X.

        Parameters
        ----------
        X : array-like of shape (n_samples_X, n_features) or list of object
            Query points where the GP is evaluated.

        n_samples : int, default=1
            Number of samples drawn from the T-process per query point.

        random_state : int, RandomState instance or None, default=0
            Determines random number generation to randomly draw samples.
            Pass an int for reproducible results across multiple function
            calls.
            See :term:`Glossary <random_state>`.

        Returns
        -------
        y_samples : ndarray of shape (n_samples_X, n_samples), or \
            (n_samples_X, n_targets, n_samples)
            Values of n_samples samples drawn from T-process and
            evaluated at query points.
        """
        rng = check_random_state(random_state)
        n = getattr(self, "_n", 0)
        degrees_of_freedom = self.degrees_of_freedom_init + n

        y_mean, y_t_shape_matrix = self.predict(X, return_t_shape_matrix=True)
        if y_mean.ndim == 0:
            y_samples = (
                multivariate_t(
                    np.array([0]), y_t_shape_matrix, degrees_of_freedom, seed=rng
                )
                .rvs(n_samples)
                .T
            )
        elif y_mean.ndim == 1:
            y_samples = (
                multivariate_t(y_mean, y_t_shape_matrix, degrees_of_freedom, seed=rng)
                .rvs(n_samples)
                .T
            )
        else:
            y_samples = [
                multivariate_t(
                    y_mean[:, target],
                    y_t_shape_matrix[..., target],
                    degrees_of_freedom,
                    seed=rng,
                )
                .rvs(n_samples)
                .T[:, np.newaxis]
                for target in range(y_mean.shape[1])
            ]
            y_samples = np.hstack(y_samples)
        return y_samples

    @_fit_context(prefer_skip_nested_validation=True)
    def _preliminary_data_check(self, X, y):
        """Checks that the data the TP is trianing on is valid.
        Then updates reused internal calculated values.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Feature vectors or other representations of training data.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values.
        """
        X, y = super()._preliminary_data_check(X, y)
        self._n = y.shape[0] if isinstance(y, np.ndarray) else 1
        self._degrees_of_freedom = self.degrees_of_freedom_init + self._n
        self._log_likelihood_dims_const = (
            gammaln(self._degrees_of_freedom / 2)
            - gammaln(self.degrees_of_freedom_init / 2)
            - self._n / 2 * np.log(self.degrees_of_freedom_init * np.pi)
        )
        return X, y

    def _log_likelihood_calc(self, y_train, alpha, L, K):
        """Returns the log-likelihood given L and the training points.

        Parameters
        ----------
        y_train : array-like of shape (n_samples,) or (n_samples, n_targets)
                  Target values.

        alpha : K^(-1) * y_train

        L : Lower cholesky decomposition of the kernel matrix K.

        K : Kernel matrix used.

        Returns
        -------
        log_likelihood : float
            Log-marginal likelihood of multivariate T distribution
            using covariance K and training data
        """
        # Log-likelihood function Eq (27) in [TW2018]
        ### Change to shape of kernel Parameter ###
        L = (
            L
            * ((self.degrees_of_freedom_init - 2) / self.degrees_of_freedom_init) ** 0.5
        )

        self._mdis = np.einsum("ik,ik->k", y_train, alpha)
        self._shape_m_dism_dis = (
            self._mdis
            * self.degrees_of_freedom_init
            / (self.degrees_of_freedom_init - 2)
        )
        log_likelihood_dims = self._log_likelihood_dims_const
        log_likelihood_dims -= (
            self._degrees_of_freedom
            / 2
            * np.log(1 + self._shape_m_dism_dis / self.degrees_of_freedom_init)
        )
        log_likelihood_dims -= np.log(np.diag(L)).sum()
        log_likelihood = log_likelihood_dims.sum(axis=-1)
        return log_likelihood

    def _log_likelihood_gradient_calc(self, alpha, L, K, K_gradient):
        """Returns the log-likelihood gradient given the required algebraic terms.

        Parameters
        ----------
        y_train : array-like of shape (n_samples,) or (n_samples, n_targets)
                  Target values.

        alpha : K^(-1) * y_train

        L : Lower cholesky decomposition of the kernel matrix K.

        K : Kernel matrix used.

        Returns
        -------
        log_likelihood_gradient : np.array
            Log-marginal likelihood gradient with respect to theta
        """
        # Derivative of the Log-likelihood function can be found in [TW2018]
        # Optimization is based of [S2024] building off algorithms from
        # [RW2006] used for Gaussian Processes (._gpr.GaussianProcessRegressor)
        ### Change to shape of kernel Parameter ###
        L = (
            L
            * ((self.degrees_of_freedom_init - 2) / self.degrees_of_freedom_init) ** 0.5
        )
        alpha = (
            alpha * self.degrees_of_freedom_init / (self.degrees_of_freedom_init - 2)
        )
        K_gradient = (
            K_gradient
            * (self.degrees_of_freedom_init - 2)
            / self.degrees_of_freedom_init
        )

        inner_term = np.einsum("ik,jk->ijk", alpha, alpha)
        inner_term = (
            self._degrees_of_freedom
            / (self.degrees_of_freedom_init + self._shape_m_dism_dis)
            * inner_term
        )
        # compute K^-1 of shape (n_samples, n_samples)
        K_inv = cho_solve(
            (L, GPR_CHOLESKY_LOWER), np.eye(K.shape[0]), check_finite=False
        )
        # create a new axis to use broadcasting between inner_term and
        # K_inv
        inner_term -= K_inv[..., np.newaxis]
        # Since we are interested about the trace of
        # inner_term @ K_gradient, we don't explicitly compute the
        # matrix-by-matrix operation and instead use an einsum. Therefore
        # it is equivalent to:
        # for param_idx in range(n_kernel_params):
        #     for output_idx in range(n_output):
        #         log_likehood_gradient_dims[param_idx, output_idx] = (
        #             inner_term[..., output_idx] @
        #             K_gradient[..., param_idx]
        #         )
        log_likelihood_gradient_dims = 0.5 * np.einsum(
            "ijl,jik->kl", inner_term, K_gradient
        )
        # the log likehood gradient is the sum-up across the outputs
        log_likelihood_gradient = log_likelihood_gradient_dims.sum(axis=-1)
        return log_likelihood_gradient
