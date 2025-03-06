"""
Various bayesian regression
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from math import log
from numbers import Integral, Real

import numpy as np
from scipy import linalg
from scipy.linalg import pinvh

from ..base import RegressorMixin, _fit_context
from ..utils import _safe_indexing
from ..utils._param_validation import Interval
from ..utils.extmath import fast_logdet
from ..utils.validation import _check_sample_weight, validate_data
from ._base import LinearModel, _preprocess_data, _rescale_data
from ..preprocessing import StandardScaler
from ..utils import check_array

###############################################################################
# BayesianRidge regression


class BayesianRidge(RegressorMixin, LinearModel):
    """Bayesian ridge regression.

    Fit a Bayesian ridge model. See the Notes section for details on this
    implementation and the optimization of the regularization parameters
    lambda (precision of the weights) and alpha (precision of the noise).

    Read more in the :ref:`User Guide <bayesian_regression>`.
    For an intuitive visualization of how the sinusoid is approximated by
    a polynomial using different pairs of initial values, see
    :ref:`sphx_glr_auto_examples_linear_model_plot_bayesian_ridge_curvefit.py`.

    Parameters
    ----------
    max_iter : int, default=300
        Maximum number of iterations over the complete dataset before
        stopping independently of any early stopping criterion.

        .. versionchanged:: 1.3

    tol : float, default=1e-3
        Stop the algorithm if w has converged.

    alpha_1 : float, default=1e-6
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the alpha parameter.

    alpha_2 : float, default=1e-6
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the alpha parameter.

    lambda_1 : float, default=1e-6
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the lambda parameter.

    lambda_2 : float, default=1e-6
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the lambda parameter.

    alpha_init : float, default=None
        Initial value for alpha (precision of the noise).
        If not set, alpha_init is 1/Var(y).

        .. versionadded:: 0.22

    lambda_init : float, default=None
        Initial value for lambda (precision of the weights).
        If not set, lambda_init is 1.

        .. versionadded:: 0.22

    compute_score : bool, default=False
        If True, compute the log marginal likelihood at each iteration of the
        optimization.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model.
        The intercept is not treated as a probabilistic parameter
        and thus has no associated variance. If set
        to False, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    verbose : bool, default=False
        Verbose mode when fitting the model.

    Attributes
    ----------
    coef_ : array-like of shape (n_features,)
        Coefficients of the regression model (mean of distribution)

    intercept_ : float
        Independent term in decision function. Set to 0.0 if
        `fit_intercept = False`.

    alpha_ : float
       Estimated precision of the noise.

    lambda_ : float
       Estimated precision of the weights.

    sigma_ : array-like of shape (n_features, n_features)
        Estimated variance-covariance matrix of the weights

    scores_ : array-like of shape (n_iter_+1,)
        If computed_score is True, value of the log marginal likelihood (to be
        maximized) at each iteration of the optimization. The array starts
        with the value of the log marginal likelihood obtained for the initial
        values of alpha and lambda and ends with the value obtained for the
        estimated alpha and lambda.

    n_iter_ : int
        The actual number of iterations to reach the stopping criterion.

    X_offset_ : ndarray of shape (n_features,)
        If `fit_intercept=True`, offset subtracted for centering data to a
        zero mean. Set to np.zeros(n_features) otherwise.

    X_scale_ : ndarray of shape (n_features,)
        Set to np.ones(n_features).

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    ARDRegression : Bayesian ARD regression.

    Notes
    -----
    There exist several strategies to perform Bayesian ridge regression. This
    implementation is based on the algorithm described in Appendix A of
    (Tipping, 2001) where updates of the regularization parameters are done as
    suggested in (MacKay, 1992). Note that according to A New
    View of Automatic Relevance Determination (Wipf and Nagarajan, 2008) these
    update rules do not guarantee that the marginal likelihood is increasing
    between two consecutive iterations of the optimization.

    References
    ----------
    D. J. C. MacKay, Bayesian Interpolation, Computation and Neural Systems,
    Vol. 4, No. 3, 1992.

    M. E. Tipping, Sparse Bayesian Learning and the Relevance Vector Machine,
    Journal of Machine Learning Research, Vol. 1, 2001.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.BayesianRidge()
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    BayesianRidge()
    >>> clf.predict([[1, 1]])
    array([1.])
    """

    _parameter_constraints: dict = {
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "tol": [Interval(Real, 0, None, closed="neither")],
        "alpha_1": [Interval(Real, 0, None, closed="left")],
        "alpha_2": [Interval(Real, 0, None, closed="left")],
        "lambda_1": [Interval(Real, 0, None, closed="left")],
        "lambda_2": [Interval(Real, 0, None, closed="left")],
        "alpha_init": [None, Interval(Real, 0, None, closed="left")],
        "lambda_init": [None, Interval(Real, 0, None, closed="left")],
        "compute_score": ["boolean"],
        "fit_intercept": ["boolean"],
        "copy_X": ["boolean"],
        "verbose": ["verbose"],
    }

    def __init__(
        self,
        *,
        max_iter=300,
        tol=1.0e-3,
        alpha_1=1.0e-6,
        alpha_2=1.0e-6,
        lambda_1=1.0e-6,
        lambda_2=1.0e-6,
        alpha_init=None,
        lambda_init=None,
        compute_score=False,
        fit_intercept=True,
        copy_X=True,
        verbose=False,
    ):
        self.max_iter = max_iter
        self.tol = tol
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.alpha_init = alpha_init
        self.lambda_init = lambda_init
        self.compute_score = compute_score
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.verbose = verbose

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """Fit the model.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training data.
        y : ndarray of shape (n_samples,)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : ndarray of shape (n_samples,), default=None
            Individual weights for each sample.

            .. versionadded:: 0.20
               parameter *sample_weight* support to BayesianRidge.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        X, y = validate_data(
            self,
            X,
            y,
            dtype=[np.float64, np.float32],
            force_writeable=True,
            y_numeric=True,
        )
        dtype = X.dtype
        n_samples, n_features = X.shape

        sw_sum = n_samples
        y_var = y.var()
        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X, dtype=dtype)
            sw_sum = sample_weight.sum()
            y_mean = np.average(y, weights=sample_weight)
            y_var = np.average((y - y_mean) ** 2, weights=sample_weight)

        X, y, X_offset_, y_offset_, X_scale_ = _preprocess_data(
            X,
            y,
            fit_intercept=self.fit_intercept,
            copy=self.copy_X,
            sample_weight=sample_weight,
        )

        if sample_weight is not None:
            # Sample weight can be implemented via a simple rescaling.
            X, y, _ = _rescale_data(X, y, sample_weight)

        self.X_offset_ = X_offset_
        self.X_scale_ = X_scale_

        # Initialization of the values of the parameters
        eps = np.finfo(np.float64).eps
        # Add `eps` in the denominator to omit division by zero
        alpha_ = self.alpha_init
        lambda_ = self.lambda_init
        if alpha_ is None:
            alpha_ = 1.0 / (y_var + eps)
        if lambda_ is None:
            lambda_ = 1.0

        # Avoid unintended type promotion to float64 with numpy 2
        alpha_ = np.asarray(alpha_, dtype=dtype)
        lambda_ = np.asarray(lambda_, dtype=dtype)

        verbose = self.verbose
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        alpha_1 = self.alpha_1
        alpha_2 = self.alpha_2

        self.scores_ = list()
        coef_old_ = None

        XT_y = np.dot(X.T, y)
        U, S, Vh = linalg.svd(X, full_matrices=False)
        eigen_vals_ = S**2

        # Convergence loop of the bayesian ridge regression
        for iter_ in range(self.max_iter):
            # update posterior mean coef_ based on alpha_ and lambda_ and
            # compute corresponding sse (sum of squared errors)
            coef_, sse_ = self._update_coef_(
                X, y, n_samples, n_features, XT_y, U, Vh, eigen_vals_, alpha_, lambda_
            )
            if self.compute_score:
                # compute the log marginal likelihood
                s = self._log_marginal_likelihood(
                    n_samples,
                    n_features,
                    sw_sum,
                    eigen_vals_,
                    alpha_,
                    lambda_,
                    coef_,
                    sse_,
                )
                self.scores_.append(s)

            # Update alpha and lambda according to (MacKay, 1992)
            gamma_ = np.sum((alpha_ * eigen_vals_) / (lambda_ + alpha_ * eigen_vals_))
            lambda_ = (gamma_ + 2 * lambda_1) / (np.sum(coef_**2) + 2 * lambda_2)
            alpha_ = (sw_sum - gamma_ + 2 * alpha_1) / (sse_ + 2 * alpha_2)

            # Check for convergence
            if iter_ != 0 and np.sum(np.abs(coef_old_ - coef_)) < self.tol:
                if verbose:
                    print("Convergence after ", str(iter_), " iterations")
                break
            coef_old_ = np.copy(coef_)

        self.n_iter_ = iter_ + 1

        # return regularization parameters and corresponding posterior mean,
        # log marginal likelihood and posterior covariance
        self.alpha_ = alpha_
        self.lambda_ = lambda_
        self.coef_, sse_ = self._update_coef_(
            X, y, n_samples, n_features, XT_y, U, Vh, eigen_vals_, alpha_, lambda_
        )
        if self.compute_score:
            # compute the log marginal likelihood
            s = self._log_marginal_likelihood(
                n_samples,
                n_features,
                sw_sum,
                eigen_vals_,
                alpha_,
                lambda_,
                coef_,
                sse_,
            )
            self.scores_.append(s)
            self.scores_ = np.array(self.scores_)

        # posterior covariance is given by 1/alpha_ * scaled_sigma_
        scaled_sigma_ = np.dot(
            Vh.T, Vh / (eigen_vals_ + lambda_ / alpha_)[:, np.newaxis]
        )
        self.sigma_ = (1.0 / alpha_) * scaled_sigma_

        self._set_intercept(X_offset_, y_offset_, X_scale_)

        return self

    def predict(self, X, return_std=False):
        """Predict using the linear model.

        In addition to the mean of the predictive distribution, also its
        standard deviation can be returned.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples.

        return_std : bool, default=False
            Whether to return the standard deviation of posterior prediction.

        Returns
        -------
        y_mean : array-like of shape (n_samples,)
            Mean of predictive distribution of query points.

        y_std : array-like of shape (n_samples,)
            Standard deviation of predictive distribution of query points.
        """
        y_mean = self._decision_function(X)
        if not return_std:
            return y_mean
        else:
            sigmas_squared_data = (np.dot(X, self.sigma_) * X).sum(axis=1)
            y_std = np.sqrt(sigmas_squared_data + (1.0 / self.alpha_))
            return y_mean, y_std

    def _update_coef_(
        self, X, y, n_samples, n_features, XT_y, U, Vh, eigen_vals_, alpha_, lambda_
    ):
        """Update posterior mean and compute corresponding sse (sum of squared errors).

        Posterior mean is given by coef_ = scaled_sigma_ * X.T * y where
        scaled_sigma_ = (lambda_/alpha_ * np.eye(n_features)
                         + np.dot(X.T, X))^-1
        """

        if n_samples > n_features:
            coef_ = np.linalg.multi_dot(
                [Vh.T, Vh / (eigen_vals_ + lambda_ / alpha_)[:, np.newaxis], XT_y]
            )
        else:
            coef_ = np.linalg.multi_dot(
                [X.T, U / (eigen_vals_ + lambda_ / alpha_)[None, :], U.T, y]
            )

        # Note: we do not need to explicitly use the weights in this sum because
        # y and X were preprocessed by _rescale_data to handle the weights.
        sse_ = np.sum((y - np.dot(X, coef_)) ** 2)

        return coef_, sse_

    def _log_marginal_likelihood(
        self, n_samples, n_features, sw_sum, eigen_vals, alpha_, lambda_, coef, sse
    ):
        """Log marginal likelihood."""
        alpha_1 = self.alpha_1
        alpha_2 = self.alpha_2
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2

        # compute the log of the determinant of the posterior covariance.
        # posterior covariance is given by
        # sigma = (lambda_ * np.eye(n_features) + alpha_ * np.dot(X.T, X))^-1
        if n_samples > n_features:
            logdet_sigma = -np.sum(np.log(lambda_ + alpha_ * eigen_vals))
        else:
            logdet_sigma = np.full(n_features, lambda_, dtype=np.array(lambda_).dtype)
            logdet_sigma[:n_samples] += alpha_ * eigen_vals
            logdet_sigma = -np.sum(np.log(logdet_sigma))

        score = lambda_1 * log(lambda_) - lambda_2 * lambda_
        score += alpha_1 * log(alpha_) - alpha_2 * alpha_
        score += 0.5 * (
            n_features * log(lambda_)
            + sw_sum * log(alpha_)
            - alpha_ * sse
            - lambda_ * np.sum(coef**2)
            + logdet_sigma
            - sw_sum * log(2 * np.pi)
        )

        return score


###############################################################################
# ARD (Automatic Relevance Determination) regression


class ARDRegression(RegressorMixin, LinearModel):
    """Bayesian ARD regression.

    Fit the weights of a regression model, using an ARD prior. The weights of
    the regression model are assumed to be in Gaussian distributions.
    Also estimate the parameters lambda (precisions of the distributions of the
    weights) and alpha (precision of the distribution of the noise).
    The estimation is done by an iterative procedures (Evidence Maximization)

    Read more in the :ref:`User Guide <bayesian_regression>`.

    Parameters
    ----------
    max_iter : int, default=300
        Maximum number of iterations.

        .. versionchanged:: 1.3

    tol : float, default=1e-3
        Stop the algorithm if w has converged.

    alpha_1 : float, default=1e-6
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the alpha parameter.

    alpha_2 : float, default=1e-6
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the alpha parameter.

    lambda_1 : float, default=1e-6
        Hyper-parameter : shape parameter for the Gamma distribution prior
        over the lambda parameter.

    lambda_2 : float, default=1e-6
        Hyper-parameter : inverse scale parameter (rate parameter) for the
        Gamma distribution prior over the lambda parameter.

    compute_score : bool, default=False
        If True, compute the objective function at each step of the model.

    min_significance : float, default=0.5
        Minimum statistical significance (|beta|/sigma) required to keep a feature.
        Default of 0.5 provides a reasonable balance between feature selection
        and model accuracy.
        This replaces the threshold_lambda parameter for more interpretable feature pruning.
        
    standardize : bool, default=True
        Whether to standardize features before fitting. Recommended for
        consistent feature selection behavior regardless of feature scales.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    verbose : bool, default=False
        Verbose mode when fitting the model.

    Attributes
    ----------
    coef_ : array-like of shape (n_features,)
        Coefficients of the regression model (mean of distribution)

    alpha_ : float
       estimated precision of the noise.

    lambda_ : array-like of shape (n_features,)
       estimated precisions of the weights.

    sigma_ : array-like of shape (n_features, n_features)
        estimated variance-covariance matrix of the weights
        
    feature_significance_ : array-like of shape (n_features,)
        Statistical significance of each feature (|beta|/sigma).
        
    selected_features_ : array-like of shape (n_features,) of bool
        Boolean mask indicating which features were selected by the model.

    scores_ : float
        if computed, value of the objective function (to be maximized)

    n_iter_ : int
        The actual number of iterations to reach the stopping criterion.

        .. versionadded:: 1.3

    intercept_ : float
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    X_offset_ : float
        If `fit_intercept=True`, offset subtracted for centering data to a
        zero mean. Set to np.zeros(n_features) otherwise.

    X_scale_ : float
        Set to np.ones(n_features).

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    BayesianRidge : Bayesian ridge regression.

    Notes
    -----
    For an example, see :ref:`examples/linear_model/plot_ard.py
    <sphx_glr_auto_examples_linear_model_plot_ard.py>`.

    References
    ----------
    D. J. C. MacKay, Bayesian nonlinear modeling for the prediction
    competition, ASHRAE Transactions, 1994.

    R. Salakhutdinov, Lecture notes on Statistical Machine Learning,
    http://www.utstat.toronto.edu/~rsalakhu/sta4273/notes/Lecture2.pdf#page=15
    Their beta is our ``self.alpha_``
    Their alpha is our ``self.lambda_``
    ARD is a little different than the slide: only dimensions/features for
    which have statistical significance above ``self.min_significance`` are kept and the rest are
    discarded.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.ARDRegression()
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    ARDRegression()
    >>> clf.predict([[1, 1]])
    array([1.])
    """

    _parameter_constraints: dict = {
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "tol": [Interval(Real, 0, None, closed="left")],
        "alpha_1": [Interval(Real, 0, None, closed="left")],
        "alpha_2": [Interval(Real, 0, None, closed="left")],
        "lambda_1": [Interval(Real, 0, None, closed="left")],
        "lambda_2": [Interval(Real, 0, None, closed="left")],
        "compute_score": ["boolean"],
        "min_significance": [Interval(Real, 0, None, closed="left")],
        "standardize": ["boolean"],
        "fit_intercept": ["boolean"],
        "copy_X": ["boolean"],
        "verbose": ["boolean"],
    }

    def __init__(
        self,
        *,
        max_iter=300,
        tol=1.0e-3,
        alpha_1=1.0e-6,
        alpha_2=1.0e-6,
        lambda_1=1.0e-6,
        lambda_2=1.0e-6,
        compute_score=False,
        min_significance=0.5,
        standardize=True,
        fit_intercept=True,
        copy_X=True,
        verbose=False,
    ):
        self.max_iter = max_iter
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.alpha_1 = alpha_1
        self.alpha_2 = alpha_2
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.compute_score = compute_score
        self.min_significance = min_significance
        self.standardize = standardize
        self.copy_X = copy_X
        self.verbose = verbose

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y):
        """Fit the model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.
        y : array-like of shape (n_samples,)
            Target values (integers). Will be cast to X's dtype if necessary.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        X, y = validate_data(
            self,
            X,
            y,
            dtype=[np.float64, np.float32],
            force_writeable=True,
            y_numeric=True,
            ensure_min_samples=2,
        )
        dtype = X.dtype

        n_samples, n_features = X.shape
        
        # Store original y scale for proper predictions 
        self.y_mean_ = np.mean(y)
        self.y_scale_ = np.std(y)
        if self.y_scale_ == 0:
            self.y_scale_ = 1.0
            
        # Make a copy of original y
        y_orig = y.copy()

        # Apply standardization if requested
        if self.standardize:
            # Standardize X
            self.scaler_ = StandardScaler()
            X = self.scaler_.fit_transform(X)
            
            # Standardize y - this is crucial for numerical stability with large y values
            y_std = (y - self.y_mean_) / self.y_scale_
            y = y_std.copy()
            
            # Store information about feature scales for diagnostics
            self.feature_scales_ = self.scaler_.scale_
            eps = np.finfo(np.float64).eps
            self.max_feature_ratio_ = np.max(self.feature_scales_) / (np.min(self.feature_scales_) + eps)
            
            # Set flag for high response scaling
            self.high_response_scale_ = self.y_scale_ > 100

        # Initialize coefficients
        coef_ = np.zeros(n_features, dtype=dtype)

        # Handle preprocessing of data
        X, y, X_offset_, y_offset_, X_scale_ = _preprocess_data(
            X, y, fit_intercept=self.fit_intercept, copy=self.copy_X
        )

        self.X_offset_ = X_offset_
        self.X_scale_ = X_scale_

        # Launch the convergence loop
        keep_features = np.ones(n_features, dtype=bool)
        keep_features_old = keep_features.copy()

        # Initialize internal variables
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        alpha_1 = self.alpha_1
        alpha_2 = self.alpha_2
        verbose = self.verbose

        # Initialization of the values of the parameters
        eps = np.finfo(np.float64).eps
        
        # Initialize alpha (noise precision) with scale-appropriate value
        alpha_ = np.asarray(1.0 / (np.var(y) + eps), dtype=dtype)
        
        # Initialize lambda based on data scale to avoid numerical issues
        data_scale = np.mean(np.abs(X)) + eps
        lambda_ = np.ones(n_features, dtype=dtype) / (data_scale**2 + eps)

        self.scores_ = list()
        coef_old_ = None
        
        # Initialize adaptive significance threshold
        self._adaptive_threshold = self.min_significance

        def update_coeff(X, y, coef_, alpha_, keep_features, sigma_):
            coef_[keep_features] = alpha_ * np.linalg.multi_dot(
                [sigma_, X[:, keep_features].T, y]
            )
            return coef_

        update_sigma = (
            self._update_sigma
            if n_samples >= n_features
            else self._update_sigma_woodbury
        )
        
        # Initialize significance array
        significance = np.zeros(n_features, dtype=dtype)

        # Iterative procedure of ARD Regression
        for iter_ in range(self.max_iter):
            sigma_ = update_sigma(X, alpha_, lambda_, keep_features)
            coef_ = update_coeff(X, y, coef_, alpha_, keep_features, sigma_)

            # Update alpha and lambda
            residuals = y - np.dot(X[:, keep_features], coef_[keep_features])
            rmse_ = np.sum(residuals ** 2)
            
            # Protection against numerical issues
            if rmse_ < eps:
                rmse_ = eps
                
            # Calculate useful quantities
            gamma_ = np.maximum(0.0, np.minimum(1.0, 1.0 - lambda_[keep_features] * np.diag(sigma_)))
            
            # Update lambda (precisions of weight distributions)
            lambda_[keep_features] = (gamma_ + 2.0 * lambda_1) / (
                (coef_[keep_features]) ** 2 + 2.0 * lambda_2
            )
            
            # Prevent extreme lambda values that cause numerical issues
            lambda_ = np.clip(lambda_, 1e-10, 1e10)

            # Update alpha (precision of noise)
            gamma_sum = gamma_.sum()
            if gamma_sum >= n_samples:  # Prevent negative values
                gamma_sum = n_samples - 0.01
                
            alpha_ = (n_samples - gamma_sum + 2.0 * alpha_1) / (rmse_ + 2.0 * alpha_2)

            # Calculate statistical significance for feature pruning
            significance = np.zeros(n_features, dtype=dtype)
            
            # First compute raw significance values
            significance[keep_features] = np.abs(coef_[keep_features]) * np.sqrt(lambda_[keep_features])
            
            # Adaptive significance threshold based on data
            if iter_ == 0:
                # On first iteration, scale threshold based on the largest significance
                max_sig = np.max(significance) if np.any(significance > 0) else 1.0
                if max_sig > 1e3:
                    # With very large significance values, use a more aggressive threshold scaling
                    scaling_factor = np.log10(max_sig) / 10
                else:
                    scaling_factor = max(0.01, min(1.0, 1.0/max_sig))
                
                # For data with high response scaling, use a more lenient threshold
                if hasattr(self, 'high_response_scale_') and self.high_response_scale_:
                    # Decrease the threshold by 1/3 for high response scaling
                    scaling_factor *= 0.33
                
                self._adaptive_threshold = self.min_significance * scaling_factor
                if self.verbose:
                    print(f"Using adaptive significance threshold: {self._adaptive_threshold:.6f}")
                    if hasattr(self, 'max_feature_ratio_') and self.max_feature_ratio_ > 100:
                        print(f"High feature scale ratio detected: {self.max_feature_ratio_:.2f}")
                    if hasattr(self, 'high_response_scale_') and self.high_response_scale_:
                        print(f"High response scale detected: {self.y_scale_:.2f}")
            
            # Determine features to keep with a scale-aware approach
            keep_features_new = significance > self._adaptive_threshold
            
            # If no features would be selected, keep the top k most significant ones
            # OR if we have a highly scaled response with too few features, select more
            high_response_with_few_features = (
                hasattr(self, 'high_response_scale_') and 
                self.high_response_scale_ and 
                np.sum(keep_features_new) < min(5, int(n_features*0.1))
            )
            
            if (not np.any(keep_features_new) or high_response_with_few_features) and np.any(significance > 0):
                # Determine how many features to keep based on scaling context
                if hasattr(self, 'high_response_scale_') and self.high_response_scale_:
                    # For high response scaling, we need more features (at least 5 or 10% of total)
                    k = max(min(5, n_features), min(int(n_features*0.1), 10))
                    if self.verbose and high_response_with_few_features:
                        print(f"High response scaling detected - keeping {k} features for better prediction")
                else:
                    # Normal case - just keep a few top features
                    k = min(5, n_features)
                    
                top_k_idx = np.argsort(significance)[-k:]
                keep_features_new[top_k_idx] = True
                
                if self.verbose and not np.any(significance > self._adaptive_threshold):
                    print(f"No features above threshold. Keeping top {k} features.")

            # Compute the objective function
            if self.compute_score:
                s = (lambda_1 * np.log(lambda_ + eps) - lambda_2 * lambda_).sum()
                s += alpha_1 * log(alpha_ + eps) - alpha_2 * alpha_
                s += 0.5 * (
                    fast_logdet(sigma_)
                    + n_samples * log(alpha_ + eps)
                    + np.sum(np.log(lambda_[keep_features] + eps))
                )
                s -= 0.5 * (
                    alpha_ * rmse_ + (lambda_[keep_features] * coef_[keep_features]**2).sum()
                )
                self.scores_.append(s)

            # Check for convergence
            # Both coefficients and feature selection should stabilize
            coef_change = np.sum(np.abs(coef_old_ - coef_)) if coef_old_ is not None else np.inf
            features_change = np.sum(keep_features != keep_features_old)
            
            # If we have extreme scaling, use relative change instead of absolute
            if hasattr(self, 'y_scale_') and self.y_scale_ > 100:
                # For scaled data, use relative coefficient change
                if coef_old_ is not None and np.any(coef_old_ != 0):
                    nonzero_mask = np.abs(coef_old_) > eps
                    if np.any(nonzero_mask):
                        rel_change = np.mean(np.abs((coef_[nonzero_mask] - coef_old_[nonzero_mask]) / coef_old_[nonzero_mask]))
                        converged = rel_change < self.tol/10  # Need tighter convergence for scaled data
                    else:
                        converged = coef_change < self.tol
                else:
                    converged = False
            else:
                # Normal convergence check for non-scaled data
                converged = coef_change < self.tol
            
            if (iter_ > 0 and converged and features_change == 0):
                if verbose:
                    print(f"Converged after {iter_} iterations")
                break
                
            # Update for next iteration
            coef_old_ = np.copy(coef_)
            keep_features_old = keep_features.copy()
            coef_[~keep_features_new] = 0.0
            keep_features = keep_features_new

            if not keep_features.any():
                if verbose:
                    print("No features selected. Stopping.")
                break

        self.n_iter_ = iter_ + 1

        if keep_features.any():
            # update sigma and mu using updated params from the last iteration
            sigma_ = update_sigma(X, alpha_, lambda_, keep_features)
            coef_ = update_coeff(X, y, coef_, alpha_, keep_features, sigma_)
        else:
            sigma_ = np.array([]).reshape(0, 0)

        # Store internal variables for prediction
        self.y_offset_ = y_offset_
        
        # Store coefficients in standardized space
        self.coef_std_ = coef_.copy()
        
        # Store results
        self.alpha_ = alpha_
        self.sigma_ = sigma_
        self.lambda_ = lambda_
        
        # Store feature significance and selected features
        self.feature_significance_ = significance
        self.selected_features_ = keep_features
        
        # Transform coefficients back to original scale if standardization was applied
        if self.standardize:
            # Get the standard deviations used for scaling
            scale_factors = self.scaler_.scale_
            means = self.scaler_.mean_
            
            # Transform coefficients back to original scale
            # For features not selected, keep coefficient at zero
            self.coef_ = np.zeros_like(coef_)
            if keep_features.any():
                self.coef_[keep_features] = coef_[keep_features] / scale_factors[keep_features]
            
            # For the original scale intercept, we include the y_mean
            self._set_intercept(X_offset_, y_offset_, X_scale_)
            self.intercept_ = self.y_mean_
            if self.fit_intercept and keep_features.any():
                self.intercept_ -= np.sum(means[keep_features] * self.coef_[keep_features])
        else:
            self.coef_ = coef_
            self._set_intercept(X_offset_, y_offset_, X_scale_)
            
        # Final check - ensure predictions are accurate
        if keep_features.any() and self.standardize:
            # Make a test prediction on the training data
            y_pred = self.predict(X)
            mse = np.mean((y_pred - y_orig) ** 2)
            if self.verbose:
                print(f"Training MSE: {mse:.6f}")
                print(f"Features selected: {np.sum(keep_features)}/{n_features}")
                
                # Add additional verification for extreme scaling scenarios
                if self.standardize and self.y_scale_ > 100:
                    # For large scale differences, verify scaling is working
                    y_pred_small = self.predict(X[:5])
                    y_orig_small = y_orig[:5]
                    print(f"Scale verification - Original y range: [{np.min(y_orig_small):.2f}, {np.max(y_orig_small):.2f}]")
                    print(f"Scale verification - Predicted y range: [{np.min(y_pred_small):.2f}, {np.max(y_pred_small):.2f}]")

        return self

    def _update_sigma_woodbury(self, X, alpha_, lambda_, keep_features):
        # See slides as referenced in the docstring note
        # this function is used when n_samples < n_features and will invert
        # a matrix of shape (n_samples, n_samples) making use of the
        # woodbury formula:
        # https://en.wikipedia.org/wiki/Woodbury_matrix_identity
        n_samples = X.shape[0]
        X_keep = X[:, keep_features]
        inv_lambda = 1 / lambda_[keep_features].reshape(1, -1)
        sigma_ = pinvh(
            np.eye(n_samples, dtype=X.dtype) / alpha_
            + np.dot(X_keep * inv_lambda, X_keep.T)
        )
        sigma_ = np.dot(sigma_, X_keep * inv_lambda)
        sigma_ = -np.dot(inv_lambda.reshape(-1, 1) * X_keep.T, sigma_)
        sigma_[np.diag_indices(sigma_.shape[1])] += 1.0 / lambda_[keep_features]
        return sigma_

    def _update_sigma(self, X, alpha_, lambda_, keep_features):
        # See slides as referenced in the docstring note
        # this function is used when n_samples >= n_features and will
        # invert a matrix of shape (n_features, n_features)
        X_keep = X[:, keep_features]
        gram = np.dot(X_keep.T, X_keep)
        eye = np.eye(gram.shape[0], dtype=X.dtype)
        sigma_inv = lambda_[keep_features] * eye + alpha_ * gram
        sigma_ = pinvh(sigma_inv)
        return sigma_

    def predict(self, X, return_std=False):
        """Predict using the linear model.

        In addition to the mean of the predictive distribution, also its
        standard deviation can be returned.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples.

        return_std : bool, default=False
            Whether to return the standard deviation of posterior prediction.

        Returns
        -------
        y_mean : array-like of shape (n_samples,)
            Mean of predictive distribution of query points.

        y_std : array-like of shape (n_samples,)
            Standard deviation of predictive distribution of query points.
        """
        X = check_array(X)
        
        # Apply standardization if it was used during training
        if hasattr(self, 'scaler_') and self.standardize:
            X = self.scaler_.transform(X)
            
        # Get features that were retained during training
        if hasattr(self, 'selected_features_'):
            keep_features = self.selected_features_
        else:
            # Fall back to threshold-based selection for backward compatibility
            keep_features = self.lambda_ < getattr(self, 'threshold_lambda', 1e4)
        
        # Calculate prediction in the space used during training
        if not keep_features.any():
            # If no features are selected, predict the mean
            if self.fit_intercept:
                y_pred_std = np.full(X.shape[0], self.y_offset_)
            else:
                y_pred_std = np.zeros(X.shape[0])
        else:
            # Make prediction using appropriate coefficients
            if hasattr(self, 'coef_std_') and self.standardize:
                # Use standardized space coefficients if available
                y_pred_std = np.dot(X[:, keep_features], self.coef_std_[keep_features])
            else:
                # Otherwise use normal coefficients
                y_pred_std = np.dot(X[:, keep_features], self.coef_[keep_features])
            
            if self.fit_intercept:
                y_pred_std += self.y_offset_
        
        # Transform prediction back to original scale if standardization was applied
        if hasattr(self, 'y_scale_') and hasattr(self, 'y_mean_') and self.standardize:
            y_mean = y_pred_std * self.y_scale_ + self.y_mean_
        else:
            # Use the standard decision function if not standardized
            y_mean = y_pred_std
            
        if return_std is False:
            return y_mean
        else:
            # Calculate prediction uncertainty
            if not keep_features.any() or not hasattr(self, 'sigma_') or self.sigma_.size == 0:
                # If no features selected, uncertainty is just based on noise precision
                if hasattr(self, 'y_scale_') and self.standardize:
                    y_std = np.full(X.shape[0], self.y_scale_ / np.sqrt(self.alpha_))
                else:
                    y_std = np.full(X.shape[0], 1.0 / np.sqrt(self.alpha_))
            else:
                X_keep = X[:, keep_features]
                
                # Calculate uncertainty based on model covariance
                sigmas_squared_data = (np.dot(X_keep, self.sigma_) * X_keep).sum(axis=1)
                
                # Apply appropriate scaling for uncertainty
                if hasattr(self, 'y_scale_') and self.standardize:
                    y_std = self.y_scale_ * np.sqrt(sigmas_squared_data + (1.0 / self.alpha_))
                else:
                    y_std = np.sqrt(sigmas_squared_data + (1.0 / self.alpha_))
                    
            return y_mean, y_std
