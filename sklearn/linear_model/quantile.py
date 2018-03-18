# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause

import numpy as np
import warnings
from scipy import optimize

from ..base import BaseEstimator, RegressorMixin
from .base import LinearModel
from ..utils import check_X_y
from ..utils import check_consistent_length
from ..utils.extmath import safe_sparse_dot


def _smooth_quantile_loss_and_gradient(
        w, X, y, quantile, alpha, l1_ratio, sample_weight, gamma=0):
    """ Smooth approximation to quantile regression loss, gradient and hessian.
    Main loss and l1 penalty are both approximated by the same trick
    from Chen & Wei, 2005
    """
    _, n_features = X.shape
    fit_intercept = (n_features + 1 == w.shape[0])
    if fit_intercept:
        intercept = w[-1]
    else:
        intercept = 0  # regardless of len(w)
    w = w[:n_features]

    # Discriminate positive, negative and small residuals
    linear_loss = y - safe_sparse_dot(X, w)
    if fit_intercept:
        linear_loss -= intercept
    positive_error = linear_loss > quantile * gamma
    negative_error = linear_loss < (quantile - 1) * gamma
    small_error = ~ (positive_error | negative_error)

    # Calculate loss due to regression error
    regression_loss = (
        positive_error * (linear_loss*quantile - 0.5*gamma*quantile**2) +
        small_error * 0.5*linear_loss**2 / (gamma if gamma != 0 else 1) +
        negative_error * (linear_loss*(quantile-1) - 0.5*gamma*(quantile-1)**2)
        ) * sample_weight
    loss = np.sum(regression_loss)

    if fit_intercept:
        grad = np.zeros(n_features + 1)
    else:
        grad = np.zeros(n_features + 0)

    # Gradient due to the regression error
    weighted_grad = (positive_error * quantile +
                     small_error * linear_loss / (gamma if gamma != 0 else 1) +
                     negative_error * (quantile-1)) * sample_weight
    grad[:n_features] -= safe_sparse_dot(weighted_grad, X)

    if fit_intercept:
        grad[-1] -= np.sum(weighted_grad)

    # Gradient and loss due to the ridge penalty
    grad[:n_features] += alpha * (1 - l1_ratio) * 2. * w
    loss += alpha * (1 - l1_ratio) * np.dot(w, w)

    # Gradient and loss due to the lasso penalty
    # for smoothness replace abs(w) with w^2/(2*gamma)+gamma/2 for abs(w)<gamma
    if gamma > 0:
        large_coef = np.abs(w) > gamma
        small_coef = ~large_coef
        loss += alpha*l1_ratio*np.sum(large_coef*np.abs(w) +
                                      small_coef*(w**2/(2*gamma) + gamma/2))
        grad[:n_features] += alpha*l1_ratio*(large_coef*np.sign(w) +
                                             small_coef*w/gamma)
    else:
        loss += alpha * l1_ratio * np.sum(np.abs(w))
        grad[:n_features] += alpha * l1_ratio * np.sign(w)

    return loss, grad


class QuantileRegressor(LinearModel, RegressorMixin, BaseEstimator):
    """Linear regression model that is robust to outliers.

    The Quantile Regressor optimizes the skewed absolute loss
     ``(y - X'w) (q - [y - X'w < 0])``, where q is the desired quantile.

    Optimization is performed as a sequence of smooth optimization problems.

    Read more in the :ref:`User Guide <quantile_regression>`

    .. versionadded:: 0.20

    Parameters
    ----------
    quantile : float, strictly between 0.0 and 1.0, default 0.5
        The quantile that the model predicts.

    max_iter : int, default 100
        Maximum number of iterations that scipy.optimize.minimize
        should run for.

    alpha : float, default 0.0001
        Constant that multiplies ElasticNet penalty term.

    l1_ratio : float, default 0.0
        The ElasticNet mixing parameter, with ``0 <= l1_ratio <= 1``. For
        ``l1_ratio = 0`` the penalty is an L2 penalty. ``For l1_ratio = 1`` it
        is an L1 penalty.  For ``0 < l1_ratio < 1``, the penalty is a
        combination of L1 and L2.

    warm_start : bool, default False
        This is useful if the stored attributes of a previously used model
        has to be reused. If set to False, then the coefficients will
        be rewritten for every call to fit.

    fit_intercept : bool, default True
        Whether or not to fit the intercept. This can be set to False
        if the data is already centered around the origin.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    normalize : boolean, optional, default False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.

    gamma : float, default 1e-2
        Starting value for smooth approximation.
        Absolute loss is replaced with quadratic for ``|error| < gamma``.
        Lasso penalty is replaced with quadratic for ``|w| < gamma``.
        ``gamma = 0`` gives exact non-smooth loss function.
        The algorithm performs consecutive optimizations with gamma
        decreasing by factor of ``gamma_decrease``,
        until ``xtol`` criterion is met,
        or until ``max_iter`` is exceeded.

    gamma_decrease: float, default 0.1
        The factor by which ``gamma`` is multiplied at each iteration.

    n_gamma_decreases: int, default 10
        Maximal number of iterations of approximation of the cost function.
        At each iteration, ``gamma`` is multiplied by a factor
        of ``gamma_decrease``

    gtol : float, default 1e-4
        The smooth optimizing iteration will stop when
        ``max{|proj g_i | i = 1, ..., n}`` <= ``gtol``
        where pg_i is the i-th component of the projected gradient.

    xtol : float, default 1e-6
        Global optimization will stop when ``|w_{t-1} - w_t|`` < ``xtol``
        where w_t is result of t'th approximated optimization.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Features got by optimizing the Huber loss.

    intercept_ : float
        Bias.

    n_iter_ : int
        Number of iterations that scipy.optimize.mimimize has run for.

    References
    ----------
    .. [1] Koenker, R., & Bassett Jr, G. (1978). Regression quantiles.
            Econometrica: journal of the Econometric Society, 33-50.

    .. [2] Chen, C., & Wei, Y. (2005).
           Computational issues for quantile regression.
           Sankhya: The Indian Journal of Statistics, 399-417.
    """

    def __init__(self, quantile=0.5, max_iter=1000, alpha=0.0001, l1_ratio=0.0,
                 warm_start=False, fit_intercept=True,
                 normalize=False, copy_X=True,
                 gamma=1e-2, gtol=1e-4, xtol=1e-6,
                 gamma_decrease=0.1, n_gamma_decreases=100):
        self.quantile = quantile
        self.max_iter = max_iter
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.warm_start = warm_start
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.normalize = normalize
        self.gtol = gtol
        self.xtol = xtol
        self.gamma = gamma
        self.gamma_decrease = gamma_decrease
        self.n_gamma_decreases = n_gamma_decreases

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like, shape (n_samples,)
            Weight given to each sample.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(
            X, y, copy=False, accept_sparse=['csr'], y_numeric=True)

        X, y, X_offset, y_offset, X_scale = self._preprocess_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        if sample_weight is not None:
            sample_weight = np.array(sample_weight)
            check_consistent_length(y, sample_weight)
        else:
            sample_weight = np.ones_like(y)

        if self.quantile >= 1.0 or self.quantile <= 0.0:
            raise ValueError(
                "Quantile should be strictly between 0.0 and 1.0, got %f"
                % self.quantile)

        if self.warm_start and hasattr(self, 'coef_'):
            parameters = np.concatenate(
                (self.coef_, [self.intercept_]))
        else:
            if self.fit_intercept:
                parameters = np.zeros(X.shape[1] + 1)
            else:
                parameters = np.zeros(X.shape[1] + 0)

        # solve sequence of optimization problems
        # with different smoothing parameter
        total_iter = []
        loss_args = (X, y, self.quantile, self.alpha, self.l1_ratio,
                     sample_weight)
        for i in range(self.n_gamma_decreases):
            gamma = self.gamma * self.gamma_decrease ** i
            result = optimize.minimize(
                _smooth_quantile_loss_and_gradient,
                parameters,
                args=loss_args + (gamma, ),
                method='L-BFGS-B',
                jac=True,
                options={
                    'gtol': self.gtol,
                    'maxiter': self.max_iter - sum(total_iter),
                }
                )
            total_iter.append(result['nit'])
            prev_parameters = parameters
            parameters = result['x']

            # for lasso, replace parameters with exact zero,
            # if this decreases the cost function
            if self.alpha * self.l1_ratio > 0:
                value, _ = _smooth_quantile_loss_and_gradient(parameters,
                                                              *loss_args,
                                                              gamma=0)
                for j in range(len(parameters)):
                    new_parameters = parameters.copy()
                    old_param = new_parameters[j]
                    new_parameters[j] = 0
                    new_value, _ = _smooth_quantile_loss_and_gradient(
                        new_parameters, *loss_args, gamma=0)
                    # check if the cost function decreases,
                    # or increases, but by little, and param is small anyway
                    if new_value <= value \
                            or np.abs(old_param) < self.xtol \
                            and new_value < value + self.gtol:
                        value = new_value
                        parameters = new_parameters

            # stop if solution does not change between subproblems
            if np.linalg.norm(prev_parameters-parameters) < self.xtol:
                break
            # stop if maximum number of iterations is exceeded
            if sum(total_iter) >= self.max_iter:
                break
            # stop if gamma is already zero
            if gamma == 0:
                break
        # do I really need to issue this warning?
        # Its reason is lineSearchError, which cannot be easily fixed
        if not result['success']:
            warnings.warn("QuantileRegressor convergence failed:" +
                          " Scipy solver terminated with %s."
                          % result['message']
                          )
        self.n_iter_ = sum(total_iter)
        self.gamma_ = gamma
        self.total_iter_ = total_iter
        self.coef_ = parameters[:X.shape[1]]
        # do not use self.set_intercept_, because it assumes intercept is zero
        # if the data is normalized, which is false in this case
        if self.fit_intercept:
            self.coef_ = self.coef_ / X_scale
            self.intercept_ = parameters[-1] + y_offset \
                              - np.dot(X_offset, self.coef_.T)
        else:
            self.intercept_ = 0.0
        return self
