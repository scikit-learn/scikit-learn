# Authors: David Dale dale.david@mail.ru
# License: TBD

import numpy as np
import warnings
from scipy import optimize, sparse

from ..base import BaseEstimator, RegressorMixin
from .base import LinearModel
from ..utils import check_X_y
from ..utils import check_consistent_length
from ..utils import axis0_safe_slice
from ..utils.extmath import safe_sparse_dot


def _quantile_loss_and_gradient(w, X, y, quantile, alpha, l1_ratio, sample_weight):
    """Returns the quantile regression loss and its gradient.

    Parameters
    ----------
    w : ndarray, shape (n_features + 0,) or (n_features + 1,)
        Feature vector.
        w[:n_features] gives the coefficients
        w[-1] gives the intercept factor.

    X : ndarray, shape (n_samples, n_features)
        Input data.

    y : ndarray, shape (n_samples,)
        Target vector.
    
    quantile : float
        Quantile to be predicted.
    
    alpha : float
        ElasticNet regularization parameter.
        
    l1_ratio: float
        Ratio of L1 vs L2 in regularization.

    sample_weight : ndarray, shape (n_samples,)
        Weight assigned to each sample.

    Returns
    -------
    loss : float
        Quantile loss.

    gradient : ndarray, shape (len(w))
        Returns the derivative of the Quantile loss with respect to each
        coefficient and intercept as a vector.
    """
    _, n_features = X.shape
    fit_intercept = (n_features + 1 == w.shape[0])
    if fit_intercept:
        intercept = w[-1]
    else:
        intercept = 0  # regardless of len(w)
    w = w[:n_features]

    # Calculate skewed absolute loss
    linear_loss = y - safe_sparse_dot(X, w)
    if fit_intercept:
        linear_loss -= intercept
    # positive_error = linear_loss > 0
    negative_error = linear_loss < 0
    # abs_error = np.abs(linear_loss)
    
    weighted_obs = (quantile - negative_error) * sample_weight
    # regression_loss = abs_error * (positive_error * quantile + negative_error * (1 - quantile))
    regression_loss = linear_loss * weighted_obs

    if fit_intercept:
        grad = np.zeros(n_features + 1)
    else:
        grad = np.zeros(n_features + 0)

    # Gradient due to the linear loss.
    grad[:n_features] -= safe_sparse_dot(weighted_obs, X)

    # Gradient due to the ridge penalty
    grad[:n_features] += alpha * (1-l1_ratio) * 2. * w
    # Gradient due to the lasso penalty
    grad[:n_features] += alpha * l1_ratio * np.sign(w)

    if fit_intercept:
        grad[-1] -= np.sum(weighted_obs)

    loss = np.sum(regression_loss) + alpha * (1-l1_ratio) * np.dot(w, w) + alpha * l1_ratio * np.sum(np.abs(w))
    return loss, grad

# Todo: make lasso precise enough, because now there are many coefficients around zero, but not exactly. Coo descent?


def _smooth_quantile_loss_and_gradient(w, X, y, quantile, alpha, l1_ratio, sample_weight, tau=0):
    """ Smooth approximation to quantile regression loss, gradient and hessian.
    Main loss and l1 penalty are both approximated by the same trick from Chen & Wei, 2005
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
    positive_error = linear_loss > quantile * tau
    negative_error = linear_loss < (quantile - 1) * tau
    small_error = ~ (positive_error | negative_error)

    # Calculate loss due to regression error
    regression_loss = (positive_error * (linear_loss * quantile - 0.5 * quantile**2 * tau)
                       + small_error * 0.5 * linear_loss**2 / (tau if tau != 0 else 1)  # Here the article LIES!
                       + negative_error * (linear_loss * (quantile-1) - 0.5 * (quantile-1)**2 * tau)
                       ) * sample_weight
    loss = np.sum(regression_loss)

    if fit_intercept:
        grad = np.zeros(n_features + 1)
    else:
        grad = np.zeros(n_features + 0)

    # Gradient due to the regression error
    weighted_grad = (positive_error * quantile
                     + small_error * linear_loss / (tau if tau != 0 else 1)  # Here the article LIES!
                     + negative_error * (quantile-1)) * sample_weight
    grad[:n_features] -= safe_sparse_dot(weighted_grad, X)

    if fit_intercept:
        grad[-1] -= np.sum(weighted_grad)

    # Gradient and loss due to the ridge penalty
    grad[:n_features] += alpha * (1 - l1_ratio) * 2. * w
    loss += alpha * (1 - l1_ratio) * np.dot(w, w)

    # Gradient and loss due to the lasso penalty
    # for smoothness, replace abs(w) with w^2/(2*tau)+tau/2 for abs(w)<tau
    if tau > 0:
        large_coef = np.abs(w) > tau
        small_coef = ~large_coef
        loss += alpha * l1_ratio * np.sum(large_coef * np.abs(w) + small_coef * (w ** 2 / (2 * tau) + tau / 2))
        grad[:n_features] += alpha * l1_ratio * (large_coef * np.sign(w) + small_coef * w / tau)
    else:
        loss += alpha * l1_ratio * np.sum(np.abs(w))
        grad[:n_features] += alpha * l1_ratio * np.sign(w)

    return loss, grad


class QuantileRegressor(LinearModel, RegressorMixin, BaseEstimator):
    """Linear regression model that is robust to outliers.

    The Quantile Regressor optimizes the skewed absolute loss 
     ``(y - X'w) (q-[y-X'w<0])``, where q is the desired quantile.

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
        
    l1_ratio : float, defaut 0.0
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

    tol : float, default 1e-5
        The iteration will stop when
        ``max{|proj g_i | i = 1, ..., n}`` <= ``tol``
        where pg_i is the i-th component of the projected gradient.

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
    .. [1] Chen, C., & Wei, Y. (2005). Computational issues for quantile regression. 
           Sankhya: The Indian Journal of Statistics, 399-417.
    """

    def __init__(self, quantile=0.5, max_iter=1000, alpha=0.0001, l1_ratio=0.0,
                 warm_start=False, fit_intercept=True, gtol=1e-4, first_tau=1e-2, xtol=1e-6):
        self.quantile = quantile
        self.max_iter = max_iter
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.warm_start = warm_start
        self.fit_intercept = fit_intercept
        self.gtol = gtol
        self.xtol = xtol
        self.first_tau = first_tau

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

        # solve sequence of optimization problems with different smoothing parameter
        total_iter = []
        loss_args = X, y, self.quantile, self.alpha, self.l1_ratio, sample_weight
        for i in range(10):
            tau = self.first_tau * 0.1**i
            result = optimize.minimize(
                _smooth_quantile_loss_and_gradient,
                parameters,
                args=loss_args + (tau, ),
                method='BFGS',
                jac=True,
                options={
                    # Gradient norm must be less than gtol before successful termination, but it is never so
                    # todo: fix convergence problem (invent criteria of true convergence)
                    'gtol': self.gtol,  # for 'BFGS'
                    # 'xtol': self.tol,  # for 'Newton-CG'
                    'maxiter': self.max_iter - sum(total_iter),
                }
                )
            total_iter.append(result['nit'])
            prev_parameters = parameters
            parameters = result['x']

            # for lasso, replace parameters with exact zero, if it increases likelihood
            if self.alpha * self.l1_ratio > 0:
                value, _ = _smooth_quantile_loss_and_gradient(parameters, *loss_args, tau=0)
                for j in range(len(parameters)):
                    new_parameters = parameters.copy()
                    new_parameters[j] = 0
                    new_value, _ = _smooth_quantile_loss_and_gradient(new_parameters, *loss_args, tau=0)
                    if new_value <= value:
                        value = new_value
                        parameters = new_parameters

            # stop if solution does not change between subproblems
            if np.linalg.norm(prev_parameters-parameters) < self.xtol:
                break
            # stop if maximum number of iterations is exceeded
            if sum(total_iter) >= self.max_iter:
                break
            # stop if tau is already zero
            if tau == 0:
                break
        # do I really need to issue this warning??? Its reason is lineSearchError, which cannot be easily fixed
        if not result['success']:
            warnings.warn("QuantileRegressor convergence failed:" +
                          " Scipy solver terminated with %s." % result['message']
                          )
        
        self.n_iter_ = sum(total_iter)
        self.tau_ = tau
        self.total_iter_ = total_iter
        if self.fit_intercept:
            self.intercept_ = parameters[-1]
        else:
            self.intercept_ = 0.0
        self.coef_ = parameters[:X.shape[1]]
        return self
