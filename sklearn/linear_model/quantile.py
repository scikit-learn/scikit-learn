# Authors: David Dale dale.david@mail.ru
# License: TBD

import numpy as np

from scipy import optimize, sparse

from ..base import BaseEstimator, RegressorMixin
from .base import LinearModel
from ..utils import check_X_y
from ..utils import check_consistent_length
from ..utils import axis0_safe_slice
from ..utils.extmath import safe_sparse_dot


def _quantile_loss_and_gradient(w, X, y, quantile, alpha, sample_weight):
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
        Quantile to be predicted
    
    alpha : float
        Ridge regularization parameter.

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
    X_is_sparse = sparse.issparse(X)
    _, n_features = X.shape
    fit_intercept = (n_features + 1 == w.shape[0])
    if fit_intercept:
        intercept = w[-1]
    w = w[:n_features]
    n_samples = np.sum(sample_weight)

    # Calculate skewed absolute loss
    linear_loss = y - safe_sparse_dot(X, w)
    if fit_intercept:
        linear_loss -= intercept
    #positive_error = linear_loss > 0
    negative_error = linear_loss < 0
    #abs_error = np.abs(linear_loss)
    
    weighted_obs = (quantile - negative_error) * sample_weight
    #regression_loss = abs_error * (positive_error * quantile + negative_error * (1 - quantile))
    regression_loss = linear_loss * weighted_obs

    if fit_intercept:
        grad = np.zeros(n_features + 1)
    else:
        grad = np.zeros(n_features + 0)

    # Gradient due to the linear loss.
    grad[:n_features] -= safe_sparse_dot(weighted_obs, X)

    # Gradient due to the ridge penalty
    grad[:n_features] += alpha * 2. * w

    if fit_intercept:
        grad[-1] -= np.sum(weighted_obs)

    loss = np.sum(regression_loss) + alpha * np.dot(w, w)
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
        Ridge regularization parameter.

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
           Sankhyā: The Indian Journal of Statistics, 399-417.
    """

    def __init__(self, quantile=0.5, max_iter=100, alpha=0.0001,
                 warm_start=False, fit_intercept=True, tol=1e-05):
        self.quantile = quantile
        self.max_iter = max_iter
        self.alpha = alpha
        self.warm_start = warm_start
        self.fit_intercept = fit_intercept
        self.tol = tol

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

        # Type Error caused in old versions of SciPy because of no
        # maxiter argument ( <= 0.9).
        result = optimize.minimize(
            _quantile_loss_and_gradient, 
            parameters,
            args=(X, y, self.quantile, self.alpha, sample_weight),
            #maxiter=self.max_iter, 
            #gtol=self.tol,
            method='BFGS',
            #iprint=0,
            jac = True,
            )
        # ToDo: issue warning instead of ValueError
        if result['success'] == False:
            raise ValueError("QuantileRegressor convergence failed:"
                             " Scipy solver terminated with %s"
                             % result['message'])
        
        self.n_iter_ = result.get('nit', None)
        if self.fit_intercept:
            self.intercept_ = result['x'][-1]
        else:
            self.intercept_ = 0.0
        self.coef_ = result['x'][:X.shape[1]]
        return self
