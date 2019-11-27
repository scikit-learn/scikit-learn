"""
Generalized Linear Models with Exponential Dispersion Family
"""

# Author: Christian Lorentzen <lorentzen.ch@googlemail.com>
# some parts and tricks stolen from other sklearn files.
# License: BSD 3 clause

# TODO: Add cross validation support, e.g. GCV?
# TODO: Should GeneralizedLinearRegressor inherit from LinearModel?
#       So far, it does not.
# TODO: Include further classes in class.rst? ExponentialDispersionModel?
#       TweedieDistribution?
# TODO: Negative values in P1 are not allowed so far. They could be used
#       for group lasso.

# Design Decisions:
# - Which name? GeneralizedLinearModel vs GeneralizedLinearRegressor.
#   Estimators in sklearn are either regressors or classifiers. A GLM can do
#   both depending on the distr (Normal => regressor, Binomial => classifier).
#   Solution: GeneralizedLinearRegressor since this is the focus.
# - Allow for finer control of penalty terms:
#   L1: ||P1*w||_1 with P1*w as element-wise product, this allows to exclude
#       factors from the L1 penalty.
#   L2: w*P2*w with P2 a positive (semi-) definite matrix, e.g. P2 could be
#   a 1st or 2nd order difference matrix (compare B-spline penalties and
#   Tikhonov regularization).
# - The link function (instance of class Link) is necessary for the evaluation
#   of deviance, score, Fisher and Hessian matrix as functions of the
#   coefficients, which is needed by optimizers.
#   Solution: link as argument in those functions
# - Which name/symbol for sample_weight in docu?
#   sklearn.linear_models uses w for coefficients, standard literature on
#   GLMs use beta for coefficients and w for (sample) weights.
#   So far, coefficients=w and sample weights=s.
# - The intercept term is the first index, i.e. coef[0]


from __future__ import division
from abc import ABCMeta, abstractmethod
import numbers
import numpy as np
from scipy import linalg, sparse, special
import scipy.sparse.linalg as splinalg
from scipy.optimize import fmin_l_bfgs_b
import warnings
from ..base import BaseEstimator, RegressorMixin
from ..exceptions import ConvergenceWarning
from ..utils import check_array, check_X_y
from ..utils.optimize import newton_cg
from ..utils.validation import check_is_fitted, check_random_state


def _check_weights(sample_weight, n_samples):
    """Check that sample weights are non-negative and have the right shape."""
    if sample_weight is None:
        weights = np.ones(n_samples)
    elif np.isscalar(sample_weight):
        if sample_weight <= 0:
            raise ValueError("Sample weights must be non-negative.")
        weights = sample_weight * np.ones(n_samples)
    else:
        _dtype = [np.float64, np.float32]
        weights = check_array(sample_weight, accept_sparse=False,
                              force_all_finite=True, ensure_2d=False,
                              dtype=_dtype)
        if weights.ndim > 1:
            raise ValueError("Sample weight must be 1D array or scalar")
        elif weights.shape[0] != n_samples:
            raise ValueError("Sample weights must have the same length as "
                             "y")
        if not np.all(weights >= 0):
            raise ValueError("Sample weights must be non-negative.")
        elif not np.sum(weights) > 0:
            raise ValueError("Sample weights must have at least one positive "
                             "element.")

    return weights


def _safe_lin_pred(X, coef):
    """Compute the linear predictor taking care if intercept is present."""
    if coef.size == X.shape[1] + 1:
        return X @ coef[1:] + coef[0]
    else:
        return X @ coef


def _safe_toarray(X):
    """Returns a numpy array."""
    if sparse.issparse(X):
        return X.toarray()
    else:
        return np.asarray(X)


def _safe_sandwich_dot(X, d, intercept=False):
    """Compute sandwich product X.T @ diag(d) @ X.

    With ``intercept=True``, X is treated as if a column of 1 were appended as
    first column of X.
    X can be sparse, d must be an ndarray. Always returns a ndarray."""
    if sparse.issparse(X):
        temp = (X.transpose() @ X.multiply(d[:, np.newaxis]))
        # for older versions of numpy and scipy, temp may be a np.matrix
        temp = _safe_toarray(temp)
    else:
        temp = (X.T * d) @ X
    if intercept:
        dim = X.shape[1] + 1
        if sparse.issparse(X):
            order = 'F' if sparse.isspmatrix_csc(X) else 'C'
        else:
            order = 'F' if X.flags['F_CONTIGUOUS'] else 'C'
        res = np.empty((dim, dim), dtype=max(X.dtype, d.dtype), order=order)
        res[0, 0] = d.sum()
        res[1:, 0] = d @ X
        res[0, 1:] = res[1:, 0]
        res[1:, 1:] = temp
    else:
        res = temp
    return res


def _min_norm_sugrad(coef, grad, P2, P1):
    """Compute the gradient of all subgradients with minimal L2-norm.

    subgrad = grad + P2 * coef + P1 * subgrad(|coef|_1)

    g_i = grad_i + (P2*coef)_i

    if coef_i > 0:   g_i + P1_i
    if coef_i < 0:   g_i - P1_i
    if coef_i = 0:   sign(g_i) * max(|g_i|-P1_i, 0)

    Parameters
    ----------
    coef : ndarray
        coef[0] may be intercept.

    grad : ndarray, shape=coef.shape

    P2 : {1d or 2d array, None}
        always without intercept, ``None`` means P2 = 0

    P1 : ndarray
        always without intercept
    """
    intercept = (coef.size == P1.size + 1)
    idx = 1 if intercept else 0  # offset if coef[0] is intercept
    # compute grad + coef @ P2 without intercept
    grad_wP2 = grad[idx:].copy()
    if P2 is None:
        pass
    elif P2.ndim == 1:
        grad_wP2 += coef[idx:] * P2
    else:
        grad_wP2 += coef[idx:] @ P2
    res = np.where(coef[idx:] == 0,
                   np.sign(grad_wP2) * np.maximum(np.abs(grad_wP2) - P1, 0),
                   grad_wP2 + np.sign(coef[idx:]) * P1)
    if intercept:
        return np.concatenate(([grad[0]], res))
    else:
        return res


class Link(metaclass=ABCMeta):
    """Abstract base class for Link functions."""

    @abstractmethod
    def link(self, mu):
        """Compute the link function g(mu).

        The link function links the mean mu=E[Y] to the so called linear
        predictor (X*w), i.e. g(mu) = linear predictor.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Usually the (predicted) mean.
        """
        pass

    @abstractmethod
    def derivative(self, mu):
        """Compute the derivative of the link g'(mu).

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Usually the (predicted) mean.
        """
        pass

    @abstractmethod
    def inverse(self, lin_pred):
        """Compute the inverse link function h(lin_pred).

        Gives the inverse relationship between linear predictor and the mean
        mu=E[Y], i.e. h(linear predictor) = mu.

        Parameters
        ----------
        lin_pred : array, shape (n_samples,)
            Usually the (fitted) linear predictor.
        """
        pass

    @abstractmethod
    def inverse_derivative(self, lin_pred):
        """Compute the derivative of the inverse link function h'(lin_pred).

        Parameters
        ----------
        lin_pred : array, shape (n_samples,)
            Usually the (fitted) linear predictor.
        """
        pass

    @abstractmethod
    def inverse_derivative2(self, lin_pred):
        """Compute 2nd derivative of the inverse link function h''(lin_pred).

        Parameters
        ----------
        lin_pred : array, shape (n_samples,)
            Usually the (fitted) linear predictor.
        """
        pass


class IdentityLink(Link):
    """The identity link function g(x)=x."""

    def link(self, mu):
        return mu

    def derivative(self, mu):
        return np.ones_like(mu)

    def inverse(self, lin_pred):
        return lin_pred

    def inverse_derivative(self, lin_pred):
        return np.ones_like(lin_pred)

    def inverse_derivative2(self, lin_pred):
        return np.zeros_like(lin_pred)


class LogLink(Link):
    """The log link function g(x)=log(x)."""

    def link(self, mu):
        return np.log(mu)

    def derivative(self, mu):
        return 1./mu

    def inverse(self, lin_pred):
        return np.exp(lin_pred)

    def inverse_derivative(self, lin_pred):
        return np.exp(lin_pred)

    def inverse_derivative2(self, lin_pred):
        return np.exp(lin_pred)


class LogitLink(Link):
    """The logit link function g(x)=logit(x)."""

    def link(self, mu):
        return special.logit(mu)

    def derivative(self, mu):
        return 1. / (mu * (1 - mu))

    def inverse(self, lin_pred):
        return special.expit(lin_pred)

    def inverse_derivative(self, lin_pred):
        ep = special.expit(lin_pred)
        return ep * (1. - ep)

    def inverse_derivative2(self, lin_pred):
        ep = special.expit(lin_pred)
        return ep * (1. - ep) * (1. - 2 * ep)


class ExponentialDispersionModel(metaclass=ABCMeta):
    r"""Base class for reproductive Exponential Dispersion Models (EDM).

    The pdf of :math:`Y\sim \mathrm{EDM}(\mu, \phi)` is given by

    .. math:: p(y| \theta, \phi) = c(y, \phi)
        \exp\left(\frac{\theta y-A(\theta)}{\phi}\right)
        = \tilde{c}(y, \phi)
            \exp\left(-\frac{d(y, \mu)}{2\phi}\right)

    with mean :math:`\mathrm{E}[Y] = A'(\theta) = \mu`,
    variance :math:`\mathrm{Var}[Y] = \phi \cdot v(\mu)`,
    unit variance :math:`v(\mu)` and
    unit deviance :math:`d(y,\mu)`.

    Attributes
    ----------
    lower_bound
    upper_bound
    include_lower_bound
    include_upper_bound

    Methods
    -------
    in_y_range
    unit_variance
    unit_variance_derivative
    variance
    variance_derivative
    unit_deviance
    unit_deviance_derivative
    deviance
    deviance_derivative
    starting_mu

    _mu_deviance_derivative
    _score
    _fisher_matrix
    _observed_information
    _eta_mu_score_fisher

    References
    ----------

    https://en.wikipedia.org/wiki/Exponential_dispersion_model.
    """
    @property
    def lower_bound(self):
        """Get the lower bound of values for Y~EDM."""
        return self._lower_bound

    @property
    def upper_bound(self):
        """Get the upper bound of values for Y~EDM."""
        return self._upper_bound

    @property
    def include_lower_bound(self):
        """Get True if lower bound for y is included: y >= lower_bound."""
        return self._include_lower_bound

    @property
    def include_upper_bound(self):
        """Get True if upper bound for y is included: y <= upper_bound."""
        return self._include_upper_bound

    def in_y_range(self, x):
        """Returns ``True`` if x is in the valid range of Y~EDM.

        Parameters
        ----------
        x : array, shape (n_samples,)
            Target values.
        """
        if self.include_lower_bound:
            if self.include_upper_bound:
                return np.logical_and(np.greater_equal(x, self.lower_bound),
                                      np.less_equal(x, self.upper_bound))
            else:
                return np.logical_and(np.greater_equal(x, self.lower_bound),
                                      np.less(x, self.upper_bound))
        else:
            if self.include_upper_bound:
                return np.logical_and(np.greater(x, self.lower_bound),
                                      np.less_equal(x, self.upper_bound))
            else:
                return np.logical_and(np.greater(x, self.lower_bound),
                                      np.less(x, self.upper_bound))

    @abstractmethod
    def unit_variance(self, mu):
        r"""Compute the unit variance function.

        The unit variance :math:`v(\mu)` determines the variance as
        a function of the mean :math:`\mu` by
        :math:`\mathrm{Var}[Y_i] = \phi/s_i*v(\mu_i)`.
        It can also be derived from the unit deviance :math:`d(y,\mu)` as

        .. math:: v(\mu) = \frac{2}{\frac{\partial^2 d(y,\mu)}{
            \partial\mu^2}}\big|_{y=\mu}

        See also :func:`variance`.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.
        """
        pass

    @abstractmethod
    def unit_variance_derivative(self, mu):
        r"""Compute the derivative of the unit variance w.r.t. mu.

        Return :math:`v'(\mu)`.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Target values.
        """
        pass

    def variance(self, mu, phi=1, weights=1):
        r"""Compute the variance function.

        The variance of :math:`Y_i \sim \mathrm{EDM}(\mu_i,\phi/s_i)` is
        :math:`\mathrm{Var}[Y_i]=\phi/s_i*v(\mu_i)`,
        with unit variance :math:`v(\mu)` and weights :math:`s_i`.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.

        phi : float (default=1)
            Dispersion parameter.

        weights : array, shape (n_samples,) (default=1)
            Weights or exposure to which variance is inverse proportional.
        """
        return phi/weights * self.unit_variance(mu)

    def variance_derivative(self, mu, phi=1, weights=1):
        r"""Compute the derivative of the variance w.r.t. mu.

        Returns
        :math:`\frac{\partial}{\partial\mu}\mathrm{Var}[Y_i]
        =phi/s_i*v'(\mu_i)`, with unit variance :math:`v(\mu)`
        and weights :math:`s_i`.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.

        phi : float (default=1)
            Dispersion parameter.

        weights : array, shape (n_samples,) (default=1)
            Weights or exposure to which variance is inverse proportional.
        """
        return phi/weights * self.unit_variance_derivative(mu)

    @abstractmethod
    def unit_deviance(self, y, mu):
        r"""Compute the unit deviance.

        The unit_deviance :math:`d(y,\mu)` can be defined by the
        log-likelihood as
        :math:`d(y,\mu) = -2\phi\cdot
        \left(loglike(y,\mu,\phi) - loglike(y,y,\phi)\right).`

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.
        """
        pass

    def unit_deviance_derivative(self, y, mu):
        r"""Compute the derivative of the unit deviance w.r.t. mu.

        The derivative of the unit deviance is given by
        :math:`\frac{\partial}{\partial\mu}d(y,\mu) = -2\frac{y-\mu}{v(\mu)}`
        with unit variance :math:`v(\mu)`.

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.
        """
        return -2 * (y - mu) / self.unit_variance(mu)

    def deviance(self, y, mu, weights=1):
        r"""Compute the deviance.

        The deviance is a weighted sum of the per sample unit deviances,
        :math:`D = \sum_i s_i \cdot d(y_i, \mu_i)`
        with weights :math:`s_i` and unit deviance :math:`d(y,\mu)`.
        In terms of the log-likelihood it is :math:`D = -2\phi\cdot
        \left(loglike(y,\mu,\frac{phi}{s})
        - loglike(y,y,\frac{phi}{s})\right)`.

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.

        weights : array, shape (n_samples,) (default=1)
            Weights or exposure to which variance is inverse proportional.
        """
        return np.sum(weights * self.unit_deviance(y, mu))

    def deviance_derivative(self, y, mu, weights=1):
        """Compute the derivative of the deviance w.r.t. mu.

        It gives :math:`\\frac{\\partial}{\\partial\\mu} D(y, \\mu; weights)`.

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.

        weights : array, shape (n_samples,) (default=1)
            Weights or exposure to which variance is inverse proportional.
        """
        return weights * self.unit_deviance_derivative(y, mu)

    def starting_mu(self, y, weights=1, ind_weight=0.5):
        """Set starting values for the mean mu.

        These may be good starting points for the (unpenalized) IRLS solver.

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        weights : array, shape (n_samples,) (default=1)
            Weights or exposure to which variance is inverse proportional.

        ind_weight : float (default=0.5)
            Must be between 0 and 1. Specifies how much weight is given to the
            individual observations instead of the mean of y.
        """
        return (ind_weight * y +
                (1. - ind_weight) * np.average(y, weights=weights))

    def _mu_deviance_derivative(self, coef, X, y, weights, link):
        """Compute mu and the derivative of the deviance w.r.t coef."""
        lin_pred = _safe_lin_pred(X, coef)
        mu = link.inverse(lin_pred)
        d1 = link.inverse_derivative(lin_pred)
        temp = d1 * self.deviance_derivative(y, mu, weights)
        if coef.size == X.shape[1] + 1:
            devp = np.concatenate(([temp.sum()], temp @ X))
        else:
            devp = temp @ X  # sampe as X.T @ temp
        return mu, devp

    def _score(self, coef, phi, X, y, weights, link):
        r"""Compute the score function.

        The score function is the derivative of the
        log-likelihood w.r.t. `coef` (:math:`w`).
        It is given by

        .. math:

            \mathbf{score}(\boldsymbol{w})
            = \frac{\partial loglike}{\partial\boldsymbol{w}}
            = \mathbf{X}^T \mathbf{D}
            \boldsymbol{\Sigma}^-1 (\mathbf{y} - \boldsymbol{\mu})\,,

        with :math:`\mathbf{D}=\mathrm{diag}(h'(\eta_1),\ldots)` and
        :math:`\boldsymbol{\Sigma}=\mathrm{diag}(\mathbf{V}[y_1],\ldots)`.
        Note: The derivative of the deviance w.r.t. coef equals -2 * score.
        """
        lin_pred = _safe_lin_pred(X, coef)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weights=weights)
        d = link.inverse_derivative(lin_pred)
        temp = sigma_inv * d * (y - mu)
        if coef.size == X.shape[1] + 1:
            score = np.concatenate(([temp.sum()], temp @ X))
        else:
            score = temp @ X  # sampe as X.T @ temp
        return score

    def _fisher_matrix(self, coef, phi, X, y, weights, link):
        r"""Compute the Fisher information matrix.

        The Fisher information matrix, also known as expected information
        matrix is given by

        .. math:

            \mathbf{F}(\boldsymbol{w}) =
            \mathrm{E}\left[-\frac{\partial\mathbf{score}}{\partial
            \boldsymbol{w}} \right]
            = \mathrm{E}\left[
            -\frac{\partial^2 loglike}{\partial\boldsymbol{w}
            \partial\boldsymbol{w}^T}\right]
            = \mathbf{X}^T W \mathbf{X} \,,

        with :math:`\mathbf{W} = \mathbf{D}^2 \boldsymbol{\Sigma}^{-1}`,
        see func:`_score`.
        """
        lin_pred = _safe_lin_pred(X, coef)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weights=weights)
        d = link.inverse_derivative(lin_pred)
        d2_sigma_inv = sigma_inv * d * d
        intercept = (coef.size == X.shape[1] + 1)
        fisher_matrix = _safe_sandwich_dot(X, d2_sigma_inv,
                                           intercept=intercept)
        return fisher_matrix

    def _observed_information(self, coef, phi, X, y, weights, link):
        r"""Compute the observed information matrix.

        The observed information matrix, also known as the negative of
        the Hessian matrix of the log-likelihood, is given by

        .. math:

            \mathbf{H}(\boldsymbol{w}) =
            -\frac{\partial^2 loglike}{\partial\boldsymbol{w}
            \partial\boldsymbol{w}^T}
            = \mathbf{X}^T \left[
            - \mathbf{D}' \mathbf{R}
            + \mathbf{D}^2 \mathbf{V} \mathbf{R}
            + \mathbf{D}^2
            \right] \boldsymbol{\Sigma}^{-1} \mathbf{X} \,,

        with :math:`\mathbf{R} = \mathrm{diag}(y_i - \mu_i)`,
        :math:`\mathbf{V} = \mathrm{diag}\left(\frac{v'(\mu_i)}{
        v(\mu_i)}
        \right)`,
        see :func:`score_` function and :func:`_fisher_matrix`.
        """
        lin_pred = _safe_lin_pred(X, coef)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weights=weights)
        dp = link.inverse_derivative2(lin_pred)
        d2 = link.inverse_derivative(lin_pred)**2
        v = self.unit_variance_derivative(mu)/self.unit_variance(mu)
        r = y - mu
        temp = sigma_inv * (-dp * r + d2 * v * r + d2)
        intercept = (coef.size == X.shape[1] + 1)
        observed_information = _safe_sandwich_dot(X, temp,
                                                  intercept=intercept)
        return observed_information

    def _eta_mu_score_fisher(self, coef, phi, X, y, weights, link,
                             diag_fisher=False):
        """Compute linear predictor, mean, score function and fisher matrix.

        It calculates the linear predictor, the mean, score function
        (derivative of log-likelihood) and Fisher information matrix
        all in one go as function of `coef` (:math:`w`) and the data.

        Parameters
        ----------
        diag_fisher : boolean, optional (default=False)
            If ``True``, returns only an array d such that
            fisher = X.T @ np.diag(d) @ X.

        Returns
        -------
        (eta, mu, score, fisher) : tuple with 4 elements
            The 4 elements are:

            * eta: ndarray, shape (X.shape[0],)
            * mu: ndarray, shape (X.shape[0],)
            * score: ndarray, shape (X.shape[0],)
            * fisher:

                * If diag_fisher is ``False``, the full fisher matrix,
                  an array of shape (X.shape[1], X.shape[1])
                * If diag_fisher is ``True`, an array of shape (X.shape[0])
        """
        intercept = (coef.size == X.shape[1] + 1)
        # eta = linear predictor
        eta = _safe_lin_pred(X, coef)
        mu = link.inverse(eta)
        sigma_inv = 1./self.variance(mu, phi=phi, weights=weights)
        d1 = link.inverse_derivative(eta)  # = h'(eta)
        # Alternatively:
        # h'(eta) = h'(g(mu)) = 1/g'(mu), note that h is inverse of g
        # d1 = 1./link.derivative(mu)
        d1_sigma_inv = d1 * sigma_inv
        temp = d1_sigma_inv * (y - mu)
        if intercept:
            score = np.concatenate(([temp.sum()], temp @ X))
        else:
            score = temp @ X

        d2_sigma_inv = d1 * d1_sigma_inv
        if diag_fisher:
            fisher_matrix = d2_sigma_inv
        else:
            fisher_matrix = _safe_sandwich_dot(X, d2_sigma_inv,
                                               intercept=intercept)
        return eta, mu, score, fisher_matrix


class TweedieDistribution(ExponentialDispersionModel):
    r"""A class for the Tweedie distribution.

    A Tweedie distribution with mean :math:`\mu=\mathrm{E}[Y]` is uniquely
    defined by it's mean-variance relationship
    :math:`\mathrm{Var}[Y] \propto \mu^power`.

    Special cases are:

    ===== ================
    Power Distribution
    ===== ================
    0     Normal
    1     Poisson
    (0,1) Compound Poisson
    2     Gamma
    3     Inverse Gaussian

    Parameters
    ----------
    power : float (default=0)
            The variance power of the `unit_variance`
            :math:`v(\mu) = \mu^{power}`.
            For ``0<power<1``, no distribution exists.
    """
    def __init__(self, power=0):
        # validate power and set _upper_bound, _include_upper_bound attrs
        self.power = power

    @property
    def power(self):
        return self._power

    @power.setter
    def power(self, power):
        if not isinstance(power, numbers.Real):
            raise TypeError('power must be a real number, input was {0}'
                            .format(power))

        self._upper_bound = np.Inf
        self._include_upper_bound = False
        if power < 0:
            # Extreme Stable
            self._lower_bound = -np.Inf
            self._include_lower_bound = False
        elif power == 0:
            # NormalDistribution
            self._lower_bound = -np.Inf
            self._include_lower_bound = False
        elif (power > 0) and (power < 1):
            raise ValueError('For 0<power<1, no distribution exists.')
        elif power == 1:
            # PoissonDistribution
            self._lower_bound = 0
            self._include_lower_bound = True
        elif (power > 1) and (power < 2):
            # Compound Poisson
            self._lower_bound = 0
            self._include_lower_bound = True
        elif power == 2:
            # GammaDistribution
            self._lower_bound = 0
            self._include_lower_bound = False
        elif (power > 2) and (power < 3):
            # Positive Stable
            self._lower_bound = 0
            self._include_lower_bound = False
        elif power == 3:
            # InverseGaussianDistribution
            self._lower_bound = 0
            self._include_lower_bound = False
        elif power > 3:
            # Positive Stable
            self._lower_bound = 0
            self._include_lower_bound = False
        else:  # pragma: no cover
            # this branch should be unreachable.
            raise ValueError

        self._power = power

    def unit_variance(self, mu):
        """Compute the unit variance of a Tweedie distribution v(mu)=mu**power.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.
        """
        return np.power(mu, self.power)

    def unit_variance_derivative(self, mu):
        """Compute the derivative of the unit variance of a Tweedie
        distribution v(mu)=power*mu**(power-1).

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.
        """
        return self.power * np.power(mu, self.power - 1)

    def unit_deviance(self, y, mu):
        p = self.power
        if p == 0:
            # NormalDistribution
            return (y - mu)**2
        if p == 1:
            # PoissonDistribution
            # 2 * (y*log(y/mu) - y + mu), with y*log(y/mu)=0 if y=0
            return 2 * (special.xlogy(y, y/mu) - y + mu)
        elif p == 2:
            # GammaDistribution
            return 2 * (np.log(mu/y) + y/mu - 1)
        else:
            # return 2 * (np.maximum(y,0)**(2-p)/((1-p)*(2-p))
            #    - y*mu**(1-p)/(1-p) + mu**(2-p)/(2-p))
            return 2 * (np.power(np.maximum(y, 0), 2-p)/((1-p)*(2-p)) -
                        y*np.power(mu, 1-p)/(1-p) + np.power(mu, 2-p)/(2-p))


class NormalDistribution(TweedieDistribution):
    """Class for the Normal (aka Gaussian) distribution"""
    def __init__(self):
        super(NormalDistribution, self).__init__(power=0)


class PoissonDistribution(TweedieDistribution):
    """Class for the scaled Poisson distribution"""
    def __init__(self):
        super(PoissonDistribution, self).__init__(power=1)


class GammaDistribution(TweedieDistribution):
    """Class for the Gamma distribution"""
    def __init__(self):
        super(GammaDistribution, self).__init__(power=2)


class InverseGaussianDistribution(TweedieDistribution):
    """Class for the scaled InverseGaussianDistribution distribution"""
    def __init__(self):
        super(InverseGaussianDistribution, self).__init__(power=3)


class GeneralizedHyperbolicSecant(ExponentialDispersionModel):
    """A class for the Generalized Hyperbolic Secant (GHS) distribution.

    The GHS distribution is for targets y in (-inf, inf).
    """
    def __init__(self):
        self._lower_bound = -np.Inf
        self._upper_bound = np.Inf
        self._include_lower_bound = False
        self._include_upper_bound = False

    def unit_variance(self, mu):
        return 1 + mu**2

    def unit_variance_derivative(self, mu):
        return 2 * mu

    def unit_deviance(self, y, mu):
        return (2 * y * (np.arctan(y) - np.arctan(mu)) +
                np.log((1 + mu**2)/(1 + y**2)))


class BinomialDistribution(ExponentialDispersionModel):
    """A class for the Binomial distribution.

    The Binomial distribution is for targets y in [0, 1].
    """
    def __init__(self):
        self._lower_bound = 0
        self._upper_bound = 1
        self._include_lower_bound = True
        self._include_upper_bound = True

    def unit_variance(self, mu):
        return mu * (1 - mu)

    def unit_variance_derivative(self, mu):
        return 1 - 2 * mu

    def unit_deviance(self, y, mu):
        return 2 * (special.xlogy(y, y/mu) + special.xlogy(1-y, (1-y)/(1-mu)))


def _irls_step(X, W, P2, z, fit_intercept=True):
    """Compute one step in iteratively reweighted least squares.

    Solve A w = b for w with
    A = (X' W X + P2)
    b = X' W z
    z = eta + D^-1 (y-mu)

    See also fit method of :class:`GeneralizedLinearRegressor`.

    Parameters
    ----------
    X : {ndarray, sparse matrix}, shape (n_samples, n_features)
        Training data (with intercept included if present)

    W : ndarray, shape (n_samples,)

    P2 : {ndarray, sparse matrix}, shape (n_features, n_features)
        The L2-penalty matrix or vector (=diagonal matrix)

    z : ndarray, shape (n_samples,)
        Working observations

    fit_intercept : boolean, optional (default=True)

    Returns
    -------
    coef : ndarray, shape (c,)
        If fit_intercept=False, shape c=X.shape[1].
        If fit_intercept=True, then c=X.shapee[1] + 1.
    """
    # Note: solve vs least squares, what is more appropriate?
    #       scipy.linalg.solve seems faster, but scipy.linalg.lstsq
    #       is more robust.
    # Note: X.T @ W @ X is not sparse, even when X is sparse.
    #      Sparse solver would splinalg.spsolve(A, b) or splinalg.lsmr(A, b)
    if fit_intercept:
        Wz = W * z
        if sparse.issparse(X):
            b = np.concatenate(([Wz.sum()], X.transpose() @ Wz))
        else:
            b = np.concatenate(([Wz.sum()], X.T @ Wz))
        A = _safe_sandwich_dot(X, W, intercept=fit_intercept)
        if P2.ndim == 1:
            idx = np.arange(start=1, stop=A.shape[0])
            A[(idx, idx)] += P2  # add to diag elements without intercept
        elif sparse.issparse(P2):
            A[1:, 1:] += P2.toarray()
        else:
            A[1:, 1:] += P2
    else:
        if sparse.issparse(X):
            XtW = X.transpose().multiply(W)
            # for older versions of numpy and scipy, A may be a np.matrix
            A = _safe_toarray(XtW @ X)
        else:
            XtW = (X.T * W)
            A = XtW @ X
        b = XtW @ z
        if P2.ndim == 1:
            A[np.diag_indices_from(A)] += P2
        elif sparse.issparse(P2):
            A += P2.toarray()
        else:
            A += P2

    coef, *_ = linalg.lstsq(A, b, overwrite_a=True, overwrite_b=True)
    return coef


def _irls_solver(coef, X, y, weights, P2, fit_intercept, family, link,
                 max_iter, tol):
    """Solve GLM with L2 penalty by IRLS algorithm.

    Note: If X is sparse, P2 must also be sparse.
    """
    # Solve Newton-Raphson (1): Obj'' (w - w_old) = -Obj'
    #   Obj = objective function = 1/2 Dev + l2/2 w P2 w
    #   Dev = deviance, s = normalized weights, variance V(mu) but phi=1
    #   D   = link.inverse_derivative(eta) = diag_matrix(h'(X w))
    #   D2  = link.inverse_derivative(eta)^2 = D^2
    #   W   = D2/V(mu)
    #   l2  = alpha * (1 - l1_ratio)
    #   Obj' = d(Obj)/d(w) = 1/2 Dev' + l2 P2 w
    #        = -X' D (y-mu)/V(mu) + l2 P2 w
    #   Obj''= d2(Obj)/d(w)d(w') = Hessian = -X'(...) X + l2 P2
    #   Use Fisher matrix instead of full info matrix -X'(...) X,
    #    i.e. E[Dev''] with E[y-mu]=0:
    #   Obj'' ~ X' W X + l2 P2
    # (1): w = (X' W X + l2 P2)^-1 X' W z,
    #      with z = eta + D^-1 (y-mu)
    # Note: P2 must be symmetrized
    # Note: ' denotes derivative, but also transpose for matrices

    eta = _safe_lin_pred(X, coef)
    mu = link.inverse(eta)
    # D = h'(eta)
    hp = link.inverse_derivative(eta)
    V = family.variance(mu, phi=1, weights=weights)

    converged = False
    n_iter = 0
    while n_iter < max_iter:
        n_iter += 1
        # coef_old not used so far.
        # coef_old = coef
        # working weights W, in principle a diagonal matrix
        # therefore here just as 1d array
        W = hp**2 / V
        # working observations
        z = eta + (y - mu) / hp
        # solve A*coef = b
        # A = X' W X + P2, b = X' W z
        coef = _irls_step(X, W, P2, z, fit_intercept=fit_intercept)
        # updated linear predictor
        # do it here for updated values for tolerance
        eta = _safe_lin_pred(X, coef)
        mu = link.inverse(eta)
        hp = link.inverse_derivative(eta)
        V = family.variance(mu, phi=1, weights=weights)

        # which tolerace? |coef - coef_old| or gradient?
        # use gradient for compliance with newton-cg and lbfgs
        # gradient = -X' D (y-mu)/V(mu) + l2 P2 w
        temp = hp * (y - mu) / V
        if sparse.issparse(X):
            gradient = -(X.transpose() @ temp)
        else:
            gradient = -(X.T @ temp)
        idx = 1 if fit_intercept else 0  # offset if coef[0] is intercept
        if P2.ndim == 1:
            gradient += P2 * coef[idx:]
        else:
            gradient += P2 @ coef[idx:]
        if fit_intercept:
            gradient = np.concatenate(([-temp.sum()], gradient))
        if (np.max(np.abs(gradient)) <= tol):
            converged = True
            break

    if not converged:
        warnings.warn("irls failed to converge. Increase the number "
                      "of iterations (currently {0})"
                      .format(max_iter), ConvergenceWarning)

    return coef, n_iter


def _cd_cycle(d, X, coef, score, fisher, P1, P2, n_cycles, inner_tol,
              max_inner_iter=1000, selection='cyclic',
              random_state=None, diag_fisher=False):
    """Compute inner loop of coordinate descent, i.e. cycles through features.

    Minimization of 1-d subproblems::

        min_z q(d+z*e_j) - q(d)
        = min_z A_j z + 1/2 B_jj z^2 + ||P1_j (w_j+d_j+z)||_1

    A = f'(w) + d*H(w) + (w+d)*P2
    B = H+P2
    Note: f'=-score and H=fisher are updated at the end of outer iteration.
    """
    # TODO: use sparsity (coefficient already 0 due to L1 penalty)
    #       => active set of features for featurelist, see paper
    #          of Improved GLMNET or Gap Safe Screening Rules
    #          https://arxiv.org/abs/1611.05780
    n_samples, n_features = X.shape
    intercept = (coef.size == X.shape[1] + 1)
    idx = 1 if intercept else 0  # offset if coef[0] is intercept
    B = fisher
    if P2.ndim == 1:
        coef_P2 = coef[idx:] * P2
        if not diag_fisher:
            idiag = np.arange(start=idx, stop=B.shape[0])
            # B[np.diag_indices_from(B)] += P2
            B[(idiag, idiag)] += P2
    else:
        coef_P2 = coef[idx:] @ P2
        if not diag_fisher:
            if sparse.issparse(P2):
                B[idx:, idx:] += P2.toarray()
            else:
                B[idx:, idx:] += P2
    A = -score
    A[idx:] += coef_P2
    # A += d @ (H+P2) but so far d=0
    # inner loop
    for inner_iter in range(1, max_inner_iter+1):
        inner_iter += 1
        n_cycles += 1
        # cycle through features, update intercept separately at the end
        if selection == 'random':
            featurelist = random_state.permutation(n_features)
        else:
            featurelist = np.arange(n_features)
        for j in featurelist:
            # minimize_z: a z + 1/2 b z^2 + c |d+z|
            # a = A_j
            # b = B_jj > 0
            # c = |P1_j| = P1_j > 0, see 1.3
            # d = w_j + d_j
            # cf. https://arxiv.org/abs/0708.1485 Eqs. (3) - (4)
            # with beta = z+d, beta_hat = d-a/b and gamma = c/b
            # z = 1/b * S(bd-a,c) - d
            # S(a,b) = sign(a) max(|a|-b, 0) soft thresholding
            jdx = j+idx  # index for arrays containing entries for intercept
            a = A[jdx]
            if diag_fisher:
                # Note: fisher is ndarray of shape (n_samples,) => no idx
                # Calculate Bj = B[j, :] = B[:, j] as it is needed later anyway
                Bj = np.zeros_like(A)
                if intercept:
                    Bj[0] = fisher.sum()
                if sparse.issparse(X):
                    Bj[idx:] = _safe_toarray(X[:, j].transpose() @
                                             X.multiply(fisher[:, np.newaxis])
                                             ).ravel()
                else:
                    Bj[idx:] = (fisher * X[:, j]) @ X

                if P2.ndim == 1:
                    Bj[idx:] += P2[j]
                else:
                    if sparse.issparse(P2):
                        # slice columns as P2 is csc
                        Bj[idx:] += P2[:, j].toarray().ravel()
                    else:
                        Bj[idx:] += P2[:, j]
                b = Bj[jdx]
            else:
                b = B[jdx, jdx]

            # those ten lines are what it is all about
            if b <= 0:
                z = 0
            elif P1[j] == 0:
                z = -a/b
            elif a + P1[j] < b * (coef[jdx] + d[jdx]):
                z = -(a + P1[j])/b
            elif a - P1[j] > b * (coef[jdx] + d[jdx]):
                z = -(a - P1[j])/b
            else:
                z = -(coef[jdx] + d[jdx])

            # update direction d
            d[jdx] += z
            # update A because d_j is now d_j+z
            # A = f'(w) + d*H(w) + (w+d)*P2
            # => A += (H+P2)*e_j z = B_j * z
            # Note: B is symmetric B = B.transpose
            if diag_fisher:
                # Bj = B[:, j] calculated above, still valid
                A += Bj * z
            else:
                # B is symmetric, C- or F-contiguous, but never sparse
                if B.flags['F_CONTIGUOUS']:
                    # slice columns like for sparse csc
                    A += B[:, jdx] * z
                else:  # B.flags['C_CONTIGUOUS'] might be true
                    # slice rows
                    A += B[jdx, :] * z
            # end of cycle over features
        # update intercept
        if intercept:
            if diag_fisher:
                Bj = np.zeros_like(A)
                Bj[0] = fisher.sum()
                Bj[1:] = fisher @ X
                b = Bj[0]
            else:
                b = B[0, 0]
            z = 0 if b <= 0 else -A[0]/b
            d[0] += z
            if diag_fisher:
                A += Bj * z
            else:
                if B.flags['F_CONTIGUOUS']:
                    A += B[:, 0] * z
                else:
                    A += B[0, :] * z
        # end of complete cycle
        # stopping criterion for inner loop
        # sum_i(|minimum of norm of subgrad of q(d)_i|)
        # subgrad q(d) = A + subgrad ||P1*(w+d)||_1
        mn_subgrad = _min_norm_sugrad(coef=coef + d, grad=A, P2=None, P1=P1)
        mn_subgrad = linalg.norm(mn_subgrad, ord=1)
        if mn_subgrad <= inner_tol:
            if inner_iter == 1:
                inner_tol = inner_tol/4.
            break
        # end of inner loop
    return d, coef_P2, n_cycles, inner_tol


def _cd_solver(coef, X, y, weights, P1, P2, fit_intercept, family, link,
               max_iter=100, max_inner_iter=1000, tol=1e-4,
               selection='cyclic ', random_state=None,
               diag_fisher=False, copy_X=True):
    """Solve GLM with L1 and L2 penalty by coordinate descent algorithm.

    The objective being minimized in the coefficients w=coef is::

        F = f + g, f(w) = 1/2 deviance, g = 1/2 w*P2*w + ||P1*w||_1

    An Improved GLMNET for L1-regularized Logistic Regression:

    1. Find optimal descent direction d by minimizing
       min_d F(w+d) = min_d F(w+d) - F(w)
    2. Quadratic approximation of F(w+d)-F(w) = q(d):
       using f(w+d) = f(w) + f'(w)*d + 1/2 d*H(w)*d + O(d^3) gives:
       q(d) = (f'(w) + w*P2)*d + 1/2 d*(H(w)+P2)*d
       + ||P1*(w+d)||_1 - ||P1*w||_1
       Then minimize q(d): min_d q(d)
    3. Coordinate descent by updating coordinate j (d -> d+z*e_j):
       min_z q(d+z*e_j)
       = min_z q(d+z*e_j) - q(d)
       = min_z A_j z + 1/2 B_jj z^2
               + ||P1_j (w_j+d_j+z)||_1 - ||P1_j (w_j+d_j)||_1
       A = f'(w) + d*H(w) + (w+d)*P2
       B = H + P2

    Repeat steps 1-3 until convergence.
    Note: Use Fisher matrix instead of Hessian for H.
    Note: f' = -score, H = Fisher matrix

    Parameters
    ----------
    coef : ndarray, shape (c,)
        If fit_intercept=False, shape c=X.shape[1].
        If fit_intercept=True, then c=X.shape[1] + 1.

    X : {ndarray, csc sparse matrix}, shape (n_samples, n_features)
        Training data (with intercept included if present). If not sparse,
        pass directly as Fortran-contiguous data to avoid
        unnecessary memory duplication.

    y : ndarray, shape (n_samples,)
        Target values.

    weights: ndarray, shape (n_samples,)
        Sample weights with which the deviance is weighted. The weights must
        bee normalized and sum to 1.

    P1 : {ndarray}, shape (n_features,)
        The L1-penalty vector (=diagonal matrix)

    P2 : {ndarray, csc sparse matrix}, shape (n_features, n_features)
        The L2-penalty matrix or vector (=diagonal matrix). If a matrix is
        passed, it must be symmetric. If X is sparse, P2 must also be sparse.

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    family : ExponentialDispersionModel

    link : Link

    max_iter : int, optional (default=100)
        Maximum numer of outer (Newton) iterations.

    max_inner_iter : int, optional (default=1000)
        Maximum number of iterations in each inner loop, i.e. max number of
        cycles over all features per inner loop.

    tol : float, optional (default=1e-4)
        Convergence criterion is
        sum_i(|minimum of norm of subgrad of objective_i|)<=tol.

    selection : str, optional (default='cyclic')
        If 'random', randomly chose features in inner loop.

    random_state : {int, RandomState instance, None}, optional (default=None)

    diag_fisher : boolean, optional (default=False)
        ``False`` calculates full fisher matrix, ``True`` only diagonal matrix
        s.t. fisher = X.T @ diag @ X. This saves storage but needs more
        matrix-vector multiplications.

    copy_X : boolean, optional (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    Returns
    -------
    coef : ndarray, shape (c,)
        If fit_intercept=False, shape c=X.shape[1].
        If fit_intercept=True, then c=X.shape[1] + 1.

    n_iter : number of outer iterations = newton iterations

    n_cycles : number of cycles over features

    References
    ----------
    Guo-Xun Yuan, Chia-Hua Ho, Chih-Jen Lin
    An Improved GLMNET for L1-regularized Logistic Regression,
    Journal of Machine Learning Research 13 (2012) 1999-2030
    https://www.csie.ntu.edu.tw/~cjlin/papers/l1_glmnet/long-glmnet.pdf
    """
    X = check_array(X, 'csc', dtype=[np.float64, np.float32],
                    order='F', copy=copy_X)
    if P2.ndim == 2:
        P2 = check_array(P2, 'csc', dtype=[np.float64, np.float32],
                         order='F', copy=copy_X)
    if sparse.issparse(X):
        if not sparse.isspmatrix_csc(P2):
            raise ValueError("If X is sparse, P2 must also be sparse csc"
                             "format. Got P2 not sparse.")
    random_state = check_random_state(random_state)
    # Note: we already set P2 = l2*P2, P1 = l1*P1
    # Note: we already symmetrized P2 = 1/2 (P2 + P2')
    n_iter = 0  # number of outer iterations
    n_cycles = 0  # number of (complete) cycles over features
    converged = False
    n_samples, n_features = X.shape
    idx = 1 if fit_intercept else 0  # offset if coef[0] is intercept
    # line search parameters
    (beta, sigma) = (0.5, 0.01)
    # some precalculations
    # Note: For diag_fisher=False, fisher = X.T @ fisher @ X and fisher is a
    #       1d array representing a diagonal matrix.
    eta, mu, score, fisher = family._eta_mu_score_fisher(
        coef=coef, phi=1, X=X, y=y, weights=weights, link=link,
        diag_fisher=diag_fisher)
    # set up space for search direction d for inner loop
    d = np.zeros_like(coef)
    # initial stopping tolerance of inner loop
    # use L1-norm of minimum of norm of subgradient of F
    inner_tol = _min_norm_sugrad(coef=coef, grad=-score, P2=P2, P1=P1)
    inner_tol = linalg.norm(inner_tol, ord=1)
    # outer loop
    while n_iter < max_iter:
        n_iter += 1
        # initialize search direction d (to be optimized) with zero
        d.fill(0)
        # inner loop = _cd_cycle
        d, coef_P2, n_cycles, inner_tol = \
            _cd_cycle(d, X, coef, score, fisher, P1, P2, n_cycles, inner_tol,
                      max_inner_iter=max_inner_iter, selection=selection,
                      random_state=random_state, diag_fisher=diag_fisher)
        # line search by sequence beta^k, k=0, 1, ..
        # F(w + lambda d) - F(w) <= lambda * bound
        # bound = sigma * (f'(w)*d + w*P2*d
        #                  +||P1 (w+d)||_1 - ||P1 w||_1)
        P1w_1 = linalg.norm(P1 * coef[idx:], ord=1)
        P1wd_1 = linalg.norm(P1 * (coef + d)[idx:], ord=1)
        # Note: coef_P2 already calculated and still valid
        bound = sigma * (-(score @ d) + coef_P2 @ d[idx:] + P1wd_1 - P1w_1)
        Fw = (0.5 * family.deviance(y, mu, weights) +
              0.5 * (coef_P2 @ coef[idx:]) + P1w_1)
        la = 1./beta
        for k in range(20):
            la *= beta  # starts with la=1
            coef_wd = coef + la * d
            mu_wd = link.inverse(_safe_lin_pred(X, coef_wd))
            Fwd = (0.5 * family.deviance(y, mu_wd, weights) +
                   linalg.norm(P1 * coef_wd[idx:], ord=1))
            if P2.ndim == 1:
                Fwd += 0.5 * ((coef_wd[idx:] * P2) @ coef_wd[idx:])
            else:
                Fwd += 0.5 * (coef_wd[idx:] @ (P2 @ coef_wd[idx:]))
            if Fwd - Fw <= sigma * la * bound:
                break
        # update coefficients
        coef += la * d
        # calculate eta, mu, score, Fisher matrix for next iteration
        eta, mu, score, fisher = family._eta_mu_score_fisher(
            coef=coef, phi=1, X=X, y=y, weights=weights, link=link,
            diag_fisher=diag_fisher)
        # stopping criterion for outer loop
        # sum_i(|minimum-norm of subgrad of F(w)_i|)
        # fp_wP2 = f'(w) + w*P2
        # Note: eta, mu and score are already updated
        mn_subgrad = _min_norm_sugrad(coef=coef, grad=-score, P2=P2, P1=P1)
        mn_subgrad = linalg.norm(mn_subgrad, ord=1)
        if mn_subgrad <= tol:
            converged = True
            break
        # end of outer loop
    if not converged:
        warnings.warn("Coordinate descent failed to converge. Increase"
                      " the maximum number of iterations max_iter"
                      " (currently {0})".format(max_iter), ConvergenceWarning)

    return coef, n_iter, n_cycles


class GeneralizedLinearRegressor(BaseEstimator, RegressorMixin):
    """Regression via a Generalized Linear Model (GLM) with penalties.

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w). Therefore,
    the fit minimizes the following objective function with combined L1 and L2
    priors as regularizer::

            1/(2*sum(s)) * deviance(y, h(X*w); s)
            + alpha * l1_ratio * ||P1*w||_1
            + 1/2 * alpha * (1 - l1_ratio) * w*P2*w

    with inverse link function h and s=sample_weight. Note that for
    ``sample_weight=None``, one has s_i=1 and sum(s)=n_samples).
    For ``P1=P2='identity'``, the penalty is the elastic net::

            alpha * l1_ratio * ||w||_1
            + 1/2 * alpha * (1 - l1_ratio) * ||w||_2^2

    If you are interested in controlling the L1 and L2 penalties
    separately, keep in mind that this is equivalent to::

            a * L1 + b * L2

    where::

            alpha = a + b and l1_ratio = a / (a + b)

    The parameter ``l1_ratio`` corresponds to alpha in the R package glmnet,
    while ``alpha`` corresponds to the lambda parameter in glmnet.
    Specifically, l1_ratio = 1 is the lasso penalty.

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    l1_ratio : float, optional (default=0)
        The elastic net mixing parameter, with ``0 <= l1_ratio <= 1``. For
        ``l1_ratio = 0`` the penalty is an L2 penalty. ``For l1_ratio = 1`` it
        is an L1 penalty.  For ``0 < l1_ratio < 1``, the penalty is a
        combination of L1 and L2.

    P1 : {'identity', array-like}, shape (n_features,), optional \
            (default='identity')
        With this array, you can exclude coefficients from the L1 penalty.
        Set the corresponding value to 1 (include) or 0 (exclude). The
        default value ``'identity'`` is the same as a 1d array of ones.
        Note that n_features = X.shape[1].

    P2 : {'identity', array-like, sparse matrix}, shape \
            (n_features,) or (n_features, n_features), optional \
            (default='identity')
        With this option, you can set the P2 matrix in the L2 penalty `w*P2*w`.
        This gives a fine control over this penalty (Tikhonov regularization).
        A 2d array is directly used as the square matrix P2. A 1d array is
        interpreted as diagonal (square) matrix. The default 'identity' sets
        the identity matrix, which gives the usual squared L2-norm. If you just
        want to exclude certain coefficients, pass a 1d array filled with 1,
        and 0 for the coefficients to be excluded.
        Note that P2 must be positive semi-definite.

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    family : {'normal', 'poisson', 'gamma', 'inverse.gaussian', 'binomial'} \
            or an instance of class ExponentialDispersionModel, \
            optional(default='normal')
        The distributional assumption of the GLM, i.e. which distribution from
        the EDM, specifies the loss function to be minimized.

    link : {'auto', 'identity', 'log', 'logit'} or an instance of class Link, \
            optional (default='auto')
        The link function of the GLM, i.e. mapping from linear predictor
        (X*coef) to expectation (mu). Option 'auto' sets the link depending on
        the chosen family as follows:

        - 'identity' for family 'normal'

        - 'log' for families 'poisson', 'gamma', 'inverse.gaussian'

        - 'logit' for family 'binomial'

    fit_dispersion : {None, 'chisqr', 'deviance'}, optional (default=None)
        Method for estimation of the dispersion parameter phi. Whether to use
        the chi squared statistic or the deviance statistic. If None, the
        dispersion is not estimated.

    solver : {'auto', 'cd', 'irls', 'lbfgs', 'newton-cg'}, \
            optional (default='auto')
        Algorithm to use in the optimization problem:

        'auto'
            Sets 'irls' if l1_ratio equals 0, else 'cd'.

        'cd'
            Coordinate descent algorithm. It can deal with L1 as well as L2
            penalties. Note that in order to avoid unnecessary memory
            duplication of X in the ``fit`` method, X should be directly passed
            as a Fortran-contiguous numpy array or sparse csc matrix.

        'irls'
            Iterated reweighted least squares.
            It is the standard algorithm for GLMs. It cannot deal with
            L1 penalties.

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer. It cannot deal with L1 penalties.

        'newton-cg', 'lbfgs'
            Newton conjugate gradient algorithm cannot deal with L1 penalties.

        Note that all solvers except lbfgs use the fisher matrix, i.e. the
        expected Hessian instead of the Hessian matrix.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the irls, newton-cg and lbfgs solvers,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function. For the cd solver, convergence is reached
        when ``sum_i(|minimum-norm of g_i|)``, where ``g_i`` is the
        subgradient of the objective and minimum-norm of ``g_i`` is the element
        of the subgradient ``g_i`` with the smallest L2-norm.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_`` (supersedes option
        ``start_params``). If set to ``True`` or if the attribute ``coef_``
        does not exit (first call to ``fit``), option ``start_params`` sets the
        start values for ``coef_`` and ``intercept_``.

    start_params : {'guess', 'zero', array of shape (n_features*, )}, \
            optional (default='guess')
        Relevant only if ``warm_start=False`` or if fit is called
        the first time (``self.coef_`` does not yet exist).

        'guess'
            Start values of mu are calculated by family.starting_mu(..). Then,
            one Newton step obtains start values for ``coef_``. If
            ``solver='irls'``, it uses one irls step, else the Newton step is
            calculated by the cd solver.
            This gives usually good starting values.

        'zero'
        All coefficients are set to zero. If ``fit_intercept=True``, the
        start value for the intercept is obtained by the weighted average of y.

        array
        The array of size n_features* is directly used as start values
        for ``coef_``. If ``fit_intercept=True``, the first element
        is assumed to be the start value for the ``intercept_``.
        Note that n_features* = X.shape[1] + fit_intercept, i.e. it includes
        the intercept in counting.

    selection : str, optional (default='cyclic')
        For the solver 'cd' (coordinate descent), the coordinates (features)
        can be updated in either cyclic or random order.
        If set to 'random', a random coefficient is updated every iteration
        rather than looping over features sequentially in the same order. This
        (setting to 'random') often leads to significantly faster convergence
        especially when tol is higher than 1e-4.

    random_state : {int, RandomState instance, None}, optional (default=None)
        The seed of the pseudo random number generator that selects a random
        feature to be updated for solver 'cd' (coordinate descent).
        If int, random_state is the seed used by the random
        number generator; if RandomState instance, random_state is the random
        number generator; if None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    diag_fisher : boolean, optional, (default=False)
        Only relevant for solver 'cd' (see also ``start_params='guess'``).
        If ``False``, the full Fisher matrix (expected Hessian) is computed in
        each outer iteration (Newton iteration). If ``True``, only a diagonal
        matrix (stored as 1d array) is computed, such that
        fisher = X.T @ diag @ X. This saves memory and matrix-matrix
        multiplications, but needs more matrix-vector multiplications. If you
        use large sparse X or if you have many features,
        i.e. n_features >> n_samples, you might set this option to ``True``.

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    check_input : boolean, optional (default=True)
        Allow to bypass several checks on input: y values in range of family,
        sample_weight non-negative, P2 positive semi-definite.
        Don't use this parameter unless you know what you do.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    dispersion_ : float
        The dispersion parameter :math:`\\phi` if ``fit_dispersion`` was set.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.

    For alpha > 0, the feature matrix X should be standardized in order to
    penalize features equally strong. Call
    :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``.

    If the target y is a ratio, appropriate sample weights s should be
    provided.
    As an example, consider Poisson distributed counts z (integers) and
    weights s=exposure (time, money, persons years, ...). Then you fit
    y = z/s, i.e. ``GeneralizedLinearModel(family='poisson').fit(X, y,
    sample_weight=s)``. The weights are necessary for the right (finite
    sample) mean.
    Consider :math:`\\bar{y} = \\frac{\\sum_i s_i y_i}{\\sum_i s_i}`,
    in this case one might say that y has a 'scaled' Poisson distributions.
    The same holds for other distributions.

    References
    ----------
    For the coordinate descent implementation:
        * Guo-Xun Yuan, Chia-Hua Ho, Chih-Jen Lin
          An Improved GLMNET for L1-regularized Logistic Regression,
          Journal of Machine Learning Research 13 (2012) 1999-2030
          https://www.csie.ntu.edu.tw/~cjlin/papers/l1_glmnet/long-glmnet.pdf
    """
    def __init__(self, alpha=1.0, l1_ratio=0, P1='identity', P2='identity',
                 fit_intercept=True, family='normal', link='auto',
                 fit_dispersion=None, solver='auto', max_iter=100,
                 tol=1e-4, warm_start=False, start_params='guess',
                 selection='cyclic', random_state=None, diag_fisher=False,
                 copy_X=True, check_input=True, verbose=0):
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.P1 = P1
        self.P2 = P2
        self.fit_intercept = fit_intercept
        self.family = family
        self.link = link
        self.fit_dispersion = fit_dispersion
        self.solver = solver
        self.max_iter = max_iter
        self.tol = tol
        self.warm_start = warm_start
        self.start_params = start_params
        self.selection = selection
        self.random_state = random_state
        self.diag_fisher = diag_fisher
        self.copy_X = copy_X
        self.check_input = check_input
        self.verbose = verbose

    def fit(self, X, y, sample_weight=None):
        """Fit a Generalized Linear Model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : {None, array-like}, shape (n_samples,),\
                optional (default=None)
            Individual weights w_i for each sample. Note that for an
            Exponential Dispersion Model (EDM), one has
            Var[Y_i]=phi/w_i * v(mu).
            If Y_i ~ EDM(mu, phi/w_i), then
            sum(w*Y)/sum(w) ~ EDM(mu, phi/sum(w)), i.e. the mean of y is a
            weighted average with weights=sample_weight.

        Returns
        -------
        self : returns an instance of self.
        """
        #######################################################################
        # 1. input validation                                                 #
        #######################################################################
        # 1.1 validate arguments of __init__ ##################################
        # Guarantee that self._family_instance is an instance of class
        # ExponentialDispersionModel
        if isinstance(self.family, ExponentialDispersionModel):
            self._family_instance = self.family
        else:
            if self.family == 'normal':
                self._family_instance = NormalDistribution()
            elif self.family == 'poisson':
                self._family_instance = PoissonDistribution()
            elif self.family == 'gamma':
                self._family_instance = GammaDistribution()
            elif self.family == 'inverse.gaussian':
                self._family_instance = InverseGaussianDistribution()
            elif self.family == 'binomial':
                self._family_instance = BinomialDistribution()
            else:
                raise ValueError(
                    "The family must be an instance of class"
                    " ExponentialDispersionModel or an element of"
                    " ['normal', 'poisson', 'gamma', 'inverse.gaussian', "
                    "'binomial']; got (family={0})".format(self.family))

        # Guarantee that self._link_instance is set to an instance of
        # class Link
        if isinstance(self.link, Link):
            self._link_instance = self.link
        else:
            if self.link == 'auto':
                if isinstance(self._family_instance, TweedieDistribution):
                    if self._family_instance.power <= 0:
                        self._link_instance = IdentityLink()
                    if self._family_instance.power >= 1:
                        self._link_instance = LogLink()
                elif isinstance(self._family_instance,
                                GeneralizedHyperbolicSecant):
                    self._link_instance = IdentityLink()
                elif isinstance(self._family_instance, BinomialDistribution):
                    self._link_instance = LogitLink()
                else:
                    raise ValueError("No default link known for the "
                                     "specified distribution family. Please "
                                     "set link manually, i.e. not to 'auto'; "
                                     "got (link='auto', family={}"
                                     .format(self.family))
            elif self.link == 'identity':
                self._link_instance = IdentityLink()
            elif self.link == 'log':
                self._link_instance = LogLink()
            elif self.link == 'logit':
                self._link_instance = LogitLink()
            else:
                raise ValueError(
                    "The link must be an instance of class Link or "
                    "an element of ['auto', 'identity', 'log', 'logit']; "
                    "got (link={0})".format(self.link))

        if not isinstance(self.alpha, numbers.Number) or self.alpha < 0:
            raise ValueError("Penalty term must be a non-negative number;"
                             " got (alpha={0})".format(self.alpha))
        if (not isinstance(self.l1_ratio, numbers.Number) or
                self.l1_ratio < 0 or self.l1_ratio > 1):
            raise ValueError("l1_ratio must be a number in interval [0, 1];"
                             " got (l1_ratio={0})".format(self.l1_ratio))
        if not isinstance(self.fit_intercept, bool):
            raise ValueError("The argument fit_intercept must be bool;"
                             " got {0}".format(self.fit_intercept))
        if self.solver not in ['auto', 'irls', 'lbfgs', 'newton-cg', 'cd']:
            raise ValueError("GeneralizedLinearRegressor supports only solvers"
                             " 'auto', 'irls', 'lbfgs', 'newton-cg' and 'cd';"
                             " got {0}".format(self.solver))
        solver = self.solver
        if self.solver == 'auto':
            if self.l1_ratio == 0:
                solver = 'irls'
            else:
                solver = 'cd'
        if (self.alpha > 0 and self.l1_ratio > 0 and solver not in ['cd']):
            raise ValueError("The chosen solver (solver={0}) can't deal "
                             "with L1 penalties, which are included with "
                             "(alpha={1}) and (l1_ratio={2})."
                             .format(solver, self.alpha, self.l1_ratio))
        if (not isinstance(self.max_iter, int)
                or self.max_iter <= 0):
            raise ValueError("Maximum number of iteration must be a positive "
                             "integer;"
                             " got (max_iter={0!r})".format(self.max_iter))
        if not isinstance(self.tol, numbers.Number) or self.tol <= 0:
            raise ValueError("Tolerance for stopping criteria must be "
                             "positive; got (tol={0!r})".format(self.tol))
        if not isinstance(self.warm_start, bool):
            raise ValueError("The argument warm_start must be bool;"
                             " got {0}".format(self.warm_start))
        if self.selection not in ['cyclic', 'random']:
            raise ValueError("The argument selection must be 'cyclic' or "
                             "'random'; got (selection={0})"
                             .format(self.selection))
        random_state = check_random_state(self.random_state)
        if not isinstance(self.diag_fisher, bool):
            raise ValueError("The argument diag_fisher must be bool;"
                             " got {0}".format(self.diag_fisher))
        if not isinstance(self.copy_X, bool):
            raise ValueError("The argument copy_X must be bool;"
                             " got {0}".format(self.copy_X))
        if not isinstance(self.check_input, bool):
            raise ValueError("The argument check_input must be bool; got "
                             "(check_input={0})".format(self.check_input))

        family = self._family_instance
        link = self._link_instance

        # 1.2 validate arguments of fit #######################################
        _dtype = [np.float64, np.float32]
        if solver == 'cd':
            _stype = ['csc']
        else:
            _stype = ['csc', 'csr']
        X, y = check_X_y(X, y, accept_sparse=_stype,
                         dtype=_dtype, y_numeric=True, multi_output=False,
                         copy=self.copy_X)
        # Without converting y to float, deviance might raise
        # ValueError: Integers to negative integer powers are not allowed.
        # Also, y must not be sparse.
        y = np.asarray(y, dtype=np.float64)

        weights = _check_weights(sample_weight, y.shape[0])

        n_samples, n_features = X.shape

        # 1.3 arguments to take special care ##################################
        # P1, P2, start_params
        if isinstance(self.P1, str) and self.P1 == 'identity':
            P1 = np.ones(n_features)
        else:
            P1 = np.atleast_1d(self.P1)
            try:
                P1 = P1.astype(np.float64, casting='safe', copy=False)
            except TypeError:
                raise TypeError("The given P1 cannot be converted to a numeric"
                                "array; got (P1.dtype={0})."
                                .format(P1.dtype))
            if (P1.ndim != 1) or (P1.shape[0] != n_features):
                raise ValueError("P1 must be either 'identity' or a 1d array "
                                 "with the length of X.shape[1]; "
                                 "got (P1.shape[0]={0}), "
                                 "needed (X.shape[1]={1})."
                                 .format(P1.shape[0], n_features))
        # If X is sparse, make P2 sparse, too.
        if isinstance(self.P2, str) and self.P2 == 'identity':
            if sparse.issparse(X):
                P2 = (sparse.dia_matrix((np.ones(n_features), 0),
                      shape=(n_features, n_features))).tocsc()
            else:
                P2 = np.ones(n_features)
        else:
            P2 = check_array(self.P2, copy=True,
                             accept_sparse=_stype,
                             dtype=_dtype, ensure_2d=False)
            if P2.ndim == 1:
                P2 = np.asarray(P2)
                if P2.shape[0] != n_features:
                    raise ValueError("P2 should be a 1d array of shape "
                                     "(n_features,) with "
                                     "n_features=X.shape[1]; "
                                     "got (P2.shape=({0},)), needed ({1},)"
                                     .format(P2.shape[0], X.shape[1]))
                if sparse.issparse(X):
                    P2 = (sparse.dia_matrix((P2, 0),
                          shape=(n_features, n_features))).tocsc()
            elif (P2.ndim == 2 and P2.shape[0] == P2.shape[1] and
                    P2.shape[0] == X.shape[1]):
                if sparse.issparse(X):
                    P2 = sparse.csc_matrix(P2)
            else:
                raise ValueError("P2 must be either None or an array of shape "
                                 "(n_features, n_features) with "
                                 "n_features=X.shape[1]; "
                                 "got (P2.shape=({0}, {1})), needed ({2}, {2})"
                                 .format(P2.shape[0], P2.shape[1], X.shape[1]))

        start_params = self.start_params
        if isinstance(start_params, str):
            if start_params not in ['guess',  'zero']:
                raise ValueError("The argument start_params must be 'guess', "
                                 "'zero' or an array of correct length; "
                                 "got(start_params={0})".format(start_params))
        else:
            start_params = check_array(start_params, accept_sparse=False,
                                       force_all_finite=True, ensure_2d=False,
                                       dtype=_dtype, copy=True)
            if ((start_params.shape[0] != X.shape[1] + self.fit_intercept) or
                    (start_params.ndim != 1)):
                raise ValueError("Start values for parameters must have the"
                                 "right length and dimension; required (length"
                                 "={0}, ndim=1); got (length={1}, ndim={2})."
                                 .format(X.shape[1] + self.fit_intercept,
                                         start_params.shape[0],
                                         start_params.ndim))

        l1 = self.alpha * self.l1_ratio
        l2 = self.alpha * (1 - self.l1_ratio)
        # P1 and P2 are now for sure copies
        P1 = l1 * P1
        P2 = l2 * P2
        # one only ever needs the symmetrized L2 penalty matrix 1/2 (P2 + P2')
        # reason: w' P2 w = (w' P2 w)', i.e. it is symmetric
        if P2.ndim == 2:
            if sparse.issparse(P2):
                if sparse.isspmatrix_csc(P2):
                    P2 = 0.5 * (P2 + P2.transpose()).tocsc()
                else:
                    P2 = 0.5 * (P2 + P2.transpose()).tocsr()
            else:
                P2 = 0.5 * (P2 + P2.T)

        # For coordinate descent, if X is sparse, P2 must also be csc
        if solver == 'cd' and sparse.issparse(X):
            P2 = sparse.csc_matrix(P2)

        # 1.4 additional validations ##########################################
        if self.check_input:
            if not np.all(family.in_y_range(y)):
                raise ValueError("Some value(s) of y are out of the valid "
                                 "range for family {0}"
                                 .format(family.__class__.__name__))
            # check if P1 has only non-negative values, negative values might
            # indicate group lasso in the future.
            if not isinstance(self.P1, str):  # if self.P1 != 'identity':
                if not np.all(P1 >= 0):
                    raise ValueError("P1 must not have negative values.")
            # check if P2 is positive semidefinite
            # np.linalg.cholesky(P2) 'only' asserts positive definite
            if not isinstance(self.P2, str):  # self.P2 != 'identity'
                # due to numerical precision, we allow eigenvalues to be a
                # tiny bit negative
                epsneg = -10 * np.finfo(P2.dtype).epsneg
                if P2.ndim == 1 or P2.shape[0] == 1:
                    p2 = P2
                    if sparse.issparse(P2):
                        p2 = P2.toarray()
                    if not np.all(p2 >= 0):
                        raise ValueError("1d array P2 must not have negative "
                                         "values.")
                elif sparse.issparse(P2):
                    # for sparse matrices, not all eigenvals can be computed
                    # efficiently, use only half of n_features
                    # k = how many eigenvals to compute
                    k = np.min([10, n_features // 10 + 1])
                    sigma = -1000 * epsneg  # start searching near this value
                    which = 'SA'  # find smallest algebraic eigenvalues first
                    eigenvalues = splinalg.eigsh(P2, k=k, sigma=sigma,
                                                 which=which,
                                                 return_eigenvectors=False)
                    if not np.all(eigenvalues >= epsneg):
                        raise ValueError("P2 must be positive semi-definite.")
                else:
                    if not np.all(linalg.eigvalsh(P2) >= epsneg):
                        raise ValueError("P2 must be positive semi-definite.")
            # TODO: if alpha=0 check that X is not rank deficient
            # TODO: what else to check?

        #######################################################################
        # 2. rescaling of weights (sample_weight)                             #
        #######################################################################
        # IMPORTANT NOTE: Since we want to minimize
        # 1/(2*sum(sample_weight)) * deviance + L1 + L2,
        # deviance = sum(sample_weight * unit_deviance),
        # we rescale weights such that sum(weights) = 1 and this becomes
        # 1/2*deviance + L1 + L2 with deviance=sum(weights * unit_deviance)
        weights_sum = np.sum(weights)
        weights = weights/weights_sum

        #######################################################################
        # 3. initialization of coef = (intercept_, coef_)                     #
        #######################################################################
        # Note: Since phi=self.dispersion_ does not enter the estimation
        #       of mu_i=E[y_i], set it to 1.

        # set start values for coef
        coef = None
        if self.warm_start and hasattr(self, 'coef_'):
            if self.fit_intercept:
                coef = np.concatenate((np.array([self.intercept_]),
                                       self.coef_))
            else:
                coef = self.coef_
        elif isinstance(start_params, str):
            if start_params == 'guess':
                # Set mu=starting_mu of the family and do one Newton step
                # If solver=cd use cd, else irls
                mu = family.starting_mu(y, weights=weights)
                eta = link.link(mu)  # linear predictor
                if solver in ['cd', 'lbfgs', 'newton-cg']:
                    # see function _cd_solver
                    sigma_inv = 1/family.variance(mu, phi=1, weights=weights)
                    d1 = link.inverse_derivative(eta)
                    temp = sigma_inv * d1 * (y - mu)
                    if self.fit_intercept:
                        score = np.concatenate(([temp.sum()], temp @ X))
                    else:
                        score = temp @ X  # same as X.T @ temp

                    d2_sigma_inv = d1 * d1 * sigma_inv
                    diag_fisher = self.diag_fisher
                    if diag_fisher:
                        fisher = d2_sigma_inv
                    else:
                        fisher = \
                            _safe_sandwich_dot(X, d2_sigma_inv,
                                               intercept=self.fit_intercept)
                    # set up space for search direction d for inner loop
                    if self.fit_intercept:
                        coef = np.zeros(n_features+1)
                    else:
                        coef = np.zeros(n_features)
                    d = np.zeros_like(coef)
                    # initial stopping tolerance of inner loop
                    # use L1-norm of minimum of norm of subgradient of F
                    # use less restrictive tolerance for initial guess
                    inner_tol = _min_norm_sugrad(coef=coef, grad=-score, P2=P2,
                                                 P1=P1)
                    inner_tol = 4 * linalg.norm(inner_tol, ord=1)
                    # just one outer loop = Newton step
                    n_cycles = 0
                    d, coef_P2, n_cycles, inner_tol = \
                        _cd_cycle(d, X, coef, score, fisher, P1, P2, n_cycles,
                                  inner_tol, max_inner_iter=1000,
                                  selection=self.selection,
                                  random_state=random_state,
                                  diag_fisher=self.diag_fisher)
                    coef += d  # for simplicity no line search here
                else:
                    # See _irls_solver
                    # h'(eta)
                    hp = link.inverse_derivative(eta)
                    # working weights W, in principle a diagonal matrix
                    # therefore here just as 1d array
                    W = (hp**2 / family.variance(mu, phi=1, weights=weights))
                    # working observations
                    z = eta + (y-mu)/hp
                    # solve A*coef = b
                    # A = X' W X + l2 P2, b = X' W z
                    coef = _irls_step(X, W, P2, z,
                                      fit_intercept=self.fit_intercept)
            else:  # start_params == 'zero'
                if self.fit_intercept:
                    coef = np.zeros(n_features+1)
                    coef[0] = link.link(np.average(y, weights=weights))
                else:
                    coef = np.zeros(n_features)
        else:  # assign given array as start values
            coef = start_params

        #######################################################################
        # 4. fit                                                              #
        #######################################################################
        # algorithms for optimization
        # TODO: Parallelize it?

        # 4.1 IRLS ############################################################
        # Note: we already set P2 = l2*P2, see above
        # Note: we already symmetrized P2 = 1/2 (P2 + P2')
        if solver == 'irls':
            coef, self.n_iter_ = \
                _irls_solver(coef=coef, X=X, y=y, weights=weights, P2=P2,
                             fit_intercept=self.fit_intercept, family=family,
                             link=link, max_iter=self.max_iter, tol=self.tol)

        # 4.2 L-BFGS ##########################################################
        elif solver == 'lbfgs':
            def func(coef, X, y, weights, P2, family, link):
                mu, devp = \
                    family._mu_deviance_derivative(coef, X, y, weights, link)
                dev = family.deviance(y, mu, weights)
                intercept = (coef.size == X.shape[1] + 1)
                idx = 1 if intercept else 0  # offset if coef[0] is intercept
                if P2.ndim == 1:
                    L2 = P2 * coef[idx:]
                else:
                    L2 = P2 @ coef[idx:]
                obj = 0.5 * dev + 0.5 * (coef[idx:] @ L2)
                objp = 0.5 * devp
                objp[idx:] += L2
                return obj, objp

            args = (X, y, weights, P2, family, link)
            coef, loss, info = fmin_l_bfgs_b(
                func, coef, fprime=None, args=args,
                iprint=(self.verbose > 0) - 1, pgtol=self.tol,
                maxiter=self.max_iter, factr=1e3)
            if info["warnflag"] == 1:
                warnings.warn("lbfgs failed to converge."
                              " Increase the number of iterations.",
                              ConvergenceWarning)
            elif info["warnflag"] == 2:
                warnings.warn("lbfgs failed for the reason: {0}"
                              .format(info["task"]))
            self.n_iter_ = info['nit']

        # 4.3 Newton-CG #######################################################
        # We use again the fisher matrix instead of the hessian. More
        # precisely, expected hessian of deviance.
        elif solver == 'newton-cg':
            def func(coef, X, y, weights, P2, family, link):
                intercept = (coef.size == X.shape[1] + 1)
                idx = 1 if intercept else 0  # offset if coef[0] is intercept
                if P2.ndim == 1:
                    L2 = coef[idx:] @ (P2 * coef[idx:])
                else:
                    L2 = coef[idx:] @ (P2 @ coef[idx:])
                mu = link.inverse(_safe_lin_pred(X, coef))
                return 0.5 * family.deviance(y, mu, weights) + 0.5 * L2

            def grad(coef, X, y, weights, P2, family, link):
                mu, devp = \
                    family._mu_deviance_derivative(coef, X, y, weights, link)
                intercept = (coef.size == X.shape[1] + 1)
                idx = 1 if intercept else 0  # offset if coef[0] is intercept
                if P2.ndim == 1:
                    L2 = P2 * coef[idx:]
                else:
                    L2 = P2 @ coef[idx:]
                objp = 0.5 * devp
                objp[idx:] += L2
                return objp

            def grad_hess(coef, X, y, weights, P2, family, link):
                intercept = (coef.size == X.shape[1] + 1)
                idx = 1 if intercept else 0  # offset if coef[0] is intercept
                if P2.ndim == 1:
                    L2 = P2 * coef[idx:]
                else:
                    L2 = P2 @ coef[idx:]
                eta = _safe_lin_pred(X, coef)
                mu = link.inverse(eta)
                d1 = link.inverse_derivative(eta)
                temp = d1 * family.deviance_derivative(y, mu, weights)
                if intercept:
                    grad = np.concatenate(([0.5 * temp.sum()],
                                           0.5 * temp @ X + L2))
                else:
                    grad = 0.5 * temp @ X + L2  # same as 0.5* X.T @ temp + L2

                # expected hessian = fisher = X.T @ diag_matrix @ X
                # calculate only diag_matrix
                diag = d1**2 / family.variance(mu, phi=1, weights=weights)
                if intercept:
                    h0i = np.concatenate(([diag.sum()], diag @ X))

                def Hs(coef):
                    # return (0.5 * fisher + P2) @ coef
                    # ret = 0.5 * (X.T @ (diag * (X @ coef)))
                    ret = 0.5 * ((diag * (X @ coef[idx:])) @ X)
                    if P2.ndim == 1:
                        ret += P2 * coef[idx:]
                    else:
                        ret += P2 @ coef[idx:]
                    if intercept:
                        ret = np.concatenate(([0.5 * (h0i @ coef)],
                                             ret + 0.5 * coef[0] * h0i[1:]))
                    return ret

                return grad, Hs

            args = (X, y, weights, P2, family, link)
            coef, self.n_iter_ = newton_cg(grad_hess, func, grad, coef,
                                           args=args, maxiter=self.max_iter,
                                           tol=self.tol)

        # 4.4 coordinate descent ##############################################
        # Note: we already set P1 = l1*P1, see above
        # Note: we already set P2 = l2*P2, see above
        # Note: we already symmetrized P2 = 1/2 (P2 + P2')
        elif solver == 'cd':
            coef, self.n_iter_, self._n_cycles = \
                _cd_solver(coef=coef, X=X, y=y, weights=weights, P1=P1,
                           P2=P2, fit_intercept=self.fit_intercept,
                           family=family, link=link,
                           max_iter=self.max_iter, tol=self.tol,
                           selection=self.selection, random_state=random_state,
                           diag_fisher=self.diag_fisher, copy_X=self.copy_X)

        #######################################################################
        # 5. postprocessing                                                   #
        #######################################################################
        if self.fit_intercept:
            self.intercept_ = coef[0]
            self.coef_ = coef[1:]
        else:
            # set intercept to zero as the other linear models do
            self.intercept_ = 0.
            self.coef_ = coef

        if self.fit_dispersion in ['chisqr', 'deviance']:
            # attention because of rescaling of weights
            self.dispersion_ = self.estimate_phi(X, y, weights)*weights_sum

        return self

    def linear_predictor(self, X):
        """Compute the linear_predictor = X*coef_ + intercept_.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        C : array, shape (n_samples,)
            Returns predicted values of linear predictor.
        """
        check_is_fitted(self, "coef_")
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'],
                        dtype='numeric', copy=True, ensure_2d=True,
                        allow_nd=False)
        return X @ self.coef_ + self.intercept_

    def predict(self, X, sample_weight=None):
        """Predict using GLM with feature matrix X.

        If sample_weight is given, returns prediction*sample_weight.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Samples.

        sample_weight : {None, array-like}, shape (n_samples,), optional \
                (default=None)

        Returns
        -------
        C : array, shape (n_samples,)
            Returns predicted values times sample_weight.
        """
        # TODO: Is copy=True necessary?
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'],
                        dtype='numeric', copy=True, ensure_2d=True,
                        allow_nd=False)
        eta = self.linear_predictor(X)
        mu = self._link_instance.inverse(eta)
        weights = _check_weights(sample_weight, X.shape[0])

        return mu*weights

    def estimate_phi(self, X, y, sample_weight=None):
        """Estimate/fit the dispersion parameter phi.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : {None, array-like}, shape (n_samples,), optional \
                (default=None)
            Sample weights.

        Returns
        -------
        phi : float
            Dispersion parameter.
        """
        check_is_fitted(self, "coef_")
        _dtype = [np.float64, np.float32]
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'],
                         dtype=_dtype, y_numeric=True, multi_output=False)
        n_samples, n_features = X.shape
        weights = _check_weights(sample_weight, n_samples)
        eta = X @ self.coef_
        if self.fit_intercept is True:
            eta += self.intercept_
            n_features += 1
        if n_samples <= n_features:
            raise ValueError("Estimation of dispersion parameter phi requires"
                             " more samples than features, got"
                             " samples=X.shape[0]={0} and"
                             " n_features=X.shape[1]+fit_intercept={1}."
                             .format(n_samples, n_features))
        mu = self._link_instance.inverse(eta)
        if self.fit_dispersion == 'chisqr':
            chisq = np.sum(weights*(y-mu)**2 /
                           self._family_instance.unit_variance(mu))
            return chisq/(n_samples - n_features)
        elif self.fit_dispersion == 'deviance':
            dev = self._family_instance.deviance(y, mu, weights)
            return dev/(n_samples - n_features)

    # Note: check_estimator(GeneralizedLinearRegressor) might raise
    # "AssertionError: -0.28014056555724598 not greater than 0.5"
    # unless GeneralizedLinearRegressor has a score which passes the test.
    def score(self, X, y, sample_weight=None):
        """Compute D^2, the percentage of deviance explained.

        D^2 is a generalization of the coefficient of determination R^2.
        R^2 uses squared error and D^2 deviance. Note that those two are equal
        for family='normal'.

        D^2 is defined as
        :math:`D^2 = 1-\\frac{D(y_{true},y_{pred})}{D_{null}}`,
        :math:`D_{null}` is the null deviance, i.e. the deviance of a model
        with intercept alone, which corresponds to :math:`y_{pred} = \\bar{y}`.
        The mean :math:`\\bar{y}` is averaged by sample_weight.
        Best possible score is 1.0 and it can be negative (because the model
        can be arbitrarily worse).

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Test samples.

        y : array-like, shape (n_samples,)
            True values of target.

        sample_weight : {None, array-like}, shape (n_samples,), optional \
                (default=None)
            Sample weights.

        Returns
        -------
        score : float
            D^2 of self.predict(X) w.r.t. y.
        """
        # Note, default score defined in RegressorMixin is R^2 score.
        # TODO: make D^2 a score function in module metrics (and thereby get
        #       input validation and so on)
        weights = _check_weights(sample_weight, y.shape[0])
        mu = self.predict(X)
        dev = self._family_instance.deviance(y, mu, weights=weights)
        y_mean = np.average(y, weights=weights)
        dev_null = self._family_instance.deviance(y, y_mean, weights=weights)
        return 1. - dev / dev_null

    def _more_tags(self):
        return {"requires_positive_y": True}


class PoissonRegressor(GeneralizedLinearRegressor):
    """Regression with the response variable y following a Poisson distribution

    GLMs based on a reproductive Exponential Dispersion Model (EDM) aim at
    fitting and predicting the mean of the target y as mu=h(X*w).
    The fit minimizes the following objective function with L2 regularization::

            1/(2*sum(s)) * deviance(y, h(X*w); s) + 1/2 * alpha * ||w||_2^2

    with inverse link function h and s=sample_weight. Note that for
    ``sample_weight=None``, one has s_i=1 and sum(s)=n_samples).

    Read more in the :ref:`User Guide <Generalized_linear_regression>`.

    Parameters
    ----------
    alpha : float, optional (default=1)
        Constant that multiplies the penalty terms and thus determines the
        regularization strength.
        See the notes for the exact mathematical meaning of this
        parameter.``alpha = 0`` is equivalent to unpenalized GLMs. In this
        case, the design matrix X must have full column rank
        (no collinearities).

    fit_intercept : boolean, optional (default=True)
        Specifies if a constant (a.k.a. bias or intercept) should be
        added to the linear predictor (X*coef+intercept).

    fit_dispersion : {None, 'chisqr', 'deviance'}, optional (default=None)
        Method for estimation of the dispersion parameter phi. Whether to use
        the chi squared statistic or the deviance statistic. If None, the
        dispersion is not estimated.

    solver : {'irls', 'lbfgs', 'newton-cg'}, optional (default='irls')
        Algorithm to use in the optimization problem:

        'irls'
            Iterated reweighted least squares. It is the standard algorithm
            for GLMs.

        'lbfgs'
            Calls scipy's L-BFGS-B optimizer.

        'newton-cg'
            Newton conjugate gradient algorithm.

        Note that all solvers except lbfgs use the fisher matrix, i.e. the
        expected Hessian instead of the Hessian matrix.

    max_iter : int, optional (default=100)
        The maximal number of iterations for solver algorithms.

    tol : float, optional (default=1e-4)
        Stopping criterion. For the irls, newton-cg and lbfgs solvers,
        the iteration will stop when ``max{|g_i|, i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative) of
        the objective function.

    warm_start : boolean, optional (default=False)
        If set to ``True``, reuse the solution of the previous call to ``fit``
        as initialization for ``coef_`` and ``intercept_`` (supersedes option
        ``start_params``). If set to ``True`` or if the attribute ``coef_``
        does not exit (first call to ``fit``), option ``start_params`` sets the
        start values for ``coef_`` and ``intercept_``.

    start_params : {'guess', 'zero', array of shape (n_features*, )}, \
            optional (default='guess')
        Relevant only if ``warm_start=False`` or if fit is called
        the first time (``self.coef_`` does not yet exist).

        'guess'
            Start values of mu are calculated by family.starting_mu(..). Then,
            one Newton step obtains start values for ``coef_``. If
            ``solver='irls'``, it uses one irls step. This gives usually good
            starting values.

        'zero'
        All coefficients are set to zero. If ``fit_intercept=True``, the
        start value for the intercept is obtained by the weighted average of y.

        array
        The array of size n_features* is directly used as start values
        for ``coef_``. If ``fit_intercept=True``, the first element
        is assumed to be the start value for the ``intercept_``.
        Note that n_features* = X.shape[1] + fit_intercept, i.e. it includes
        the intercept in counting.

    random_state : {int, RandomState instance, None}, optional (default=None)
        If int, random_state is the seed used by the random
        number generator; if RandomState instance, random_state is the random
        number generator; if None, the random number generator is the
        RandomState instance used by `np.random`. Used when ``selection`` ==
        'random'.

    copy_X : boolean, optional, (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    verbose : int, optional (default=0)
        For the lbfgs solver set verbose to any positive number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (n_features,)
        Estimated coefficients for the linear predictor (X*coef_+intercept_) in
        the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    dispersion_ : float
        The dispersion parameter :math:`\\phi` if ``fit_dispersion`` was set.

    n_iter_ : int
        Actual number of iterations used in solver.

    Notes
    -----
    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments to be :math:`E[Y_i]=\\mu_i=h((Xw)_i)` and
    :math:`Var[Y_i]=\\frac{\\phi}{s_i} v(\\mu_i)`. The unit variance function
    :math:`v(\\mu_i)` is a property of and given by the specific EDM, see
    :ref:`User Guide <Generalized_linear_regression>`.

    The parameters :math:`w` (`coef_` and `intercept_`) are estimated by
    minimizing the deviance plus penalty term, which is equivalent to
    (penalized) maximum likelihood estimation.

    For alpha > 0, the feature matrix X should be standardized in order to
    penalize features equally strong.

    If the target y is a ratio, appropriate sample weights s should be
    provided.
    As an example, consider Poisson distributed counts z (integers) and
    weights s=exposure (time, money, persons years, ...). Then you fit
    y = z/s, i.e. ``PoissonRegressor().fit(X, y, sample_weight=s)``.
    The weights are necessary for the right (finite sample) mean.
    Consider :math:`\\bar{y} = \\frac{\\sum_i s_i y_i}{\\sum_i s_i}`,
    in this case one might say that y has a 'scaled' Poisson distributions.

    References
    ----------
    For the coordinate descent implementation:
        * Guo-Xun Yuan, Chia-Hua Ho, Chih-Jen Lin
          An Improved GLMNET for L1-regularized Logistic Regression,
          Journal of Machine Learning Research 13 (2012) 1999-2030
          https://www.csie.ntu.edu.tw/~cjlin/papers/l1_glmnet/long-glmnet.pdf
    """
    def __init__(self, alpha=1.0, fit_intercept=True, fit_dispersion=None,
                 solver='irls', max_iter=100,
                 tol=1e-4, warm_start=False, start_params='guess',
                 random_state=None, copy_X=True, check_input=True, verbose=0):

        super().__init__(alpha=alpha, fit_intercept=fit_intercept,
                         family="poisson", link='log',
                         fit_dispersion=fit_dispersion, solver=solver,
                         max_iter=max_iter, tol=tol, warm_start=warm_start,
                         start_params=start_params, random_state=random_state,
                         copy_X=copy_X, verbose=verbose)
