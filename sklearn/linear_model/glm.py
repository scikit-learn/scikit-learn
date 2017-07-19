"""
Generalized Linear Models with Exponential Dispersion Family
"""

# Author: Christian Lorentzen <lorentzen.ch@googlemail.ch>
# License: BSD 3 clause

# TODO: Which name/symbol for coefficients and weights in docu?
#       sklearn.linear_models uses w for coefficients.
#       So far, coefficients=beta and weight=w (as standard literature)
# TODO: Add l2-penalty
# TODO: Add l1-penalty (elastic net)
# TODO: Add cross validation
# TODO: Write docu and examples

# Design Decisions:
# - Which name? GeneralizedLinearModel vs GeneralizedLinearRegressor.
#   So far, it is GeneralizedLinearModel, since it could very easily
#   extended by Bernoulli/Binomial distribution.
#   Solution: GeneralizedLinearRegressor
# - The link funtion (instance of class Link) is necessary for the evaluation
#   of deviance, score, Fisher and Hessian matrix as functions of the
#   coefficients, which is needed by optimizers.
#   Solution: link as argument in those functions

from __future__ import division
from abc import ABCMeta, abstractmethod, abstractproperty
import numbers
import numpy as np
from scipy import linalg, optimize, sparse
import warnings
from .base import LinearRegression
from ..base import BaseEstimator, RegressorMixin
from ..exceptions import ConvergenceWarning
from ..externals import six
from ..utils import check_X_y
from ..utils.extmath import safe_sparse_dot
from ..utils.optimize import newton_cg
from ..utils.validation import check_is_fitted


class Link(six.with_metaclass(ABCMeta)):
    """Abstract base class for Link funtions
    """

    @abstractmethod
    def link(self, mu):
        """The link function g(mu) with argument mu=E[Y] returns the
        linear predictor.
        """
        raise NotImplementedError

    @abstractmethod
    def derivative(self, mu):
        """Derivative of the link g'(mu).
        """
        raise NotImplementedError

    @abstractmethod
    def inverse(self, lin_pred):
        """The inverse link function h(lin_pred) with the linear predictor as
        argument returns mu=E[Y].
        """
        raise NotImplementedError

    @abstractmethod
    def inverse_derivative(self, lin_pred):
        """Derivative of the inverse link function h'(lin_pred).
        """
        raise NotImplementedError

    @abstractmethod
    def inverse_derivative2(self, lin_pred):
        """Second derivative of the inverse link function h''(lin_pred).
        """
        raise NotImplementedError


class IdentityLink(Link):
    """The identity link function g(x)=x.
    """

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
    """The log link function g(x)=log(x).
    """

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


class ExponentialDispersionModel(six.with_metaclass(ABCMeta)):
    """Base class for reproductive Exponential Dispersion Models (EDM).

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

    _score
    _fisher_matrix
    _observed_information
    _deviance
    _deviance_derivative
    _deviance_hessian

    References
    ----------
    See https://en.wikipedia.org/wiki/Exponential_dispersion_model.
    """

    @abstractproperty
    def lower_bound(self):
        """The lower bound of values of Y~EDM.
        """
        raise NotImplementedError()

    @abstractproperty
    def upper_bound(self):
        """The upper bound of values of Y~EDM.
        """
        raise NotImplementedError()

    @abstractmethod
    def in_y_range(self, x):
        """Returns true if x is in the valid range of Y~EDM.
        """
        raise NotImplementedError()

    @abstractmethod
    def unit_variance(self, mu):
        """The unit variance :math:`v(mu)` determines the variance as
        a function of the mean mu by
        :math:`\mathrm{Var}[Y_i] = \phi/w_i*v(\mu_i)`.
        It can also be derived from the unit deviance :math:`d(y,\mu)` as

        .. math:: v(\mu) = \frac{2}{\frac{\partial^2 d(y,\mu)}{
            \partial\mu^2}}\big|_{y=\mu}
        """
        raise NotImplementedError()

    @abstractmethod
    def unit_variance_derivative(self, mu):
        """The derivative of the unit variance w.r.t. mu, :math:`v'(\mu)`.
        """
        raise NotImplementedError()

    def variance(self, mu, phi=1, weight=1):
        """The variance of :math:`Y \sim \mathrm{EDM}(\mu,\phi)` is
        :math:`\mathrm{Var}[Y_i]=\phi/w_i*v(\mu_i)`,
        with unit variance v(mu).
        """
        return phi/weight * self.unit_variance(mu)

    def variance_derivative(self, mu, phi=1, weight=1):
        """The derivative of the variance w.r.t. mu,
        :math:`\frac{\partial}{\partial\mu}\mathrm{Var}[Y_i]
        =phi/w_i*v'(\mu_i)`, with unit variance v(mu).
        """
        return phi/weight * self.unit_variance_derivative(mu)

    @abstractmethod
    def unit_deviance(self, y, mu):
        """The unit_deviance :math:`d(y,\mu)`.
        In terms of the log-likelihood it is given by
        :math:`d(y,\mu) = -2\phi\cdot
        \left(loglike(y,\mu,phi) - loglike(y,y,phi)\right).`
        """
        raise NotImplementedError()

    def unit_deviance_derivative(self, y, mu):
        """The derivative w.r.t. mu of the unit_deviance
        :math:`\frac{d}{d\mu}d(y,\mu) = -2\frac{y-\mu}{v(\mu)}`
        with unit variance :math:`v(\mu)`.

        Returns
        -------
        derivative: array, shape = (n_samples,)
        """
        return -2*(y-mu)/self.unit_variance(mu)

    def deviance(self, y, mu, weight=1):
        """The deviance is given by :math:`D = \sum_i w_i \cdot d(y, \mu)
        with weight :math:`w_i` and unit_deviance :math:`d(y,mu)`.
        In terms of the likelihood it is :math:`D = -2\phi\cdot
        \left(loglike(y,\mu,\frac{phi}{w})
        - loglike(y,y,\frac{phi}{w})\right).`
        """
        return np.sum(weight*self.unit_deviance(y, mu))

    def _deviance(self, coef, X, y, weight, link):
        """The deviance as a function of the coefficients ``coef``
        (:math:`beta`).
        """
        lin_pred = safe_sparse_dot(X, coef, dense_output=True)
        mu = link.inverse(lin_pred)
        return self.deviance(y, mu, weight)

    def deviance_derivative(self, y, mu, weight=1):
        """The derivative w.r.t. mu of the deviance.`
        """
        return weight*self.unit_deviance_derivative(y, mu)

    def _score(self, coef, phi, X, y, weight, link):
        """The score function :math:`s` is the derivative of the
        log-likelihood w.r.t. the ``coef`` (:math:`\beta`).
        It is given by

        .. math:

            \mathbf{s}(\boldsymbol{\beta}) = \mathbf{X}^T \mathbf{D}
            \boldsymbol{\Sigma}^-1 (\mathbf{y} - \boldsymbol{\mu})\,,

        with :math:`\mathbf{D}=\mathrm{diag}(h'(\eta_1),\ldots)` and
        :math:`\boldsymbol{\Sigma}=\mathrm{diag}(\mathbf{V}(y_1),\ldots)`.
        """
        n_samples = X.shape[0]
        lin_pred = safe_sparse_dot(X, coef, dense_output=True)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weight=weight)
        d = link.inverse_derivative(lin_pred)
        d_sigma_inv = sparse.dia_matrix((sigma_inv*d, 0),
                                        shape=(n_samples, n_samples))
        temp = safe_sparse_dot(d_sigma_inv, (y-mu), dense_output=False)
        score = safe_sparse_dot(X.T, temp, dense_output=False)
        return score

    def _fisher_matrix(self, coef, phi, X, y, weight, link):
        """The Fisher information matrix, also known as expected
        information matrix. It is given by

        .. math:

            \mathbf{F}(\boldsymbol{\beta}) = \mathrm{E}\left[
            -\frac{\partial^2 loglike}{\partial\boldsymbol{\beta}
            \partial\boldsymbol{\beta}^T}\right]
            = \mathbf{X}^T W \mathbf{X} \,,

        with :math:`\mathbf{W} = \mathbf{D}^2 \boldsymbol{\Sigma}^{-1}`,
        see score function.
        """
        n_samples = X.shape[0]
        lin_pred = safe_sparse_dot(X, coef, dense_output=True)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weight=weight)
        d2 = link.inverse_derivative(lin_pred)**2
        d2_sigma_inv = sparse.dia_matrix((sigma_inv*d2, 0),
                                         shape=(n_samples, n_samples))
        temp = safe_sparse_dot(d2_sigma_inv, X, dense_output=False)
        fisher_matrix = safe_sparse_dot(X.T, temp, dense_output=False)
        return fisher_matrix

    def _observed_information(self, coef, phi, X, y, weight, link):
        """The observed information matrix, also known as the negative of
        the Hessian matrix of the log-likelihood. It is given by

        .. math:

            \mathbf{H}(\boldsymbol{\beta}) =
            -\frac{\partial^2 loglike}{\partial\boldsymbol{\beta}
            \partial\boldsymbol{\beta}^T}
            = \mathbf{X}^T \legt[
            - \mathbf{D}' \mathbf{R}
            + \mathbf{D}^2 \mathbf{V} \mathbf{R}
            + \mathbf{D}^2
            \right] \boldsymbol{\Sigma}^{-1} \mathbf{X} \,,

        with :math:`\mathbf{R} = \mathrm{diag}(y_i - \mu_i)`,
        :math:`\mathbf{V} = \mathrm{diag}\left(\frac{v'(\mu_i)}{
        v(\mu_i)}
        \right)`,
        see score function and Fisher matrix.
        """
        n_samples = X.shape[0]
        lin_pred = safe_sparse_dot(X, coef, dense_output=True)
        mu = link.inverse(lin_pred)
        sigma_inv = 1/self.variance(mu, phi=phi, weight=weight)
        dp = link.inverse_derivative2(lin_pred)
        d2 = link.inverse_derivative(lin_pred)**2
        v = self.unit_variance_derivative(mu)/self.unit_variance(mu)
        r = y - mu
        temp = sparse.dia_matrix((sigma_inv*(-dp*r+d2*v*r+d2), 0),
                                 shape=(n_samples, n_samples))
        temp = safe_sparse_dot(temp, X, dense_output=False)
        observed_information = safe_sparse_dot(X.T, temp, dense_output=False)
        return observed_information

    def _deviance_derivative(self, coef, X, y, weight, link):
        """The derivative w.r.t. ``coef`` (:math:`\beta`) of the deviance as a
        function of the coefficients ``coef``.
        This is equivalent to :math:`-2\phi` times the score function
        :math:`s` (derivative of the log-likelihood).
        """
        score = self._score(coef=coef, phi=1, X=X, y=y, weight=weight,
                            link=link)
        return -2*score

    def _deviance_hessian(self, coef, X, y, weight, link):
        """The hessian matrix w.r.t. ``coef`` (:math:`\beta`) of the deviance
        as a function of the coefficients ``coef``.
        This is equivalent to :math:`+2\phi` times the observed information
        matrix.
        """
        info_matrix = self._observed_information(coef=coef, phi=1, X=X, y=y,
                                                 weight=weight, link=link)
        return 2*info_matrix

    def starting_mu(self, y, weight=1):
        """Starting values for the mean mu_i in IRLS."""
        return (weight*y+np.mean(weight*y))/(2.*np.sum(np.ones_like(y)*weight))


class TweedieDistribution(ExponentialDispersionModel):
    """A class for the Tweedie distribution.
    They have mu=E[X] and Var[X] \propto mu**power.

    Attributes
    ----------
    power : float
            The variance power of the unit_variance
            :math:`v(mu) = mu^{power}`.
    """
    def __init__(self, power=0):
        self.power = power
        self._upper_bound = np.Inf
        self._upper_compare = lambda x: np.less(x, self.upper_bound)
        if power < 0:
            # Extreme Stable
            self._lower_bound = -np.Inf
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)
        elif power == 0:
            # GaussianDistribution
            self._lower_bound = -np.Inf
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)
        elif (power > 0) and (power < 1):
            raise ValueError('For 0<power<1, no distribution exists.')
        elif power == 1:
            # PoissonDistribution
            self._lower_bound = 0
            self._lower_compare = (
                lambda x: np.greater_equal(x, self.lower_bound))
        elif (power > 1) and (power < 2):
            # Compound Poisson
            self._lower_bound = 0
            self._lower_compare = (
                lambda x: np.greater_equal(x, self.lower_bound))
        elif power == 2:
            # GammaDistribution
            self._lower_bound = 0
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)
        elif (power > 2) and (power < 3):
            # Positive Stable
            self._lower_bound = 0
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)
        elif power == 3:
            # InverseGaussianDistribution
            self._lower_bound = 0
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)
        elif power > 3:
            # Positive Stable
            self._lower_bound = 0
            self._lower_compare = lambda x: np.greater(x, self.lower_bound)

    @property
    def power(self):
        return self._power

    @power.setter
    def power(self, power):
        if not isinstance(power, numbers.Real):
            raise TypeError('power must be a real number, input was {0}'
                            .format(power))
        self._power = power

    @property
    def lower_bound(self):
        return self._lower_bound

    @property
    def upper_bound(self):
        return self._upper_bound

    def in_y_range(self, x):
        return np.logical_and(self._lower_compare(x), self._upper_compare(x))

    def unit_variance(self, mu):
        """The unit variance of a Tweedie distribution is v(mu)=mu**power.
        """
        return np.power(mu, self.power)

    def unit_variance_derivative(self, mu):
        """The derivative of the unit variance of a Tweedie distribution is
        v(mu)=power*mu**(power-1).
        """
        return self.power*np.power(mu, self.power-1)

    def unit_deviance(self, y, mu):
        p = self.power
        if p == 0:
            # NormalDistribution
            return (y-mu)**2
        if p == 1:
            # PoissonDistribution
            return 2 * (np.where(y == 0, 0, y*np.log(y/mu))-y+mu)
        elif p == 2:
            # GammaDistribution
            return 2 * (np.log(mu/y)+y/mu-1)
        else:
            # return 2 * (np.maximum(y,0)**(2-p)/((1-p)*(2-p))
            #    - y*mu**(1-p)/(1-p) + mu**(2-p)/(2-p))
            return 2 * (np.power(np.maximum(y, 0), 2-p)/((1-p)*(2-p)) -
                        y*np.power(mu, 1-p)/(1-p) + np.power(mu, 2-p)/(2-p))

    def likelihood(self, y, X, beta, phi, weight=1):
        raise NotImplementedError('This function is not (yet) implemented.')


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


class GeneralizedHyperbolicSecand(ExponentialDispersionModel):
    """A class for the von Generalized Hyperbolic Secand (GHS) distribution.

    The GHS distribution is for data y in (-inf, inf).
    """
    def __init__(self):
        self._lower_bound = -np.Inf
        self._upper_bound = np.Inf

    @property
    def lower_bound(self):
        return self._lower_bound

    @property
    def upper_bound(self):
        return self._upper_bound

    def in_y_range(self, x):
        np.logical_and(
            np.greater(x, self.lower_bound),
            np.less(x, self.lower_bound)
            )

    def unit_variance(self, mu):
        return 1 + mu**2

    def unit_variance_derivative(self, mu):
        return 2*mu

    def unit_deviance(self, y, mu):
        return (2*y*(np.arctan(y) - np.arctan(mu)) +
                np.log((1+mu**2)/(1+y**2)))


class GeneralizedLinearRegressor(BaseEstimator, RegressorMixin):
    """
    Class to fit a Generalized Linear Model (GLM) based on reproductive
    Exponential Dispersion Models (EDM).

    Assumptions:

    - The target values y_i are realizations of random variables
      :math:`Y_i \sim \mathrm{EDM}(\mu_i, \frac{\phi}{w_i})` with dispersion
      parameter :math:`\phi` and weights :math:`w_i`.
    - The expectation of :math:`Y_i` is :math:`\mu_i=\mathrm{E}[Y]=h(\eta_i)`
      whith the linear predictor :math:`\eta=X*\beta`, inverse link function
      :math:`h(\eta)`, design matrix :math:`X` and parameters :math:`\beta`
      to be estimated.

    Note that the first assumption implies
    :math:`\mathrm{Var}[Y_i]=\frac{\phi}{w_i} v(\mu_i)` with uni variance
    function :math:`v(\mu)`.

    The fit itself does not need Y to be from an EDM, but only assumes
    the first two moments :math:`E[Y_i]=\mu_i=h(\eta_i)` and
    :math:`Var[Y_i]=\frac{\phi}{w_i} v(\mu_i)`

    The parameters :math:`\beta` are estimated by maximum likelihood which is
    equivalent to minimizing the deviance.

    TODO: Estimation of the dispersion parameter phi.

    TODO: Notes on weights and 'scaled' Poisson, e.g. fit y = x/w with
    with x=counts and w=exposure (time, money, persons, ...) => y is a
    ratio with weights w.

    Parameters
    ----------
    fit_intercept : boolean, optional, default True
        whether to calculate the intercept for this model. If set
        to False, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    family : {'normal', 'poisson', 'gamma', 'inverse.gaussian'} or an instance
        of a subclass of ExponentialDispersionModel, optional, default 'normal'
        the distributional assumption of the GLM.

    link : {'identity', 'log'} or an instance of a subclass of Link,
        optional, default IdentityLink()
        the link function (class) of the GLM

    fit_dispersion : {None, 'chisqr', 'deviance'}, defaul 'chisqr'
        method for estimation of the dispersion parameter phi. Whether to use
        the chi squared statisic or the deviance statistic. If None, the
        dispersion is not estimated.

    solver : {'irls', 'newton-cg', 'lbfgs'}, defaul 'irls'
        Algorithm to use in the optimization problem.

        - 'irls' is iterated reweighted least squares. It is the standard
            algorithm for GLMs.

        - 'newton-cg', 'lbfgs'

    max_iter : int, default 100
        TODO

    tol : float
        Stopping criterion. For the irls, newton-cg and lbfgs solvers,
        the iteration will stop when ``max{|g_i | i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient (derivative of
        the deviance).

    start_params : {array shape (n_features, ), 'ols'}, default None
        sets the start values for coef_ in the fit.
        If None, default values are taken.
        If 'ols' the result of an ordinary least squares in the link space
        (linear predictor) is taken.
        If an array is given, these values are taken as coef_ to start with.
        If fit_intercept is true, the first value is assumed to be the start
        value for the intercept_.

    verbose : int, default: 0
        For the lbfgs solver set verbose to any positive
        number for verbosity.

    Attributes
    ----------
    coef_ : array, shape (1, n_features)
        Estimated coefficients for the linear predictor (X*coef_) in the GLM.

    intercept_ : float
        Intercept (a.k.a. bias) added to linear predictor.

    dispersion_ : float
        The dispersion parameter :math:`\phi` if fit_dispersion is set.

    n_iter_ : int
        Actual number of iterations of the solver.

    Notes
    -----

    References
    ----------
    TODO
    """

    def __init__(self, fit_intercept=True, family=NormalDistribution(),
                 link=IdentityLink(), fit_dispersion='chisqr', solver='irls',
                 max_iter=100, tol=1e-4, start_params=None, verbose=0):
        self.fit_intercept = fit_intercept
        self.family = family
        self.link = link
        self.fit_dispersion = fit_dispersion
        self.solver = solver
        self.max_iter = 100
        self.tol = tol
        self.start_params = start_params
        self.verbose = verbose

    def fit(self, X, y, weight=None):
        """Fit a generalized linear model.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        weight : numpy array of shape [n_samples]
            Individual weights for each sample.
            Var[Y_i]=phi/weight_i * v(mu)
            If Y_i ~ EDM(mu, phi/w_i) then
            sum(w*Y)/sum(w) ~ EDM(mu, phi/sum(w))

        Returns
        -------
        self : returns an instance of self.
        """
        if not isinstance(self.family, ExponentialDispersionModel):
            if self.family == 'normal':
                self.family = NormalDistribution()
            elif self.family == 'poisson':
                self.family = PoissonDistribution()
            elif self.family == 'gamma':
                self.family = GammaDistribution()
            elif self.family == 'inverse.gaussian':
                self.family = InverseGaussianDistribution()
            else:
                raise ValueError(
                    "The argument family must be an instance of class"
                    " ExponentialDispersionModel or an element of"
                    " ['normal', 'poisson', 'gamma', 'inverse.gaussian'].")
        if not isinstance(self.link, Link):
            if self.link == 'identity':
                self.link = IdentityLink()
            if self.link == 'log':
                self.link = LogLink()
            else:
                raise ValueError(
                    "The argument link must be an instance of class Link or"
                    " an element of ['identity', 'log'].")
        if not isinstance(self.fit_intercept, bool):
            raise ValueError("The argument fit_intercept must be bool,"
                             " got {0}".format(self.fit_intercept))
        if self.solver not in ['irls', 'lbfgs', 'newton-cg']:
            raise ValueError("GLM Regression supports only irls, lbfgs and"
                             "newton-cg solvers, got {0}".format(self.solver))
        if not isinstance(self.max_iter, numbers.Number) or self.max_iter < 0:
            raise ValueError("Maximum number of iteration must be positive;"
                             " got (max_iter={0!r})".format(self.max_iter))
        if not isinstance(self.tol, numbers.Number) or self.tol < 0:
            raise ValueError("Tolerance for stopping criteria must be "
                             "positive; got (tol={0!r})".format(self.tol))
        start_params = self.start_params
        if start_params is not None and start_params is not 'ols':
            start_params = np.atleast_1d(start_params)
            if start_params.shape[0] != X.shape[1] + self.fit_intercept:
                raise ValueError("Start values for parameters must have the"
                                 "right length; required length {0}, got {1}"
                                 .format(X.shape[1] + self.fit_intercept,
                                         start_params.shape[0]))

        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'],
                         y_numeric=True, multi_output=False)
        y = y.astype(np.float64)

        if not np.all(self.family.in_y_range(y)):
            raise ValueError("Some value(s) of y are out of the valid "
                             "range for family {0}"
                             .format(self.family.__class__.__name__))

        if weight is None:
            weight = np.ones_like(y)
        elif np.isscalar(weight):
            weight = weight*np.ones_like(y)
        else:
            weight = np.atleast_1d(weight)
            if weight.ndim > 1:
                raise ValueError("Weights must be 1D array or scalar")
            elif weight.shape[0] != y.shape[0]:
                raise ValueError("Weights must have the same length as y")

        if self.fit_intercept:
            # intercept is first column <=> coef[0] is for intecept
            if sparse.issparse(X):
                Xnew = sparse.hstack([np.ones([X.shape[0], 1]), X])
            else:
                Xnew = np.concatenate((np.ones((X.shape[0], 1)), X), axis=1)
        else:
            Xnew = X

        n_samples, n_features = Xnew.shape

        # Note: Since dispersion_ alias phi does not enter the estimation
        #       of mu_i=E[y_i] set it to 1 where convenient.

        # set start values for coef
        coef = None
        if start_params is None:
            # Use mu_start and apply one irls step to calculate coef
            mu = self.family.starting_mu(y, weight)
            # linear predictor
            eta = self.link.link(mu)
            # h'(eta)
            hp = self.link.inverse_derivative(eta)
            # working weights w, in principle a diagonal matrix
            # therefore here just as 1d array
            w = (hp**2 / self.family.variance(mu, phi=1, weight=weight))
            wroot = np.sqrt(w)
            # working observations
            yw = eta + (y-mu)/hp
            # least squares rescaled with wroot
            wroot = sparse.dia_matrix((wroot, 0), shape=(n_samples, n_samples))
            X_rescale = safe_sparse_dot(wroot, Xnew, dense_output=True)
            yw_rescale = safe_sparse_dot(wroot, y, dense_output=True)
            coef = linalg.lstsq(X_rescale, yw_rescale)[0]
        elif start_params is 'ols':
            reg = LinearRegression(copy_X=False, fit_intercept=False)
            reg.fit(Xnew, self.link.link(y))
            coef = reg.coef_
        else:
            coef = start_params

        # algorithms for optimiation
        # TODO: Parallelize it
        self.n_iter_ = 0
        converged = False
        if self.solver == 'irls':
            # linear predictor
            eta = safe_sparse_dot(Xnew, coef, dense_output=True)
            mu = self.link.inverse(eta)
            while self.n_iter_ < self.max_iter:
                self.n_iter_ += 1
                # coef_old not used so far.
                # coef_old = coef
                # h'(eta)
                hp = self.link.inverse_derivative(eta)
                # working weights w, in principle a diagonal matrix
                # therefore here just as 1d array
                w = (hp**2 / self.family.variance(mu, phi=1, weight=weight))
                wroot = np.sqrt(w)
                # working observations
                yw = eta + (y-mu)/hp
                # least squares rescaled with wroot
                wroot = sparse.dia_matrix((wroot, 0),
                                          shape=(n_samples, n_samples))
                X_rescale = safe_sparse_dot(wroot, Xnew, dense_output=True)
                yw_rescale = safe_sparse_dot(wroot, yw, dense_output=True)
                coef, residues, rank, singular_ = (
                    linalg.lstsq(X_rescale, yw_rescale))

                # updated linear predictor
                # do it here for updated values for tolerance
                eta = safe_sparse_dot(Xnew, coef, dense_output=True)
                mu = self.link.inverse(eta)

                # which tolerace? |coef - coef_old| or gradient?
                # use gradient for compliance with newton-cg and lbfgs
                # TODO: faster computation of gradient, use mu and eta directly
                gradient = self.family._deviance_derivative(
                    coef=coef, X=Xnew, y=y, weight=weight, link=self.link)
                if (np.max(np.abs(gradient)) <= self.tol):
                    converged = True
                    break

            if not converged:
                warnings.warn("irls failed to converge. Increase the number "
                              "of iterations (currently {0})"
                              .format(self.max_iter), ConvergenceWarning)

        # TODO: performance: make one function return both deviance and
        #       gradient of deviance
        elif self.solver == 'lbfgs':
            func = self.family._deviance
            fprime = self.family._deviance_derivative
            args = (Xnew, y, weight, self.link)
            coef, loss, info = optimize.fmin_l_bfgs_b(
                func, coef, fprime=fprime,
                args=args,
                iprint=(self.verbose > 0) - 1, pgtol=self.tol,
                maxiter=self.max_iter)
            if self.verbose > 0:
                if info["warnflag"] == 1:
                    warnings.warn("lbfgs failed to converge."
                                  " Increase the number of iterations.",
                                  ConvergenceWarning)
                elif info["warnflag"] == 2:
                    warnings.warn("lbfgs failed for the reason: {0}".format(
                        info["task"]))
            self.n_iter_ = info['nit']
        elif self.solver == 'newton-cg':
            func = self.family._deviance
            grad = self.family._deviance_derivative

            def grad_hess(coef, X, y, weight, link):
                grad = (self.family._deviance_derivative(
                    coef, X, y, weight, link))
                hessian = (self.family._deviance_hessian(
                    coef, X, y, weight, link))

                def Hs(s):
                    ret = np.dot(hessian, s)
                    return ret
                return grad, Hs
            hess = grad_hess
            args = (Xnew, y, weight, self.link)
            coef, n_iter_i = newton_cg(hess, func, grad, coef, args=args,
                                       maxiter=self.max_iter, tol=self.tol)
            self.coef_ = coef

        if self.fit_intercept is True:
            self.intercept_ = coef[0]
            self.coef_ = coef[1:]
        else:
            self.coef_ = coef

        if self.fit_dispersion in ['chisqr', 'deviance']:
            self.dispersion_ = self.estimate_phi(y, X, weight)

        return self

    def predict(self, X, weight=1):
        """Prediction with features X.
        If weights are given, returns prediction*weights.
        """
        check_is_fitted(self, "coef_")
        eta = safe_sparse_dot(X, self.coef_, dense_output=True)
        if self.fit_intercept is True:
            eta += self.intercept_
        mu = self.link.inverse(eta)
        return mu*weight

    def estimate_phi(self, y, X, weight):
        n_samples, n_features = X.shape
        eta = safe_sparse_dot(X, self.coef_, dense_output=True)
        if self.fit_intercept is True:
            eta += self.intercept_
        mu = self.link.inverse(eta)
        if self.fit_dispersion == 'chisqr':
            chisq = np.sum(weight*(y-mu)**2/self.family.unit_variance(mu))
            return chisq/(n_samples - n_features)
        elif self.fit_dispersion == 'deviance':
            dev = self.family.deviance(y, mu, weight)
            return dev/(n_samples - n_features)

    def score(self, X, y, weight=1):
        """The natural score for a GLM is -deviance.
        Returns the weight averaged negitive deviance (the better the score,
        the better the fit). Maximum score is therefore 0.
        """
        # RegressorMixin has R^2 score.
        # TODO: Make it more compatible with the score function in
        #      sklearn.metrics.regression.py
        eta = safe_sparse_dot(X, self.coef_, dense_output=True)
        if self.fit_intercept is True:
            eta += self.intercept_
        mu = self.link.inverse(eta)
        output_errors = self.family.unit_deviance(y, mu)
        weight = weight * np.ones_like(y)
        return np.average(output_errors, weights=weight)
