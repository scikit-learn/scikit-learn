"""
Distribution functions used in GLM
"""

# Author: Christian Lorentzen <lorentzen.ch@googlemail.com>
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod
from collections import namedtuple
import numbers

import numpy as np
from scipy.special import xlogy


def _safe_lin_pred(X, coef):
    """Compute the linear predictor taking care if intercept is present."""
    if coef.size == X.shape[1] + 1:
        return X @ coef[1:] + coef[0]
    else:
        return X @ coef


DistributionBoundary = namedtuple("DistributionBoundary",
                                  ("value", "inclusive"))


class ExponentialDispersionModel(metaclass=ABCMeta):
    """Base class for reproductive Exponential Dispersion Models (EDM).

    The pdf of Y∼EDM(μ, φ) is given by::

        p(y| θ, φ) = c1(y, φ) * exp((θy-A(θ))/φ)
        = c2(y, φ) * exp(-d(y, μ)/(2φ))

    with mean E[Y] = A'(θ) = μ, variance Var[Y] = φ * v(μ),
    unit variance v(μ), unit deviance d(y,μ) and dispersion parameter φ.

    Methods
    -------
    deviance
    deviance_derivative
    in_y_range
    unit_deviance
    unit_deviance_derivative
    unit_variance
    unit_variance_derivative

    References
    ----------
    https://en.wikipedia.org/wiki/Exponential_dispersion_model.
    """

    def in_y_range(self, y):
        """Returns ``True`` if y is in the valid range of Y∼EDM.

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.
        """
        if hasattr(self, '_upper_bound'):
            # All currently supported distributions have an upper bound at
            # +inf, however this may need to be implemented for other
            # distributions
            raise NotImplementedError

        if not isinstance(self._lower_bound, DistributionBoundary):
            raise TypeError('_lower_bound attribute must be of type '
                            'DistributionBoundary')

        if self._lower_bound.inclusive:
            return np.greater_equal(y, self._lower_bound.value)
        else:
            return np.greater(y, self._lower_bound.value)

    @abstractmethod
    def unit_variance(self, mu):
        """Compute the unit variance function.

        The unit variance v(μ) determines the variance as a function of the
        mean μ by Var[Y_i] = φ/s_i * v(μ_i).
        It can also be derived from the unit deviance d(y,μ) as::

            v(μ) = 2/(∂^2 d(y,μ)/(∂ μ^2))|_{y=μ}

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.
        """
        pass  # pragma: no cover

    @abstractmethod
    def unit_variance_derivative(self, mu):
        """Compute the derivative of the unit variance w.r.t. mu.

        Return v'(μ).

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Target values.
        """
        pass  # pragma: no cover

    @abstractmethod
    def unit_deviance(self, y, mu, check_input=False):
        """Compute the unit deviance.

        The unit_deviance d(y,μ) can be defined by the log-likelihood as::

        d(y,μ) = -2φ * (loglike(y,μ,φ) - loglike(y,y,φ))

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.

        check_input : bool, default=False
            If True raise an exception on invalid y or mu values, otherwise
            they will be propagated as NaN.
        Returns
        -------
        deviance: array, shape (n_samples,)
            Computed deviance
        """
        pass  # pragma: no cover

    def unit_deviance_derivative(self, y, mu):
        """Compute the derivative of the unit deviance w.r.t. mu.

        The derivative of the unit deviance is given by
        ∂ d(y,μ)/(∂ μ) = -2(y-μ)/v(μ) with unit variance v(μ).

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.
        """
        return -2 * (y - mu) / self.unit_variance(mu)

    def deviance(self, y, mu, weights=1):
        """Compute the deviance.

        The deviance is a weighted sum of the per sample unit deviances,
        D = sum_i s_i * d(y_i,μ_i)
        with weights s_i and unit deviance d(y,μ).
        In terms of the log-likelihood it is
        D = -2φ * (loglike(y,μ,φ/s) - loglike(y,y,φ/s)).

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

        It gives ∂ D(y, μ; weights)/(∂ μ).

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

    def _mu_deviance_derivative(self, coef, X, y, weights, link):
        """Compute mu and the derivative of the deviance w.r.t coef."""
        lin_pred = _safe_lin_pred(X, coef)
        mu = link.inverse(lin_pred)
        d1 = link.inverse_derivative(lin_pred)
        temp = d1 * self.deviance_derivative(y, mu, weights)
        if coef.size == X.shape[1] + 1:
            devp = np.concatenate(([temp.sum()], temp @ X))
        else:
            devp = temp @ X  # same as X.T @ temp
        return mu, devp


class TweedieDistribution(ExponentialDispersionModel):
    """A class for the Tweedie distribution.

    A Tweedie distribution with mean μ=E[Y] is uniquely defined by it's
    mean-variance relationship Var[Y] ∝ μ^power.

    Special cases are:

    ===== ================
    Power Distribution
    ===== ================
    0     Normal
    1     Poisson
    (1,2) Compound Poisson
    2     Gamma
    3     Inverse Gaussian

    Parameters
    ----------
    power : float (default=0)
            The variance power of the unit variance v(μ) = μ^power.
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

        if power <= 0:
            # Extreme Stable or Normal distribution
            self._lower_bound = DistributionBoundary(-np.Inf, inclusive=False)
        elif 0 < power < 1:
            raise ValueError('Tweedie distribution is only defined for '
                             'power<=0 and p>=1.')
        elif 1 <= power < 2:
            # Poisson or Compound Poisson distribution
            self._lower_bound = DistributionBoundary(0, inclusive=True)
        elif power >= 2:
            # Gamma, Positive Stable, Inverse Gaussian distributions
            self._lower_bound = DistributionBoundary(0, inclusive=False)
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
        distribution v(mu)=power * mu**(power-1).

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Predicted mean.
        """
        return self.power * np.power(mu, self.power - 1)

    def unit_deviance(self, y, mu, check_input=False):
        """Compute the unit deviance.

        The unit deviance d(y,μ) can be defined by the log-likelihood as
        d(y,μ) = -2φ * (loglike(y,μ,φ) - loglike(y,y,φ)).

        Parameters
        ----------
        y : array, shape (n_samples,)
            Target values.

        mu : array, shape (n_samples,)
            Predicted mean.

        check_input : bool, default=False
            If True raise an exception on invalid y or mu values, otherwise
            they will be propagated as NaN.
        Returns
        -------
        deviance: array, shape (n_samples,)
            Computed deviance
        """
        p = self.power

        if check_input:
            message = ("Mean Tweedie deviance error with p={} can only be "
                       "used on ".format(p))
            if p < 0:
                # 'Extreme stable', y any realy number, mu > 0
                if (mu <= 0).any():
                    raise ValueError(message + "strictly positive mu.")
            elif p == 0:
                # Normal, y and mu can be any real number
                pass
            elif 0 < p < 1:
                raise ValueError("Tweedie deviance is only defined for p<=0 "
                                 "and p>=1.")
            elif 1 <= p < 2:
                # Poisson and Compound poisson distribution, y >= 0, mu > 0
                if (y < 0).any() or (mu <= 0).any():
                    raise ValueError(message + "non-negative y and strictly "
                                     "positive mu.")
            elif p >= 2:
                # Gamma and Extreme stable distribution, y and mu > 0
                if (y <= 0).any() or (mu <= 0).any():
                    raise ValueError(message + "strictly positive y and mu.")
            else:  # pragma: nocover
                # Unreachable statement
                raise ValueError

        if p < 0:
            # 'Extreme stable', y any realy number, mu > 0
            dev = 2 * (np.power(np.maximum(y, 0), 2-p)/((1-p) * (2-p)) -
                       y * np.power(mu, 1-p)/(1-p) +
                       np.power(mu, 2-p)/(2-p))

        elif p == 0:
            # Normal distribution, y and mu any real number
            dev = (y - mu)**2
        elif p < 1:
            raise ValueError("Tweedie deviance is only defined for p<=0 and "
                             "p>=1.")
        elif p == 1:
            # Poisson distribution
            dev = 2 * (xlogy(y, y/mu) - y + mu)
        elif p == 2:
            # Gamma distribution
            dev = 2 * (np.log(mu/y) + y/mu - 1)
        else:
            dev = 2 * (np.power(y, 2-p)/((1-p) * (2-p)) -
                       y * np.power(mu, 1-p)/(1-p) +
                       np.power(mu, 2-p)/(2-p))
        return dev


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


EDM_DISTRIBUTIONS = {
    'normal': NormalDistribution,
    'poisson': PoissonDistribution,
    'gamma': GammaDistribution,
    'inverse-gaussian': InverseGaussianDistribution,
}
