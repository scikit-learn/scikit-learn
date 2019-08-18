"""
Link functions used in GLM
"""

# Author: Christian Lorentzen <lorentzen.ch@googlemail.com>
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.special import expit, logit


class Link(metaclass=ABCMeta):
    """Abstract base class for Link functions."""

    @abstractmethod
    def __call__(self, mu):
        """Compute the link function g(mu).

        The link function links the mean mu=E[Y] to the so called linear
        predictor (X*w), i.e. g(mu) = linear predictor.

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Usually the (predicted) mean.
        """
        pass  # pragma: no cover

    @abstractmethod
    def derivative(self, mu):
        """Compute the derivative of the link g'(mu).

        Parameters
        ----------
        mu : array, shape (n_samples,)
            Usually the (predicted) mean.
        """
        pass  # pragma: no cover

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
        pass  # pragma: no cover

    @abstractmethod
    def inverse_derivative(self, lin_pred):
        """Compute the derivative of the inverse link function h'(lin_pred).

        Parameters
        ----------
        lin_pred : array, shape (n_samples,)
            Usually the (fitted) linear predictor.
        """
        pass  # pragma: no cover

    @abstractmethod
    def inverse_derivative2(self, lin_pred):
        """Compute 2nd derivative of the inverse link function h''(lin_pred).

        Parameters
        ----------
        lin_pred : array, shape (n_samples,)
            Usually the (fitted) linear predictor.
        """
        pass  # pragma: no cover


class IdentityLink(Link):
    """The identity link function g(x)=x."""

    def __call__(self, mu):
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

    def __call__(self, mu):
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

    def __call__(self, mu):
        return logit(mu)

    def derivative(self, mu):
        return 1. / (mu * (1 - mu))

    def inverse(self, lin_pred):
        return expit(lin_pred)

    def inverse_derivative(self, lin_pred):
        ep = expit(lin_pred)
        return ep * (1. - ep)

    def inverse_derivative2(self, lin_pred):
        ep = expit(lin_pred)
        return ep * (1. - ep) * (1. - 2 * ep)
