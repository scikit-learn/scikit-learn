# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

r"""Utilities for orthogonal polynomials.

An orthogonal polynomial sequence is a family of polynomials such that any two
different polynomials in the sequence are orthogonal to each other under some
inner product. This file contains utility classes and functions that can be
used to construct orthogonal polynomials. The main goal of the utilities in
this file is to provide a unified framework for evaluating sequences of
(one-dimensional) orthogonal polynomials up to a certain degree in a given set
of points (aka constructing the Vandermonde matrix).

The main (abstract) class is `Polynomial`. It provides the method
`vandermonde`, which evaluates the Vandermonde matrix in a given set of points
:math:`\mathbf x` up to a given degree of the orthogonal polynomial :math:`k`.
The abstract class `Polynomial` also provides the factory method
`from_distribution`, which generates an instance of the orhogonal polynomial
that corresponds to a given `scipy.stats` distribution. Finally, the class
provides the method `scale_features_from_distribution`, which can be used to
scale a given set of input features, distributed according to a given
distribution, to the support of the corresponding orthogonal polynomial.

Currently, this class supports the following orthogonal polynomial types:
- `Hermite` (`scipy.stats.norm`)
- `Jacobi` (`scipy.stats.beta`)
- `Laguerre` (`scipy.stats.expon`)
- `Legendre` (`scipy.stats.uniform`)

New orthogonal polynomials can easily be added by extending the `Polynomial`
base clase, and providing an implementation for the methods `_vandermonde`,
`_distribution` and `_norm_squared`.

Example
-------
>>> import numpy as np
>>> from sklearn.utils._orthogonal_polynomial import Legendre
>>> X = np.linspace(0, 1, num=4)
>>> leg = Legendre()
>>> leg.vandermonde(X, 3)
array([[ 1.        ,  0.        , -0.5       , -0.        ],
       [ 1.        ,  0.33333333, -0.33333333, -0.40740741],
       [ 1.        ,  0.66666667,  0.16666667, -0.25925926],
       [ 1.        ,  1.        ,  1.        ,  1.        ]])
"""

from abc import ABC, abstractmethod
from math import factorial

import numpy as np
from numpy.polynomial.hermite_e import hermevander
from numpy.polynomial.laguerre import lagvander
from numpy.polynomial.legendre import legvander
from scipy.special import eval_jacobi, gamma

from sklearn.utils._param_validation import Integral, Real
from sklearn.utils.validation import check_array, check_scalar, column_or_1d


class Polynomial(ABC):
    """An abstract base class for polynomials.

    Parameters
    ----------
    normalize : bool
        If True, return an orthonormal polynomial.
    """

    def __init__(self, normalize: bool = False):
        check_scalar(normalize, "normalize", bool)
        self.normalize = normalize

    def vandermonde(self, points, degree):
        r"""Returns the generalized Vandermonde matrix associated with this
        orthogonal polynomial.

        For a given set of :math:`n` points :math:`\mathbf x`, this function
        constructs the generalized Vandermonde matrix :math:`V` of dimensions
        :math:`n \times (k + 1)`, where :math:`k` is the degree of this
        polynomial, with entries

        .. math::
            V_{i, j} = \phi_j(x_i)

        Parameters
        ----------
        points : array-like of length (`n_points`)
            The points :math:`\mathbf x` at which to evaluate this orthogonal
            polynomial. The points must be interpretable as a 1d array.

        degree : int
            The maximum degree :math:`k` of this orthogonal polynomial to
            consider. The number of columns in the Vandermonde matrix will be
            actual to :math:`k + 1`.

        Returns
        -------
        V : array-like of shape (`n_points`, `degree + 1`)
            The generalized Vandermonde matrix.
        """
        # Check points
        try:
            points = column_or_1d(points)
        except ValueError:
            raise ValueError("could not interpret points as 1d array")

        # Check degree
        if not isinstance(degree, Integral) or degree < 0:
            raise ValueError(f"degree must be a non-negative int, got '{degree}'")

        # Compute orthogonal Vandermonde
        V = self._vandermonde(points, degree)

        # Normalize columns if requested
        if self.normalize:
            for j in range(degree + 1):
                V[:, j] /= np.sqrt(self._norm_squared(j))

        return V

    @abstractmethod
    def _vandermonde(self, points, degree):
        """This method must be overridden by concrete classes."""

    # This is an alternative construction method that returns an instance of an
    # orthogonal polynomial, given a `scipy.stats` distribution. This method is
    # used, for example, in Polynomial Chaos Expansions: for a given
    # distribution, we want to find the corresponding orthogonal polynomial
    # type. An alternative implementation would be to use an if/else or or
    # match/case to find the corresponding polynomial type. Instead, this
    # function uses `__subclasses__` to automatically loop over all implemented
    # orthogonal polynomial types that inherit from this class, and selects the
    # appropriate polynomial type based on the distribution name (returned by
    # the `_distribution` method). The advantage of this implementation is
    # extensibility: users can define their own distribution and matching
    # orthogonal polynomial (inheriting from this abstract class), and
    # everything else should still work out of the box.
    @classmethod
    def from_distribution(cls, distribution):
        """Return an orthogonal polynomial that corresponds to the given
        distribution.

        Parameters
        ----------
        distribution : scipy.stats frozen distribution
            The distribution for which the corresponding orthogonal polynomial
            should be returned. This distribution should be given as a *frozen*
            distribution from `scipy.stats`, generated using e.g., `uniform()`
            or `norm()`.

        Returns
        -------
        polynomial : Polynomial
            An instance of the corresponding orthogonal polynomial.

        Note
        ----
        The polynomials are orthogonal with respect to a certain weight
        function. When that weight function is equal to the probability
        density function of a random variable, the polynomials can be used
        as a basis to represent this random variable.

        When extending the Polynomial class with a new orthogonal polynomial
        type, this method should be able to automatically detect that new
        polynomial, if the concrete class implements the `_distribution`
        method.
        """
        # Check distribution
        if not hasattr(distribution, "dist"):
            raise ValueError(
                "this distribution does not have a 'dist' attribute and "
                "cannot be interpreted as a 'scipy.stats' frozen distribution"
            )

        # Match subclass
        name = distribution.dist.name
        for polynomial in Polynomial.__subclasses__():
            if name == polynomial._distribution():
                # Find the shape parameters of this distribution ...
                shape, _, _ = distribution.dist._parse_args(
                    *distribution.args, **distribution.kwds
                )
                # ... and pass them as arguments in the constructor
                if len(shape) > 0:
                    return polynomial(shape_params=shape)
                else:
                    return polynomial()
        raise ValueError(f"no polynomial type matches distribution {name}'")

    # Returns the `scipy.stats` distribution name
    @staticmethod
    @abstractmethod
    def _distribution():
        """This method must be overridden by concrete classes."""

    def scale_features_from_distribution(self, X, distribution):
        """Scale the given features, distributed according to the given
        distribution, to the support of this orthogonal polynomial. For
        example, features distributed according to a uniform distribution on
        `[a, b]` are mapped to `[-1, 1]`, since that is the support of the
        Lagrange polynomials, which are orthogonal with respect to the uniform
        distribution on `[-1, 1]`.

        Parameters
        ----------
        X : array-like of shape (`n_samples`, `n_features`)
            The features to scale.

        distribution : scipy.stats frozen distribution
            The distribution of the features.
        """
        # Check distribution
        if not hasattr(distribution, "dist"):
            raise ValueError(
                "this distribution does not have a 'dist' attribute and "
                "cannot be interpreted as a 'scipy.stats' frozen distribution"
            )

        # Check features
        try:
            X = check_array(X, ensure_2d=False)
        except ValueError:
            raise ValueError("could not interpret feature matrix as 2d array")

        # Scale features
        if distribution.dist.name != self._distribution():
            raise ValueError(
                "this orthogonal polynomial can only scale features "
                f"according to '{self._distribution()}' distributions, got "
                f"'{distribution.dist.name}' distribution instead"
            )
        _, loc, scale = distribution.dist._parse_args(
            *distribution.args, **distribution.kwds
        )
        if scale != 0:
            return (X - loc) / scale
        return X

    # Not an abstract static method because polynomials may be parametrized
    # (e.g. Jacobi polynomials)
    def norm(self, degree):
        """Returns the norm of this orthogonal polynomial of the given degree.

        Parameters
        ----------
        degree : int
            The degree of this orthogonal polynomial to consider.

        Returns
        -------
        norm : float
            The norm of the orthogonal polynomial of the given degree.
        """
        # Check degree
        if not isinstance(degree, Integral) or degree < 0:
            raise ValueError(f"degree must be a non-negative int, got '{degree}'")

        # Compute norm
        return 1.0 if self.normalize else np.sqrt(self._norm_squared(degree))

    @abstractmethod
    def _norm_squared(self, degree):
        """This method must be overridden by concrete classes."""

    def __repr__(self):
        return self.__class__.__name__


class Hermite(Polynomial):
    """A class representing Hermite polynomials."""

    def _vandermonde(self, points, degree):
        return hermevander(points, degree)

    @staticmethod
    def _distribution():
        return "norm"

    def _norm_squared(self, degree):
        return factorial(degree)


class Jacobi(Polynomial):
    """A class representing Jacobi polynomials."""

    def __init__(self, alpha=None, beta=None, shape_params=None, normalize=False):
        super().__init__(normalize=normalize)

        # Get alpha and beta from shape parameters, if provided
        if shape_params is not None:
            if alpha is not None or beta is not None:
                raise ValueError(
                    "both 'shape_params' and 'alpha' and 'beta' are specified"
                )
            if len(shape_params) != 2:
                raise ValueError(
                    "need exactly 2 parameters to generate Jacobi polynomial "
                    "from distribution shape parameters, got "
                    f"{len(shape_params)}"
                )
            # For a `Beta(alpha, beta)` distribution, the associated Jacobi
            # polynomial has parameters `beta - 1` and `alpha - 1`
            alpha = shape_params[1] - 1
            beta = shape_params[0] - 1

        # Set alpha
        if not isinstance(alpha, Real):
            raise ValueError(f"alpha must be float or int, got '{type(alpha)}'")
        if not (alpha > 0):
            raise ValueError(f"alpha must be > 0, got '{alpha}'")
        self.alpha = alpha

        # Set beta
        if not isinstance(beta, Real):
            raise ValueError(f"beta must be float or int, got '{type(beta)}'")
        if not (beta > 0):
            raise ValueError(f"beta must be > 0, got '{beta}'")
        self.beta = beta

    def _vandermonde(self, points, degree):
        V = np.zeros((len(points), degree + 1))
        for j in range(degree + 1):
            V[:, j] = eval_jacobi(j, self.alpha, self.beta, points)
        return V

    @staticmethod
    def _distribution():
        return "beta"

    # Jacobi polynomials are associated with the beta distribution,
    # and are supported on [-1, 1]. Since the beta distribution is defined on
    # [0, 1], we additionally need to scale the features from [0, 1] to
    # [-1, 1].
    def scale_features_from_distribution(self, X, distribution):
        X = super().scale_features_from_distribution(X, distribution)
        return 2 * X - 1

    def _norm_squared(self, degree):
        # The alpha and beta here are parameters for the Jacobi polynomial
        term1 = 1.0 / (2 * degree + self.alpha + self.beta + 1)
        term2 = (
            gamma(self.alpha + self.beta + 2)
            / gamma(self.alpha + 1)
            / gamma(self.beta + 1)
        )
        term3 = (
            gamma(degree + self.alpha + 1)
            * gamma(degree + self.beta + 1)
            / factorial(degree)
            / gamma(degree + self.alpha + self.beta + 1)
        )
        return term1 * term2 * term3

    def __repr__(self):
        return f"Jacobi_{self.alpha}_{self.beta}"


class Laguerre(Polynomial):
    """A class representing Laguerre polynomials."""

    def _vandermonde(self, points, degree):
        return lagvander(points, degree)

    @staticmethod
    def _distribution():
        return "expon"

    def _norm_squared(self, degree):
        return 1


class Legendre(Polynomial):
    """A class representing Legendre polynomials."""

    def _vandermonde(self, points, degree):
        return legvander(points, degree)

    @staticmethod
    def _distribution():
        return "uniform"

    # Legendre polynomials are associated with the uniform distribution,
    # and are supported on [-1, 1]. Since the standard uniform distribution
    # is defined on [0, 1], we additionally need to scale the features from
    # [0, 1] to [-1, 1].
    def scale_features_from_distribution(self, X, distribution):
        X = super().scale_features_from_distribution(X, distribution)
        return 2 * X - 1

    def _norm_squared(self, degree):
        return 1 / (2 * degree + 1)
