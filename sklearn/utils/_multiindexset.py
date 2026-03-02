# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Utilities for constructing multiindex sets.

A multiindex is a tuple `(i, j, k, ...)` where the integers `i`, `j` and `k`,
all non-negative, represent the generalization of the concept of an integer
index to an ordered tuple of indices. A multiindex set is a collection of
multiindices. This file contains utility classes and functions that can be used
to construct such multiindex set. A multiindex set of a given type can be
constructed from two required inputs: :math:`d`, the dimension of the
multiindex set, and :math:`k`, the degree of the multiindex set. The degree
governs the size of the multiindex set. Optionally, a *weighted* multiindex set
can be constructed, which allows one to select certain preferential directions
in the space of multiindices.

The main (abstract) class is `MultiIndexSet`. This class provides a general
implementation of a method to generate all indices of the multiindex set of the
desired type. The type of the index set is determined by the concrete class
that inherits from the abstract multiindex set class. These concrete classes
only need to implement the (private) method `_contains(self, index)`, that
determines if a given multiindex belongs to the desired multiindex set.

Currently, the following multiindex set shapes are implemented (illustrated for
2 dimensions and degree 4)

- `FullTensor`
  x x x x x
  x x x x x
  x x x x x
  x x x x x
  x x x x x

- `TotalDegree`
  x
  x x
  x x x
  x x x x
  x x x x x

- `HyperbolicCross`
  x
  x
  x
  x x
  x x x x x

- `ZarembaCross`
  x x
  x x
  x x x
  x x x x x
  x x x x x

New multiindexset types can easily be added by extending the `MultiIndexSet`
base clase, and providing an implementation for the `_contains` method.

Example
-------
>> from sklearn.utils._multiindices import TotalDegree
>> m = TotalDegree(2, 3) # dimension 2, degree 3
>> indices = list(m.indices())
>> indices
[[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [0, 3]
]
"""

# qualified import statements
from abc import ABC, abstractmethod  # for abstract classes
from math import prod  # make prod peer to sum
from re import sub

# sklearn imports
from sklearn.utils._param_validation import Integral, Iterable, Real


class MultiIndexSet(ABC):
    """An abstract base class for multiindex sets.

    Parameters
    ----------
    dimension : int
        The dimension :math:`d` of this multiindex set.

    degree : int
        The degree :math:`k` of this multiindex set.

    weights : array-like of length (`dimension`)
        Optional weights :math:`w_j, 1 < j <d` associated with this
        multiindex set.
    """

    def __init__(self, dimension, degree, weights=None):
        # check dimension
        if not isinstance(dimension, Integral) or dimension < 1:
            raise ValueError(f"dimension must be a non-negative int, got '{dimension}'")
        self.dimension = dimension

        # check degree
        if not isinstance(degree, Integral) or degree < 0:
            raise ValueError(f"degree must be a non-negative int, got '{degree}'")
        self.degree = degree

        # check weights
        if weights is None:
            self.weights = [1] * self.dimension
        elif isinstance(weights, Iterable):
            weights = list(weights)  # weights can be a generator
            for j, weight in enumerate(weights):
                if not isinstance(weight, Real):
                    raise ValueError(
                        "weights must be numbers, but the weight at position "
                        f"{j} has type '{type(weight).__name__}'"
                    )
                if not (weight > 0):
                    raise ValueError(
                        "weights must be > 0, but the weight at position "
                        f"{j} has value '{weight}'"
                    )
            if len(weights) != self.dimension:
                raise ValueError(
                    "number of weights must be equal to the dimension of the "
                    f"multiindex set, expected {self.dimension}, got "
                    f"{len(weights)}"
                )
            self.weights = weights
        else:
            raise ValueError(
                "could not interpret weights of type "
                f"'{type(weights).__name__}' as valid weights"
            )

    def indices(self):
        """Returns all indices in this multiindex set.

        Returns
        -------
        indices : generator
            A generator that supplies the multiindices in this multiindex set.
        """
        j = 0  # parameter that loops over dimensions
        index = [0 for _ in range(self.dimension)]  # [0, 0, ..., 0]
        yield index.copy()  # don't forget to return the zero index
        while True:
            index[j] += 1
            if not self._contains(index):  # index not a part of the index set
                if j == self.dimension - 1:  # all dimensions have been visited
                    break
                index[j] = 0
                j += 1  # move to the next dimension
            else:
                j = 0  # restart with first index
                yield index.copy()

    @abstractmethod
    def _contains(self, index):
        """This method must be overridden by concrete classes."""

    # This is an alternative construction method that returns a multiindex set
    # based on a given name (`str` argument). It uses the `__subclasses__`
    # method to find the matching multiindex set. This method is used, for
    # example, in the `OrthogonalPolynomialFeatures` preprocessing class. An
    # alternative implementation would be to use an if/else or match/case
    # statement to generate the appropriate multiindex set based on a given
    # name. The advantage of this implementation is extensibility: users can
    # define their own multiindex set (inheriting from this abstract class),
    # and everything else should still work out of the box.
    @staticmethod
    def from_string(name):
        """Returns the multiindex set with the given name, if it exists.

        Parameters
        ----------
        name : str
            The name of the multiindex set (lower case and using underscore
            separators).

        Returns
        -------
        multiindex_set : MultiIndexSet
            The multiindex set with the given name. If no multiindex set with
            the given name can be found, a `ValueError` is raised.
        """
        name_ = "".join(word.capitalize() for word in name.split("_"))
        for multiindex_set in MultiIndexSet.__subclasses__():
            if multiindex_set.__name__ == name_:
                return multiindex_set
        raise ValueError(f"unknown multiindex set type '{name}'")

    def __repr__(self):
        out = sub(r"\B([A-Z])", r" \1", self.__class__.__name__)
        out = f"<{out} index set of degree {self.degree} in "
        out += f"{self.dimension} dimensions"
        if any(weight != 1 for weight in self.weights):
            out += f" with weights {self.weights}"
        out += ">"
        return out


class FullTensor(MultiIndexSet):
    r"""A full tensor multiindex set.

    The indices :math:`\boldsymbol{\ell} = \{\ell_j\}_{j=1}^d` in this index
    set satisfy

    .. math:
        \frac{\ell_j}{w_j} < k

    for :math:`1 < j < d`, where :math:`w_j` is the :math:`j`th weight and
    :math:`k` is the degree of the index set.
    """

    def _contains(self, index):
        return all(
            idx / weight <= self.degree for idx, weight in zip(index, self.weights)
        )


class TotalDegree(MultiIndexSet):
    r"""A total degree multiindex set.

    The indices :math:`\boldsymbol{\ell} = \{\ell_j\}_{j=1}^d` in this index
    set satisfy

    .. math:
        \sum_{j} \frac{\ell_j}{w_j} < k

    where :math:`w_j` is the :math:`j`th weight and :math:`k` is the degree of
    the index set.
    """

    def _contains(self, index):
        return (
            sum(idx / weight for idx, weight in zip(index, self.weights)) <= self.degree
        )


class HyperbolicCross(MultiIndexSet):
    r"""A hyperbolic cross multiindex set.

    The indices :math:`\boldsymbol{\ell} = \{\ell_j\}_{j=1}^d` in this index
    set satisfy

    .. math:
        \prod_{j} \left(\frac{\ell_j}{w_j} + 1\right) - 1 < k

    where :math:`w_j` is the :math:`j`th weight and :math:`k` is the degree of
    the index set.
    """

    def _contains(self, index):
        return (
            prod(idx / weight + 1 for idx, weight in zip(index, self.weights))
            <= self.degree + 1
        )


class ZarembaCross(MultiIndexSet):
    r"""A Zaremba cross multiindex set.

    The indices :math:`\boldsymbol{\ell} = \{\ell_j\}_{j=1}^d` in this index
    set satisfy

    .. math:
        \prod_{j} \max\left(\frac{\ell_j}{w_j}, 1\right) < k

    where :math:`w_j` is the :math:`j`th weight and :math:`k` is the degree of
    the index set.
    """

    def _contains(self, index):
        if not self.degree:
            return False
        else:
            return (
                prod(max(1, i / weight) for i, weight in zip(index, self.weights))
                <= self.degree
            )
