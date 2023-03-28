"""
The :mod:`sklearn.mixture` module implements mixture modeling algorithms.
"""

from ._gaussian_mixture import GaussianMixture
from ._bayesian_mixture import BayesianGaussianMixture
from ._base import get_responsibilities


__all__ = ["GaussianMixture", "BayesianGaussianMixture", "get_responsibilities"]
