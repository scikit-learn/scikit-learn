"""
The :mod:`sklearn.mixture` module implements mixture modeling algorithms.
"""

from .gaussian_mixture import GaussianMixture
from .bayesian_mixture import BayesianGaussianMixture
from .bernoulli_mixture import BernoulliMixture

__all__ = ['GaussianMixture',
           'BayesianGaussianMixture',
           'BernoulliMixture']
