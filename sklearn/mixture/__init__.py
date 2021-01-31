"""
The :mod:`sklearn.mixture` module implements mixture modeling algorithms.
"""

from ._gaussian_mixture import GaussianMixture
from ._bayesian_mixture import BayesianGaussianMixture
from ._phase_constrained_gaussian_mixture import PhaseConstrainedGaussianMixture
from ._strict_bayesian_mixture import StrictBayesianGaussianMixture


__all__ = ['GaussianMixture',
           'BayesianGaussianMixture',
           'PhaseConstrainedGaussianMixture',
           'StrictBayesianGaussianMixture']
