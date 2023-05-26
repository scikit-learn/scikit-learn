"""
The :mod:`sklearn.mixture` module implements mixture modeling algorithms.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_bayesian_mixture": ["BayesianGaussianMixture"],
        "_gaussian_mixture": ["GaussianMixture"],
    },
)
