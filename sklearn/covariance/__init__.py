"""
The :mod:`sklearn.covariance` module includes methods and algorithms to
robustly estimate the covariance of features given a set of points. The
precision matrix defined as the inverse of the covariance is also estimated.
Covariance estimation is closely related to the theory of Gaussian Graphical
Models.
"""

from .empirical_covariance_ import empirical_covariance, EmpiricalCovariance, \
    log_likelihood
from .shrunk_covariance_ import shrunk_covariance, ShrunkCovariance, \
    ledoit_wolf, ledoit_wolf_shrinkage, \
    LedoitWolf, oas, OAS
from .robust_covariance import fast_mcd, MinCovDet
from .graph_lasso_ import graphical_lasso, GraphicalLasso, GraphicalLassoCV
from .elliptic_envelope import EllipticEnvelope


__all__ = ['EllipticEnvelope',
           'EmpiricalCovariance',
           'GraphicalLasso',
           'GraphicalLassoCV',
           'LedoitWolf',
           'MinCovDet',
           'OAS',
           'ShrunkCovariance',
           'empirical_covariance',
           'fast_mcd',
           'graphical_lasso',
           'ledoit_wolf',
           'ledoit_wolf_shrinkage',
           'log_likelihood',
           'oas',
           'shrunk_covariance']
