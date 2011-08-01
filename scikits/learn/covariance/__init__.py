"""
Covariance estimators
=====================

:mod:`scikits.learn.covariance` is a module to fit to estimate
robustly the covariance of features given a set of points.
The precision matrix defined as the inverse of the covariance
is also estimated. Covariance estimation is closely related
to the theory of Gaussian Graphical Models.

"""

from .empirical_covariance_ import empirical_covariance, EmpiricalCovariance, \
    log_likelihood
from .shrunk_covariance_ import shrunk_covariance, ShrunkCovariance, \
    ledoit_wolf, LedoitWolf, oas, OAS
