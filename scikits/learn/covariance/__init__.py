"""
Covariance estimators
=====================

scikits.learn.covariance is a module to fit to estimate
robustly the covariance of features given a set of points.
The precision matrix defined as the inverse of the covariance
is also estimated. Covariance estimation is closely related
to the theory of Gaussian Graphical Models.

"""

from .base_covariance_ import base_covariance, BaseCovariance, log_likelihood
from .shrunk_covariance_ import shrunk_covariance, ShrunkCovariance
from .shrunk_covariance_ import ledoit_wolf, LedoitWolf
from .shrunk_covariance_ import OAS, oas
