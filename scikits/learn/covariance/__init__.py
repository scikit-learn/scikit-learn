"""
Covariance estimators
=====================

scikits.learn.covariance is a module to fit to estimate
robustly the covariance of features given a set of points.
The precision matrix defined as the inverse of the covariance
is also estimated. Covariance estimation is closely related
to the theory of Gaussian Graphical Models.

"""

from .covariance import Covariance
from .shrunk_covariance import shrunk_covariance, ShrunkCovariance
from .ledoit_wolf import ledoit_wolf, LedoitWolf
