"""Gaussian process based regression and classification."""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# License: BSD 3 clause

from . import kernels
from ._gpc import GaussianProcessClassifier
from ._gpr import GaussianProcessRegressor

__all__ = ["GaussianProcessRegressor", "GaussianProcessClassifier", "kernels"]
