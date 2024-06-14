"""Gaussian process based regression and classification."""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# SPDX-License-Identifier: BSD-3-Clause

from . import kernels
from ._gpc import GaussianProcessClassifier
from ._gpr import GaussianProcessRegressor

__all__ = ["GaussianProcessRegressor", "GaussianProcessClassifier", "kernels"]
