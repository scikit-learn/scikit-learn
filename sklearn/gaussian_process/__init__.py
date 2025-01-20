"""Gaussian process based regression and classification."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from . import kernels
from ._gpc import GaussianProcessClassifier
from ._gpr import GaussianProcessRegressor
from ._tpr import TProcessRegressor

__all__ = [
    "GaussianProcessRegressor",
    "GaussianProcessClassifier",
    "TProcessRegressor",
    "kernels",
]
