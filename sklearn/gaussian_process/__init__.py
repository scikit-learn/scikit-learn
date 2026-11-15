"""Gaussian process based regression and classification."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.gaussian_process import kernels
from sklearn.gaussian_process._gpc import GaussianProcessClassifier
from sklearn.gaussian_process._gpr import GaussianProcessRegressor

__all__ = ["GaussianProcessClassifier", "GaussianProcessRegressor", "kernels"]
