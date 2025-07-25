"""Mixture modeling algorithms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._bayesian_mixture import BayesianGaussianMixture
from ._gaussian_mixture import GaussianMixture

__all__ = ["BayesianGaussianMixture", "GaussianMixture"]
