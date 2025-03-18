"""Algorithms for cross decomposition."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._pls import CCA, PLSSVD, PLSCanonical, PLSRegression

__all__ = ["CCA", "PLSSVD", "PLSCanonical", "PLSRegression"]


def __dir__():
    return __all__
