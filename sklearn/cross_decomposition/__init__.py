"""Algorithms for cross decomposition."""

from ._pls import CCA, PLSSVD, PLSCanonical, PLSRegression

__all__ = ["PLSCanonical", "PLSRegression", "PLSSVD", "CCA"]
