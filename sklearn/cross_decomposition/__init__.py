"""
The :mod:`sklearn.cross_decomposition` module includes algorithms for cross
decomposition techniques.
"""

from ._pls import CCA, PLSSVD, PLSCanonical, PLSRegression

__all__ = ["PLSCanonical", "PLSRegression", "PLSSVD", "CCA"]
