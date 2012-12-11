"""
The :mod:`sklearn.random_projection` module implements several random
projection matrix.
"""

from sklearn.random_projection.random_projection import (
    BernouilliRandomProjection,
    GaussianRandomProjection,
    johnson_lindenstrauss_min_dim)

__all__ = [
    "BernouilliRandomProjection",
    "GaussianRandomProjection"
    "johnson_lindenstrauss_min_dim",
]
