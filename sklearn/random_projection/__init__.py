"""
The :mod:`sklearn.random_projection` module implements several random
projection operators.
"""

from sklearn.random_projection.random_projection import (
    BernoulliRandomProjection,
    GaussianRandomProjection,
    johnson_lindenstrauss_min_dim)

__all__ = [
    "BernoulliRandomProjection",
    "GaussianRandomProjection"
    "johnson_lindenstrauss_min_dim",
]
