"""Methods and algorithms to robustly estimate covariance.

They estimate the covariance of features at given sets of points, as well as the
precision matrix defined as the inverse of the covariance. Covariance estimation is
closely related to the theory of Gaussian graphical models.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._elliptic_envelope import EllipticEnvelope
from ._empirical_covariance import (
    EmpiricalCovariance,
    empirical_covariance,
    log_likelihood,
)
from ._graph_lasso import GraphicalLasso, GraphicalLassoCV, graphical_lasso
from ._robust_covariance import MinCovDet, fast_mcd
from ._shrunk_covariance import (
    OAS,
    LedoitWolf,
    ShrunkCovariance,
    ledoit_wolf,
    ledoit_wolf_shrinkage,
    oas,
    shrunk_covariance,
)

__all__ = [
    "OAS",
    "EllipticEnvelope",
    "EmpiricalCovariance",
    "GraphicalLasso",
    "GraphicalLassoCV",
    "LedoitWolf",
    "MinCovDet",
    "ShrunkCovariance",
    "empirical_covariance",
    "fast_mcd",
    "graphical_lasso",
    "ledoit_wolf",
    "ledoit_wolf_shrinkage",
    "log_likelihood",
    "oas",
    "shrunk_covariance",
]
