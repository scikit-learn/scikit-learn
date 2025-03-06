"""The k-nearest neighbors algorithms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._ball_tree import BallTree
from ._base import VALID_METRICS, VALID_METRICS_SPARSE, sort_graph_by_row_values
from ._classification import KNeighborsClassifier, RadiusNeighborsClassifier
from ._graph import (
    KNeighborsTransformer,
    RadiusNeighborsTransformer,
    kneighbors_graph,
    radius_neighbors_graph,
)
from ._kd_tree import KDTree
from ._kde import KernelDensity
from ._lof import LocalOutlierFactor
from ._nca import NeighborhoodComponentsAnalysis
from ._nearest_centroid import NearestCentroid
from ._regression import KNeighborsRegressor, RadiusNeighborsRegressor
from ._unsupervised import NearestNeighbors

__all__ = [
    "VALID_METRICS",
    "VALID_METRICS_SPARSE",
    "BallTree",
    "KDTree",
    "KNeighborsClassifier",
    "KNeighborsRegressor",
    "KNeighborsTransformer",
    "KernelDensity",
    "LocalOutlierFactor",
    "NearestCentroid",
    "NearestNeighbors",
    "NeighborhoodComponentsAnalysis",
    "RadiusNeighborsClassifier",
    "RadiusNeighborsRegressor",
    "RadiusNeighborsTransformer",
    "kneighbors_graph",
    "radius_neighbors_graph",
    "sort_graph_by_row_values",
]
