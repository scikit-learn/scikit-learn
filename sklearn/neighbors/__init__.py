"""The k-nearest neighbors algorithms."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.neighbors._ball_tree import BallTree
from sklearn.neighbors._base import (
    VALID_METRICS,
    VALID_METRICS_SPARSE,
    sort_graph_by_row_values,
)
from sklearn.neighbors._classification import (
    KNeighborsClassifier,
    RadiusNeighborsClassifier,
)
from sklearn.neighbors._graph import (
    KNeighborsTransformer,
    RadiusNeighborsTransformer,
    kneighbors_graph,
    radius_neighbors_graph,
)
from sklearn.neighbors._kd_tree import KDTree
from sklearn.neighbors._kde import KernelDensity
from sklearn.neighbors._lof import LocalOutlierFactor
from sklearn.neighbors._nca import NeighborhoodComponentsAnalysis
from sklearn.neighbors._nearest_centroid import NearestCentroid
from sklearn.neighbors._regression import KNeighborsRegressor, RadiusNeighborsRegressor
from sklearn.neighbors._unsupervised import NearestNeighbors

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
