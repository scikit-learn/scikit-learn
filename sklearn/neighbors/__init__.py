"""
The :mod:`sklearn.neighbors` module implements the k-nearest neighbors
algorithm.
"""

from .ball_tree import BallTree
from .kd_tree import KDTree
from .dist_metrics import DistanceMetric
from .graph import kneighbors_graph, radius_neighbors_graph
from .graph import KNeighborsTransformer, RadiusNeighborsTransformer
from .unsupervised import NearestNeighbors
from .classification import KNeighborsClassifier, RadiusNeighborsClassifier
from .regression import KNeighborsRegressor, RadiusNeighborsRegressor
from .nearest_centroid import NearestCentroid
from .kde import KernelDensity
from .lof import LocalOutlierFactor
from .nca import NeighborhoodComponentsAnalysis
from .base import VALID_METRICS, VALID_METRICS_SPARSE

__all__ = ['BallTree',
           'DistanceMetric',
           'KDTree',
           'KNeighborsClassifier',
           'KNeighborsRegressor',
           'KNeighborsTransformer',
           'NearestCentroid',
           'NearestNeighbors',
           'RadiusNeighborsClassifier',
           'RadiusNeighborsRegressor',
           'RadiusNeighborsTransformer',
           'kneighbors_graph',
           'radius_neighbors_graph',
           'KernelDensity',
           'LocalOutlierFactor',
           'NeighborhoodComponentsAnalysis',
           'VALID_METRICS',
           'VALID_METRICS_SPARSE']
