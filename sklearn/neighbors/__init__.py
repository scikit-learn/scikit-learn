"""
The :mod:`sklearn.neighbors` module implements the k-nearest neighbors
algorithm.
"""

__all__ = ['BallTree',
           'kneighbors_graph', 'radius_neighbors_graph',
           'NearestNeighbors',
           'KNeighborsClassifier', 'RadiusNeighborsClassifier',
           'KNeighborsRegressor', 'RadiusNeighborsRegressor']

from .ball_tree import BallTree
from .graph import kneighbors_graph, radius_neighbors_graph
from .unsupervised import NearestNeighbors
from .classification import KNeighborsClassifier, RadiusNeighborsClassifier
from .regression import KNeighborsRegressor, RadiusNeighborsRegressor
from .nearest_centroid import NearestCentroid
