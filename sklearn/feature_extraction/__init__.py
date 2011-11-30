"""
The :mod:`sklearn.feature_extraction` module deals with feature extraction from
raw data. It currently methods to extract features from text and images.
"""

from .image import img_to_graph, grid_to_graph
from . import text
