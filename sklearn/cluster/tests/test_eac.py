"""
Tests for DBSCAN clustering algorithm
"""

import pickle

import numpy as np
from scipy.spatial import distance

from sklearn.utils.testing import assert_equal
from sklearn.cluster import EAC
from .common import generate_clustered_data


n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)


def test_eac():
    """Tests the EAC algorithm runs with basic data."""
    clusterer = EAC().fit(X)
    
