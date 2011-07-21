"""
Tests for DBSCAN clustering algorithm
"""

import numpy as np
from numpy.testing import assert_equal
from scipy.spatial import distance

from ..dbscan_ import DBSCAN, dbscan
from .common import generate_clustered_data


n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)


def test_dbscan_similarity():
    """
    Tests the DBSCAN algorithm with a similarity array

    """
    # Compute similarities
    D = distance.squareform(distance.pdist(X))
    S = 1 - (D / np.max(D))

    # Compute DBSCAN
    # parameters chosen for task
    core_points, labels = dbscan(S, metric="precomputed",
                                 eps=0.85, min_points=10)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - (1 if -1 in labels else 0)

    assert_equal(n_clusters, n_clusters_1)

    db = DBSCAN(metric="precomputed")
    labels = db.fit(S, eps=0.85, min_points=10).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters, n_clusters_2)


def test_dbscan_feature():
    """
    Tests the DBSCAN algorithm with a feature vector array

    """
    metric = 'euclidean'
    # Compute DBSCAN
    # parameters chosen for task
    core_points, labels = dbscan(X, metric=metric,
                                 eps=0.85, min_points=10)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters, n_clusters_1)

    db = DBSCAN()
    labels = db.fit(X, metric=metric,
                    eps=0.85, min_points=10).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters, n_clusters_2)
