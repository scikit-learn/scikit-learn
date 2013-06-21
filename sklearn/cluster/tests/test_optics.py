"""
Tests for OPTICS clustering algorithm
"""

import pickle

import numpy as np
from scipy.spatial import distance

from sklearn.utils.testing import assert_equal
from sklearn.cluster.optics_ import OPTICS, optics
from .common import generate_clustered_data

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)

def test_optics_similarity():
    """Tests the OPTICS algorithm with a similarity array."""
    # Parameters chosen specifically for this task.
    eps = 0.6
    # Compute similarities
    D = distance.squareform(distance.pdist(X))
    # Compute OPTICS
    _, _, _, labels = optics(D, metric="precomputed", eps=eps)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - (1 if -1 in labels else 0)

    assert_equal(n_clusters_1, n_clusters)

    opt = OPTICS(metric="precomputed", eps=eps)
    labels = opt.fit(D).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)

def test_optics_feature():
    """Tests the OPTICS algorithm with a feature vector array."""
    # Parameters chosen specifically for this task.
    eps = 0.6
    metric = 'euclidean'
    # Compute OPTICS
    _, _, _, labels = optics(X, metric=metric, eps=eps)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    opt = OPTICS(metric=metric, eps=eps)
    labels = opt.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)

#def test_optics_callable():
#    """Tests the OPTICS algorithm with a callable metric."""
#    # Parameters chosen specifically for this task.
#    eps = 0.6
#    # metric is the function reference, not the string key.
#    metric = distance.euclidean
#    # Compute OPTICS
#    _, _, _, labels = optics(X, metric=metric, eps=eps)
#
#    # number of clusters, ignoring noise if present
#    n_clusters_1 = len(set(labels)) - int(-1 in labels)
#    assert_equal(n_clusters_1, n_clusters)
#
#    opt = OPTICS(metric=metric, eps=eps)
#    labels = opt.fit(X).labels_
#
#    n_clusters_2 = len(set(labels)) - int(-1 in labels)
#    assert_equal(n_clusters_2, n_clusters)

def test_pickle():
    obj = OPTICS()
    s = pickle.dumps(obj)
    assert_equal(type(pickle.loads(s)), obj.__class__)

