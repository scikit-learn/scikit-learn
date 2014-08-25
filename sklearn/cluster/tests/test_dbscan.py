"""
Tests for DBSCAN clustering algorithm
"""

import pickle

import numpy as np
from numpy.testing import assert_raises

from scipy.spatial import distance

from sklearn.utils.testing import assert_equal
from sklearn.cluster.dbscan_ import DBSCAN, dbscan
from .common import generate_clustered_data
from sklearn.metrics.pairwise import pairwise_distances


n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)


def test_dbscan_similarity():
    """Tests the DBSCAN algorithm with a similarity array."""
    # Parameters chosen specifically for this task.
    eps = 0.15
    min_samples = 10
    # Compute similarities
    D = distance.squareform(distance.pdist(X))
    D /= np.max(D)
    # Compute DBSCAN
    core_samples, labels = dbscan(D, metric="precomputed", eps=eps,
                                  min_samples=min_samples)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - (1 if -1 in labels else 0)

    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric="precomputed", eps=eps, min_samples=min_samples)
    labels = db.fit(D).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)


def test_dbscan_feature():
    """Tests the DBSCAN algorithm with a feature vector array."""
    # Parameters chosen specifically for this task.
    # Different eps to other test, because distance is not normalised.
    eps = 0.8
    min_samples = 10
    metric = 'euclidean'
    # Compute DBSCAN
    # parameters chosen for task
    core_samples, labels = dbscan(X, metric=metric, eps=eps,
                                  min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric=metric, eps=eps, min_samples=min_samples)
    labels = db.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)


def test_dbscan_callable():
    """Tests the DBSCAN algorithm with a callable metric."""
    # Parameters chosen specifically for this task.
    # Different eps to other test, because distance is not normalised.
    eps = 0.8
    min_samples = 10
    # metric is the function reference, not the string key.
    metric = distance.euclidean
    # Compute DBSCAN
    # parameters chosen for task
    core_samples, labels = dbscan(X, metric=metric, eps=eps,
                                  min_samples=min_samples,
                                  algorithm='ball_tree')

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric=metric, eps=eps, min_samples=min_samples,
                algorithm='ball_tree')
    labels = db.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)


def test_dbscan_balltree():
    """Tests the DBSCAN algorithm with balltree for neighbor calculation."""
    eps = 0.8
    min_samples = 10

    D = pairwise_distances(X)
    core_samples, labels = dbscan(D, metric="precomputed", eps=eps,
                                  min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(p=2.0, eps=eps, min_samples=min_samples, algorithm='ball_tree')
    labels = db.fit(X).labels_

    n_clusters_2 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_2, n_clusters)

    db = DBSCAN(p=2.0, eps=eps, min_samples=min_samples, algorithm='kd_tree')
    labels = db.fit(X).labels_

    n_clusters_3 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_3, n_clusters)

    db = DBSCAN(p=1.0, eps=eps, min_samples=min_samples, algorithm='ball_tree')
    labels = db.fit(X).labels_

    n_clusters_4 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_4, n_clusters)

    db = DBSCAN(leaf_size=20, eps=eps, min_samples=min_samples,
                algorithm='ball_tree')
    labels = db.fit(X).labels_

    n_clusters_5 = len(set(labels)) - int(-1 in labels)
    assert_equal(n_clusters_5, n_clusters)


def test_input_validation():
    """DBSCAN.fit should accept a list of lists."""
    X = [[1., 2.], [3., 4.]]
    DBSCAN().fit(X)             # must not raise exception


def test_dbscan_badargs():
    """Test bad argument values: these should all raise ValueErrors"""
    assert_raises(ValueError,
                  dbscan,
                  X, eps=-1.0)
    assert_raises(ValueError,
                  dbscan,
                  X, algorithm='blah')
    assert_raises(ValueError,
                  dbscan,
                  X, metric='blah')
    assert_raises(ValueError,
                  dbscan,
                  X, leaf_size=-1)
    assert_raises(ValueError,
                  dbscan,
                  X, p=-1)


def test_pickle():
    obj = DBSCAN()
    s = pickle.dumps(obj)
    assert_equal(type(pickle.loads(s)), obj.__class__)
