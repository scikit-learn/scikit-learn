"""
Tests for DBSCAN clustering algorithm
"""

import pickle
import warnings

import numpy as np
from scipy.spatial import distance

from sklearn.utils.testing import assert_equal
from sklearn.cluster.dbscan_ import DBSCAN, dbscan
from .common import generate_clustered_data


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
    core_samples, classes = dbscan(D, metric="precomputed", eps=eps,
                                   min_samples=min_samples)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(classes)) - (1 if -1 in classes else 0)

    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric="precomputed", eps=eps, min_samples=min_samples)
    classes = db.fit(D).classes_

    #Test for deprecations of labels_
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _ = db.fit(D).labels_
        # Verify some things
        print w
        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert "deprecated" in str(w[-1].message)

    n_clusters_2 = len(set(classes)) - int(-1 in classes)
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
    core_samples, classes = dbscan(X, metric=metric, eps=eps,
                                   min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(classes)) - int(-1 in classes)
    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric=metric, eps=eps, min_samples=min_samples)
    classes = db.fit(X).classes_

    n_clusters_2 = len(set(classes)) - int(-1 in classes)
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
    core_samples, classes = dbscan(X, metric=metric, eps=eps,
                                   min_samples=min_samples)

    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(classes)) - int(-1 in classes)
    assert_equal(n_clusters_1, n_clusters)

    db = DBSCAN(metric=metric, eps=eps, min_samples=min_samples)
    classes = db.fit(X).classes_

    n_clusters_2 = len(set(classes)) - int(-1 in classes)
    assert_equal(n_clusters_2, n_clusters)


def test_pickle():
    obj = DBSCAN()
    s = pickle.dumps(obj)
    assert_equal(type(pickle.loads(s)), obj.__class__)
