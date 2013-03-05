"""
Testing for Clustering methods

"""

import numpy as np
import warnings

from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_raises)
from sklearn.cluster.affinity_propagation_ import AffinityPropagation
from sklearn.cluster.affinity_propagation_ import affinity_propagation
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import euclidean_distances

n_clusters = 3
centers = np.array([[1, 1], [-1, -1], [1, -1]]) + 10
X, _ = make_blobs(n_samples=60, n_features=2, centers=centers,
                  cluster_std=0.4, shuffle=True, random_state=0)


def test_affinity_propagation():
    """Affinity Propagation algorithm
    """
    # Compute similarities
    S = -euclidean_distances(X, squared=True)
    preference = np.median(S) * 10
    # Compute Affinity Propagation
    cluster_centers_indices, classes = affinity_propagation(
        S, preference=preference)

    n_clusters_ = len(cluster_centers_indices)

    assert_equal(n_clusters, n_clusters_)

    af = AffinityPropagation(preference=preference, affinity="precomputed")
    classes_precomputed = af.fit(S).classes_

    af = AffinityPropagation(preference=preference, verbose=True)
    classes = af.fit(X).classes_

    assert_array_equal(classes, classes_precomputed)

    #Test for deprecations of labels_
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _ = af.fit(X).labels_
        # Verify some things
        assert len(w) == 1
        assert issubclass(w[-1].category, DeprecationWarning)
        assert "deprecated" in str(w[-1].message)

    cluster_centers_indices = af.cluster_centers_indices_

    n_clusters_ = len(cluster_centers_indices)
    assert_equal(np.unique(classes).size, n_clusters_)
    assert_equal(n_clusters, n_clusters_)

    # Test also with no copy
    _, classes_no_copy = affinity_propagation(S, preference=preference,
                                              copy=False)
    assert_array_equal(classes, classes_no_copy)

    # Test input validation
    assert_raises(ValueError, affinity_propagation, S[:, :-1])
    assert_raises(ValueError, affinity_propagation, S, damping=0)
    af = AffinityPropagation(affinity="unknown")
    assert_raises(ValueError, af.fit, X)
