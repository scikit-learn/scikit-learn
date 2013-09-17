"""
Tests for SingleLinkageCluster clustering algorithm
"""

import pickle

import numpy as np
from scipy import sparse
from collections import defaultdict

from sklearn.utils.testing import assert_array_equal
from sklearn.cluster import SingleLinkageCluster
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import pairwise_distances


def test_single_linkage_cluster_cluster():
    """
    Tests that SingleLinkageCluster has same results for different input.
    """
    n_samples = 150
    threshold = 2.0
    centers = np.array([
        [0.0, 5.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 4.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 5.0, 1.0],
    ])
    X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                                cluster_std=1., random_state=42)
    D = pairwise_distances(X, metric='euclidean')
    print np.histogram(D)
    # Using precomputed distance matrix
    model_p = SingleLinkageCluster(threshold=threshold, metric='precomputed')
    labels_p = model_p.fit(D).labels_
    # Computing distance matrix
    model_d = SingleLinkageCluster(threshold=threshold)
    labels_d = model_d.fit(X).labels_
    assert_array_equal(labels_p, labels_d)
    # Using a sparse distance matrix
    Dsp = sparse.csr_matrix(D)
    model_sp = SingleLinkageCluster(threshold=threshold, metric='precomputed')
    labels_sp = model_sp.fit(Dsp).labels_
    assert_array_equal(labels_p, labels_sp)


def test_trivial_mst():
    """Tests a trivial SingleLinkageCluster example (from scipy docs)."""
    X = sparse.csr_matrix([[0, 8, 0, 3],
                           [0, 0, 2, 5],
                           [0, 0, 0, 6],
                           [0, 0, 0, 0]])
    model = SingleLinkageCluster(threshold=4.9, metric='precomputed')
    y_pred = model.fit(X).labels_
    y_true = np.array([0, 1, 1, 0], dtype='int32')
    assert_array_equal(y_pred, y_true)


test_trivial_mst()
