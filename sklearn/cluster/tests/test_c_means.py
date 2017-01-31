"""Testing for C-means"""
import warnings

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import (assert_array_almost_equal,
                                   assert_array_less)

from sklearn.cluster import CMeans
from sklearn.cluster.c_means_ import (_cmeans_single_probabilistic,
                                      _init_centroids)
from sklearn.cluster._c_means import (_calculate_memberships,
                                      _calculate_centers)

from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics.pairwise import euclidean_distances

warnings.filterwarnings('ignore', category=RuntimeWarning)

centers = np.array([
    [0.0, 5.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 4.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 5.0, 1.0],
])
n_samples = 100
n_clusters, n_features = centers.shape
X, true_labels = make_blobs(n_samples=n_samples, centers=centers,
                            cluster_std=1., random_state=42)


def test_n_init_error():
    cm = CMeans(n_init=0)
    assert_raise_message(ValueError,
                         'Number of initializations should be a positive '
                         'number, got 0 instead.',
                         cm.fit, X)
    cm = CMeans(n_init=-1)
    assert_raise_message(ValueError,
                         'Number of initializations should be a positive '
                         'number, got -1 instead.',
                         cm.fit, X)


def test_max_iter_error():
    cm = CMeans(max_iter=0)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number,'
                         ' got 0 instead.',
                         cm.fit, X)
    cm = CMeans(max_iter=-1)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number, '
                         'got -1 instead.',
                         cm.fit, X)


def test_copyx():
    # Check if copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    CMeans(copy_x=False, n_clusters=n_clusters, random_state=42).fit(my_X)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)


def test_init_centroids_bounds():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    assert_array_almost_equal(c, np.nan_to_num(c))


def test_calculate_memberships_shape():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = _calculate_memberships(distances, 2)
    assert_equal(m.shape, (n_samples, n_clusters))


def test_calculate_memberships_bounds():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = _calculate_memberships(distances, 2)
    assert_array_less(m, 1e-6 + np.ones_like(m))
    assert_array_less(np.zeros_like(m) - 1e-6, m)


def test_calculate_centers_shape():
    distances = euclidean_distances(X, centers)
    m = _calculate_memberships(distances, 2)
    c = _calculate_centers(X, m, 2)
    assert_equal(c.shape, (n_clusters, n_features))


def test_probabilistic_membership_shape():
    m, i, c, _ = _cmeans_single_probabilistic(X, n_clusters)
    assert_equal(m.shape, (n_samples, n_clusters))


def test_probabilistic_center_shape():
    m, i, c, _ = _cmeans_single_probabilistic(X, n_clusters)
    assert_equal(c.shape, (n_clusters, n_features))
