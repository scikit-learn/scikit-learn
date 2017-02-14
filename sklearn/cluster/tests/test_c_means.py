"""Testing for C-means"""
import warnings

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import (assert_array_almost_equal,
                                   assert_array_less)

from sklearn.cluster.c_means_ import (probabilistic_single,
                                      possibilistic_single,
                                      ProbabilisticCMeans,
                                      PossibilisticCMeans,
                                      _init_centroids)
from sklearn.cluster._c_means import Probabilistic, Possibilistic

from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.cluster import v_measure_score

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
    cm = ProbabilisticCMeans(n_init=0)
    assert_raise_message(ValueError,
                         'Number of initializations should be a positive '
                         'number, got 0 instead.',
                         cm.fit, X)
    cm = ProbabilisticCMeans(n_init=-1)
    assert_raise_message(ValueError,
                         'Number of initializations should be a positive '
                         'number, got -1 instead.',
                         cm.fit, X)


def test_max_iter_error():
    cm = ProbabilisticCMeans(max_iter=0)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number,'
                         ' got 0 instead.',
                         cm.fit, X)
    cm = ProbabilisticCMeans(max_iter=-1)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number, '
                         'got -1 instead.',
                         cm.fit, X)


def test_copyx():
    # Check if copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    ProbabilisticCMeans(copy_x=False, n_clusters=n_clusters, random_state=42).fit(my_X)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)


def test_init_centroids_bounds():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    assert_array_almost_equal(c, np.nan_to_num(c))


def test_probabilistic_update_memberships_shape():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = Probabilistic.memberships(distances, 2)
    assert_equal(m.shape, (n_samples, n_clusters))


def test_probabilistic_update_memberships_bounds():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = Probabilistic.memberships(distances, 2)
    assert_array_less(m, 1e-6 + np.ones_like(m))
    assert_array_less(np.zeros_like(m) - 1e-6, m)


def test_probabilistic_update_centers_shape():
    distances = euclidean_distances(X, centers)
    m = Probabilistic.memberships(distances, 2)
    c = Probabilistic.centers(X, m, 2)
    assert_equal(c.shape, (n_clusters, n_features))


def test_probabilistic_membership_shape():
    m = probabilistic_single(X, n_clusters)['memberships']
    assert_equal(m.shape, (n_samples, n_clusters))


def test_probabilistic_center_shape():
    c = probabilistic_single(X, n_clusters)['centers']
    assert_equal(c.shape, (n_clusters, n_features))


def test_possibilistic_update_memberships_shape():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = Possibilistic.memberships(distances, 2)
    assert_equal(m.shape, (n_samples, n_clusters))


def test_possibilistic_update_memberships_bounds():
    c = _init_centroids(X, n_clusters, 'random', random_state=42)
    distances = euclidean_distances(X, c)
    m = Possibilistic.memberships(distances, 2)
    assert_array_less(m, 1e-6 + np.ones_like(m))
    assert_array_less(np.zeros_like(m) - 1e-6, m)


def test_possibilistic_update_centers_shape():
    distances = euclidean_distances(X, centers)
    m = Possibilistic.memberships(distances, 2)
    c = Possibilistic.centers(X, m, 2)
    assert_equal(c.shape, (n_clusters, n_features))


def test_possibilistic_membership_shape():
    m = possibilistic_single(X, n_clusters)['memberships']
    assert_equal(m.shape, (n_samples, n_clusters))


def test_possibilistic_center_shape():
    c = possibilistic_single(X, n_clusters)['centers']
    assert_equal(c.shape, (n_clusters, n_features))
    

def test_labels():
    cm = ProbabilisticCMeans()
    cm.memberships_ = np.array([
        [1.0, 0.0],
        [0.7, 0.3],
        [0.4, 0.6],
        [0.2, 0.3],
        [0.5, 0.5],
    ])
    labels_true = np.array([0, 0, 1, 1, 0])
    assert_array_almost_equal(cm.labels_, labels_true)


def _check_fitted_model(cm):
    assert_equal(v_measure_score(true_labels, cm.labels_), 1.0)
    assert_greater(cm.inertia_, 0.0)


def test_probabilistic_results():
    cm = ProbabilisticCMeans(n_clusters=3, algorithm="probabilistic", random_state=4)
    cm.fit(X)
    _check_fitted_model(cm)


def test_possibilistic_results():
    cm = PossibilisticCMeans(n_clusters=3, algorithm="possibilistic", random_state=4)
    cm.fit(X)
    _check_fitted_model(cm)


def test_probabilistic_predict():
    cm = ProbabilisticCMeans(n_clusters=3, algorithm="probabilistic", random_state=4)
    cm.fit(X)
    pred = cm.predict(centers)
    assert_equal(pred.shape, (n_clusters, n_clusters))
    assert_equal(v_measure_score(np.arange(n_clusters), np.argmax(pred, axis=1)), 1.0)


def test_possibilistic_predict():
    cm = PossibilisticCMeans(n_clusters=3, algorithm="possibilistic", random_state=4)
    cm.fit(X)
    pred = cm.predict(centers)
    assert_equal(pred.shape, (n_clusters, n_clusters))
    assert_equal(v_measure_score(np.arange(n_clusters), np.argmax(pred, axis=1)), 1.0)


