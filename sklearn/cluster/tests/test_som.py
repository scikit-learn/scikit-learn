"""
Testing for SOM

"""
import numpy as np
from numpy.testing import assert_equal

from ..som_ import SelfOrganizingMap
from .common import generate_clustered_data

n_clusters = 4
n_features = 2
X2 = generate_clustered_data(n_clusters=n_clusters, n_features=2, std=.1)
X3 = generate_clustered_data(n_clusters=8, n_features=3, std=.1,
                             n_samples_per_cluster=10)


def test_som():
    np.random.seed(1)
    som = SelfOrganizingMap(adjacency=(2, 2), n_iterations=10, learning_rate=1)
    som.fit(X2)
    labels = som.labels_

    assert_equal(np.unique(labels).shape[0], 4)
    assert_equal(np.unique(labels[:20]).shape[0], 1)
    assert_equal(np.unique(labels[20:40]).shape[0], 1)
    assert_equal(np.unique(labels[40:60]).shape[0], 1)
    assert_equal(np.unique(labels[60:]).shape[0], 1)


def test_som_init_matrix():
    np.random.seed(1)
    random_ind = np.random.randint(0, X2.shape[0], size=n_clusters)
    init_nodes = X2[random_ind]

    som = SelfOrganizingMap(adjacency=(2, 2), init=init_nodes,
                            learning_rate=0.1)

    som.fit(X2)
    labels = som.labels_
    assert_equal(np.unique(labels).shape[0], 4)
    assert_equal(np.unique(labels[:20]).shape[0], 1)
    assert_equal(np.unique(labels[20:40]).shape[0], 1)
    assert_equal(np.unique(labels[40:60]).shape[0], 1)
    assert_equal(np.unique(labels[60:]).shape[0], 1)


def test_som_one_dimensional():
    np.random.seed(1)

    som = SelfOrganizingMap(adjacency=(4,))

    som.fit(X2)
    labels = som.labels_
    assert_equal(np.unique(labels).shape[0], 4)
    assert_equal(np.unique(labels[:20]).shape[0], 1)
    assert_equal(np.unique(labels[20:40]).shape[0], 1)
    assert_equal(np.unique(labels[40:60]).shape[0], 1)
    assert_equal(np.unique(labels[60:]).shape[0], 1)


def test_som_3_dimensional():
    np.random.seed(1)

    som = SelfOrganizingMap(adjacency=(2, 2, 2), n_iterations=200)

    som.fit(X3)
    labels = som.labels_
    assert_equal(np.unique(labels).shape[0], 8)
    assert_equal(np.unique(labels[:10]).shape[0], 1)
    assert_equal(np.unique(labels[10:20]).shape[0], 1)
    assert_equal(np.unique(labels[20:30]).shape[0], 1)
    assert_equal(np.unique(labels[30:40]).shape[0], 1)
    assert_equal(np.unique(labels[40:50]).shape[0], 1)
    assert_equal(np.unique(labels[50:60]).shape[0], 1)
    assert_equal(np.unique(labels[60:70]).shape[0], 1)
    assert_equal(np.unique(labels[70:]).shape[0], 1)
