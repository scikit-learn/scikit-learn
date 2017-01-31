"""Testing for C-means"""

import numpy as np

from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_array_almost_equal

from sklearn.cluster import CMeans, c_means
from sklearn.datasets.samples_generator import make_blobs


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
                         'Number of initializations should be a positive number,'
                         ' got 0 instead.',
                         cm.fit, X)
    cm = CMeans(n_init=-1)
    assert_raise_message(ValueError,
                         'Number of initializations should be a positive number, got -1 instead.',
                         cm.fit, X)

def test_max_iter_error():
    cm = CMeans(max_iter=0)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number,'
                         ' got 0 instead.',
                         cm.fit, X)
    cm = CMeans(max_iter=-1)
    assert_raise_message(ValueError,
                         'Number of iterations should be a positive number, got -1 instead.',
                         cm.fit, X)


def test_copyx():
    # Check if copy_x=False returns nearly equal X after de-centering.
    my_X = X.copy()
    cm = CMeans(copy_x=False, n_clusters=n_clusters, random_state=42).fit(my_X)

    # check if my_X is centered
    assert_array_almost_equal(my_X, X)
