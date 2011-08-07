"""
Testing for mean shift clustering methods

"""

import numpy as np
from numpy.testing import assert_equal

from .. import MeanShift, mean_shift, estimate_bandwidth
from ...datasets.samples_generator import make_blobs

n_clusters = 3

centers = np.array([[ 1,  1, 1, 0],
                    [-1, -1, 0, 1],
                    [ 1, -1, 1, 1]]) + 10

X, _ = make_blobs(n_samples=500, n_features=2, centers=centers,
                  cluster_std=0.4, shuffle=True, random_state=0)


def test_mean_shift():
    """ Test MeanShift algorithm
    """
    bandwidth = 1.2

    bandwidth_ = estimate_bandwidth(X, n_samples=300)
    assert 1.1 <= bandwidth_ <= 1.5

    ms = MeanShift(bandwidth=bandwidth)
    labels = ms.fit(X).labels_
    cluster_centers = ms.cluster_centers_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    assert_equal(n_clusters_, n_clusters)

    cluster_centers, labels = mean_shift(X, bandwidth=bandwidth)
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    assert_equal(n_clusters_, n_clusters)
