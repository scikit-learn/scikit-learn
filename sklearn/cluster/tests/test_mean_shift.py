"""
Testing for mean shift clustering methods

"""

import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_true

from .. import MeanShift, mean_shift, estimate_bandwidth, get_bin_seeds
from ...datasets.samples_generator import make_blobs

n_clusters = 3
centers = np.array([[1, 1], [-1, -1], [1, -1]]) + 10
X, _ = make_blobs(n_samples=500, n_features=2, centers=centers,
                  cluster_std=0.4, shuffle=True, random_state=0)


def test_mean_shift():
    """ Test MeanShift algorithm
    """
    bandwidth = 1.2

    bandwidth_ = estimate_bandwidth(X, n_samples=300)
    assert_true(0.9 <= bandwidth_ <= 1.5)

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


def test_bin_seeds():
    """
    Test the bin seeding technique which can be used in the mean shift
    algorithm
    """
    # Data is just 6 points in the plane
    X = np.array([[1., 1.], [1.5, 1.5], [1.8, 1.2],
                  [2., 1.], [2.1, 1.1], [0., 0.]])

    # With a bin coarseness of 1.0 and min_bin_freq of 1, 3 bins should be
    # found
    ground_truth = set([(1., 1.), (2., 1.), (0., 0.)])
    test_bins = get_bin_seeds(X, 1, 1)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(ground_truth.symmetric_difference(test_result)) == 0)

    # With a bin coarseness of 1.0 and min_bin_freq of 2, 2 bins should be
    # found
    ground_truth = set([(1., 1.), (2., 1.)])
    test_bins = get_bin_seeds(X, 1, 2)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(ground_truth.symmetric_difference(test_result)) == 0)

    # With a bin size of 0.01 and min_bin_freq of 1, 6 bins should be found
    test_bins = get_bin_seeds(X, 0.01, 1)
    test_result = set([tuple(p) for p in test_bins])
    assert_true(len(test_result) == 6)
