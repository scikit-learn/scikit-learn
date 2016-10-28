"""
Testing for mean shift clustering methods

"""

import numpy as np
import warnings

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message

from sklearn.cluster import MeanShift
from sklearn.cluster import mean_shift
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import get_bin_seeds
from sklearn.datasets.samples_generator import make_blobs


n_clusters = 3
centers = np.array([[1, 1], [-1, -1], [1, -1]]) + 10
X, _ = make_blobs(n_samples=300, n_features=2, centers=centers,
                  cluster_std=0.4, shuffle=True, random_state=11)


def test_estimate_bandwidth():
    # Test estimate_bandwidth
    bandwidth = estimate_bandwidth(X, n_samples=200)
    assert_true(0.9 <= bandwidth <= 1.5)


def test_mean_shift():
    # Test MeanShift algorithm
    bandwidth = 1.2

    ms = MeanShift(bandwidth=bandwidth)
    labels = ms.fit(X).labels_
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    assert_equal(n_clusters_, n_clusters)

    cluster_centers, labels = mean_shift(X, bandwidth=bandwidth)
    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)
    assert_equal(n_clusters_, n_clusters)


def test_parallel():
    ms1 = MeanShift(n_jobs=2)
    ms1.fit(X)

    ms2 = MeanShift()
    ms2.fit(X)

    assert_array_equal(ms1.cluster_centers_, ms2.cluster_centers_)
    assert_array_equal(ms1.labels_, ms2.labels_)


def test_meanshift_predict():
    # Test MeanShift.predict
    ms = MeanShift(bandwidth=1.2)
    labels = ms.fit_predict(X)
    labels2 = ms.predict(X)
    assert_array_equal(labels, labels2)


def test_meanshift_all_orphans():
    # init away from the data, crash with a sensible warning
    ms = MeanShift(bandwidth=0.1, seeds=[[-9, -9], [-10, -10]])
    msg = "No point was within bandwidth=0.1"
    assert_raise_message(ValueError, msg, ms.fit, X,)


def test_unfitted():
    # Non-regression: before fit, there should be not fitted attributes.
    ms = MeanShift()
    assert_false(hasattr(ms, "cluster_centers_"))
    assert_false(hasattr(ms, "labels_"))


def test_bin_seeds():
    # Test the bin seeding technique which can be used in the mean shift
    # algorithm
    # Data is just 6 points in the plane
    X = np.array([[1., 1.], [1.4, 1.4], [1.8, 1.2],
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
    # we bail and use the whole data here.
    with warnings.catch_warnings(record=True):
        test_bins = get_bin_seeds(X, 0.01, 1)
    assert_array_equal(test_bins, X)

    # tight clusters around [0, 0] and [1, 1], only get two bins
    X, _ = make_blobs(n_samples=100, n_features=2, centers=[[0, 0], [1, 1]],
                      cluster_std=0.1, random_state=0)
    test_bins = get_bin_seeds(X, 1)
    assert_array_equal(test_bins, [[0, 0], [1, 1]])
