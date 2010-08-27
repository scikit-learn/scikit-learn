"""
Testing for mean shift clustering methods

"""

import numpy as np
from numpy.testing import assert_equal

from .. import MeanShift, mean_shift
from .common import generate_clustered_data

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters)

def test_mean_shift():
    """ Test MeanShift algorithm
    """
    bandwidth = 1.2

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

