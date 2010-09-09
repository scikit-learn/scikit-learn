"""
Testing for K-means. 

"""

import numpy as np
from numpy.testing import assert_equal

from ..k_means_ import KMeans
from .common import generate_clustered_data

n_clusters = 3
X = generate_clustered_data(n_clusters=n_clusters, std=.1)


def test_k_means():
    np.random.seed(1)
    k_means = KMeans().fit(X, k=3) 
    cluster_centers_indices = k_means.cluster_centers_
    labels = k_means.labels_

    assert_equal(np.unique(labels[:20]).size, 1)
    assert_equal(np.unique(labels[20:40]).size, 1)
    assert_equal(np.unique(labels[40:]).size, 1)
