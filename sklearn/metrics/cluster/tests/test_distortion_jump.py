import numpy as np

from sklearn.utils.testing import (assert_almost_equal, assert_equal)

from sklearn.cluster.k_means_ import KMeans
from sklearn.metrics.cluster.distortion_jump import distortion_jump
from sklearn.datasets import make_blobs


def test_distortion_jump():
    X, _ = make_blobs(90, centers=np.array([[-2, -2], [2, 0], [-2, 2]]),
                      random_state=0)
    cluster_estimator = KMeans()
    assert_equal(distortion_jump(X, cluster_estimator, k_max=6), 3)
