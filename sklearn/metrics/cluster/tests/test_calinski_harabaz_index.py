import numpy as np

from sklearn.utils.testing import (assert_almost_equal, assert_equal)

from sklearn.cluster.k_means_ import KMeans
from sklearn.metrics.cluster.calinski_harabaz_index import (calinski_harabaz_index,
                                                            max_CH_index)
from sklearn.datasets import make_blobs


def test_calinsk_harabaz_index():
    X = np.array([[0, 0]] * 50 + [[1, 1]] * 50
                 + [[3, 3]] * 50 + [[4, 4]] * 50)
    labels = [0] * 100 + [1] * 100
    assert_almost_equal(calinski_harabaz_index(X, labels), 4.5 * (200 - 2) / .5)


def test_Calinsk_Harabasz():
    X, _ = make_blobs(90, centers=np.array([[-2, -2], [2, 0], [-2, 2]]), random_state=0)
    cluster_estimator = KMeans()
    assert_equal(max_CH_index(X, cluster_estimator, k_max=6), 3)
