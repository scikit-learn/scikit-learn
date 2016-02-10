import numpy as np
from math import sqrt

from sklearn.utils.testing import assert_almost_equal, assert_equal

from sklearn.cluster.k_means_ import KMeans
from sklearn.metrics.cluster.stability import (
    _one_stability_measure, stability, function_cluster_similarity,
    fowlkes_mallows_index)
from sklearn.datasets import make_blobs


def test_one_stability_measure():
    X = np.arange(10) < 5
    X.reshape((10, 1))

    # test perfect clustering has 1 stability
    class SameCluster(object):
        def set_params(self, *args, **kwargs):
            pass

        def fit_predict(self, X):
            return X
    same_clusters = SameCluster()
    assert_almost_equal(
        _one_stability_measure(same_clusters, X, .8, fowlkes_mallows_index), 1)


def test_stability():
    X, _ = make_blobs(90, centers=np.array([[-2, -2], [2, 0], [-2, 2]]),
                      random_state=0)
    cluster_estimator = KMeans()
    assert_equal(stability(X, cluster_estimator, k_max=6,
                           nb_draw=10, random_state=0), 3)


def test_function_cluster_dist():
    clustering_1 = [0, 0, 0, 1, 1, 1]
    clustering_2 = [0, 0, 1, 1, 2, 2]
    dist_fun = function_cluster_similarity(metric='fowlkes-mallows')
    assert_almost_equal(dist_fun(clustering_1, clustering_2),
                        10 / sqrt(18 * 12))
    dist_fun = function_cluster_similarity(metric='cityblock')
    assert_almost_equal(dist_fun(clustering_1, clustering_2), -10)
    dist_fun = function_cluster_similarity(metric='euclidean')
    assert_almost_equal(dist_fun(clustering_1, clustering_2), -sqrt(10))
