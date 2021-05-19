import unittest
import numpy as np
import scipy.sparse as sp

from .._bisect_k_means import BisectKMeans
from numpy.testing import assert_array_equal, assert_array_almost_equal

# import numpy.random
# import time
# import statistics
# from .._kmeans import KMeans


class MyTestCase(unittest.TestCase):
    def test_two_clusters(self):
        X = np.array([[1, 2], [1, 4], [1, 0],
                      [10, 2], [10, 4], [10, 0]])

        bisect_means = BisectKMeans(n_clusters=2, random_state=0)
        bisect_means.fit(X)

        expected_centers = [[1, 2], [10, 2]]
        assert_array_equal(expected_centers, bisect_means.cluster_centers_)
        assert_array_equal(bisect_means.predict([[0, 0], [12, 3]]), [0, 1])

    def test_three_clusters(self):
        X = np.array([[1, 2], [1, 4], [1, 0],
                      [10, 2], [10, 4], [10, 0],
                      [10, 6], [10, 8], [10, 10]])
        bisect_means = BisectKMeans(n_clusters=3, random_state=0)
        bisect_means.fit(X)

        expected_centers = [[1, 2], [10, 8], [10, 2]]
        assert_array_equal(expected_centers, bisect_means.cluster_centers_)
        assert_array_equal(bisect_means.predict([[0, 0], [12, 3]]), [0, 2])

    # def test_performance(self):
    #     X = np.random.randint(10, size=(100, 2))
    #     results = []
    #     times = []
    #     bisect_means = BisectKMeans(n_clusters=4, random_state=0)
    #     kmeans = KMeans(n_clusters=4, random_state=0)
    #
    #     for i in range(10):
    #         start_time = time.time()
    #
    #         kmeans.fit(X)
    #
    #         point1 = (time.time() - start_time)
    #
    #         start_time = time.time()
    #
    #         bisect_means.fit(X)
    #
    #         point2 = (time.time() - start_time)
    #
    #         msg = "k2" if point2 < point1 else "k1"
    #         times.append(point1 - point2)
    #         results.append(msg)
    #
    #     k1 = len([x for x in results if x == "k1"])
    #     k2 = len(results) - k1
    #     print("Biscect K-Means:{}, K-Means:{}".format(k1, k2))
    #     print("AVG:", statistics.mean(times))
    #     print("Bisect" if statistics.mean(times) < 0 else "Kmeans")

    def test_sparse(self):
        rng = np.random.RandomState(0)

        X = rng.rand(40, 2)
        X[X < .8] = 0
        X_csr = sp.csr_matrix(X)

        bisect_means = BisectKMeans(n_clusters=3, random_state=0)
        bisect_means.fit(X_csr)

        sparse_centers = bisect_means.cluster_centers_

        bisect_means.fit(X)
        assert_array_almost_equal(bisect_means.cluster_centers_, sparse_centers)


if __name__ == '__main__':
    unittest.main()
