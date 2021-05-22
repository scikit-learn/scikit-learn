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
        """ Test Bisecting K-Means for two clusters
            This should also result with
        """

        X = np.array([[1, 2], [1, 4], [1, 0],
                      [10, 2], [10, 4], [10, 0]])

        bisect_means = BisectKMeans(n_clusters=2, random_state=0)
        bisect_means.fit(X)

        expected_centers = [[1, 2], [10, 2]]
        assert_array_equal(expected_centers, bisect_means.cluster_centers_)
        assert_array_equal(bisect_means.predict([[0, 0], [12, 3]]), [0, 1])

    def test_three_clusters(self):
        """ Tries to perform bisect k-means for three clusters to check
            if splitting data is performed correctly
        """

        X = np.array([[1, 2], [1, 4], [1, 0],
                      [10, 2], [10, 4], [10, 0],
                      [10, 6], [10, 8], [10, 10]])
        bisect_means = BisectKMeans(n_clusters=3, random_state=0)
        bisect_means.fit(X)

        expected_centers = [[1, 2], [10, 8], [10, 2]]
        assert_array_equal(expected_centers, bisect_means.cluster_centers_)
        assert_array_equal(bisect_means.predict([[0, 0], [12, 3]]), [0, 2])

    def test_sparse(self):
        """ Test Bisecting K-Means with sparse data """

        rng = np.random.RandomState(0)

        X = rng.rand(40, 2)
        X[X < .8] = 0
        X_csr = sp.csr_matrix(X)

        bisect_means = BisectKMeans(n_clusters=3, random_state=0)
        bisect_means.fit(X_csr)

        sparse_centers = bisect_means.cluster_centers_

        bisect_means.fit(X)
        normal_centers = bisect_means.cluster_centers_
        assert_array_almost_equal(normal_centers, sparse_centers)


if __name__ == '__main__':
    unittest.main()
