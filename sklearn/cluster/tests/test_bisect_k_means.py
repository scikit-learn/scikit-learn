import unittest
import numpy as np
from .._bisect_k_means import BisectKMeans
from numpy.testing import assert_array_equal


class MyTestCase(unittest.TestCase):
    def test_bicluster(self):
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



if __name__ == '__main__':
    unittest.main()
