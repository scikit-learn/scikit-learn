import numpy as np

from sklearn.metrics.cluster.adjacency_matrix import adjacency_matrix
from sklearn.utils.testing import assert_array_equal


def test_adjacency_matrix():
    assi = [0, 0, 1, 1]
    adj_matrix = np.array([[1, 1, 0, 0], [1, 1, 0, 0],
                           [0, 0, 1, 1], [0, 0, 1, 1]])
    found_adj = adjacency_matrix(assi)
    assert_array_equal(found_adj, adj_matrix)
