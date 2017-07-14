"""
Testing for Clustering methods
"""

from sklearn.cluster import MarkovClustering
from sklearn.utils.testing import assert_array_equal

def test_cosine_dissimilarity():
    A = [[1, 2, 0], [2, 1, 0], [0, 0, 4], [0, 0, 9]]
    predicted = MarkovClustering().fit(A, verbose=True).labels_
    assert_array_equal([1, 1, 0, 0], predicted)
