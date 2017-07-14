"""
Testing for Clustering methods
"""

from sklearn.cluster import MarkovClustering
from sklearn.utils.testing import assert_array_equal

def test_cosine_dissimilarity():
    A = [[1, 2, 0], [2, 1, 0], [0, 0, 4], [0, 0, 9]]
    predicted = MarkovClustering().fit(A).labels_
    #Making sure that element 0 & 1 are in the same cluster, and 2 & 3 are in an other
    assert_array_equal([predicted[0], predicted[0], predicted[3], predicted[3]], predicted)
