"""
Testing for Clustering methods
"""

from sklearn.cluster import MarkovClustering
from sklearn.utils.testing import assert_array_equal


def test_cosine_dissimilarity():
    A = [[1, 2, 0], [2, 1, 0], [0, 0, 4], [0, 0, 9]]
    l = MarkovClustering().fit(A).labels_
    # Making sure that element 0 & 1 are in the same cluster,
    # and that elements 2 & 3 are in a second one
    assert_array_equal([l[0], l[0], l[3], l[3]], l)
