from numpy.testing import assert_array_almost_equal

from ..pairwise import euclidean_distances


def test_euclidean_distances():
    """Check that the pairwise euclidian distances computation"""
    X = [[0]]
    Y = [[1], [2]]
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])
   
