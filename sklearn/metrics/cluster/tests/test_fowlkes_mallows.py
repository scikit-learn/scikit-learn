from __future__ import division
from math import sqrt

from sklearn.utils.testing import assert_almost_equal

from sklearn.metrics.cluster.fowlkes_mallows import fowlkes_mallows_index


def test_fowlkes_mallows_index():
    clustering_1 = [0, 0, 0, 1, 1, 1]
    clustering_2 = [0, 0, 1, 1, 2, 2]
    assert_almost_equal(
        fowlkes_mallows_index(clustering_1, clustering_2),
        10 / sqrt(18 * 12))
