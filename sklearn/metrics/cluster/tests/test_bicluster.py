"""Testing for bicluster metrics module"""

import numpy as np

from sklearn.utils.testing import assert_equal

from sklearn.metrics.cluster.bicluster import _jaccard
from sklearn.metrics import consensus_score


def test_jaccard():
    a1 = np.array([True, True, False, False])
    a2 = np.array([True, True, True, True])
    a3 = np.array([False, True, True, False])
    a4 = np.array([False, False, True, True])

    assert_equal(_jaccard(a1, a1, a1, a1), 1)
    assert_equal(_jaccard(a1, a1, a2, a2), 0.25)
    assert_equal(_jaccard(a1, a1, a3, a3), 1.0 / 7)
    assert_equal(_jaccard(a1, a1, a4, a4), 0)


def test_consensus_score():
    a = [[True, True, False, False],
         [False, False, True, True]]
    b = a[::-1]

    assert_equal(consensus_score((a, a), (a, a)), 1)
    assert_equal(consensus_score((a, a), (b, b)), 1)
    assert_equal(consensus_score((a, b), (a, b)), 1)
    assert_equal(consensus_score((a, b), (b, a)), 1)

    assert_equal(consensus_score((a, a), (b, a)), 0)
    assert_equal(consensus_score((a, a), (a, b)), 0)
    assert_equal(consensus_score((b, b), (a, b)), 0)
    assert_equal(consensus_score((b, b), (b, a)), 0)
