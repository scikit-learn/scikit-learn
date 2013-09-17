"""Testing for bicluster metrics module"""

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true

from ..bicluster_metrics import _jaccard
from ..bicluster_metrics import _make_similarity
from ..bicluster_metrics import consensus_score


def test_jaccard():
    a1 = np.array([True, True, False, False])
    a2 = np.array([True, True, True, True])
    a3 = np.array([False, True, True, False])
    a4 = np.array([False, False, True, True])

    f = _make_similarity(_jaccard, None)

    assert_equal(f(a1, a1, a1, a1), 1)
    assert_equal(f(a1, a1, a2, a2), 0.25)
    assert_equal(f(a1, a1, a3, a3), 1.0 / 7)
    assert_equal(f(a1, a1, a4, a4), 0)


def test_consensus_score():
    a = [[True, True, False, False],
         [False, False, True, True]]
    b = a[::-1]

    for method in ('jaccard', 'goodness', 'dice'):
        for correction in (None, 16):
            assert_equal(consensus_score((a, a), (a, a),
                                         method, correction), 1)
            assert_equal(consensus_score((a, a), (b, b),
                                         method, correction), 1)
            assert_equal(consensus_score((a, b), (a, b),
                                         method, correction), 1)
            assert_equal(consensus_score((a, b), (b, a),
                                         method, correction), 1)

            assert_true(consensus_score((a, a), (b, a),
                                        method, correction) <= 0)
            assert_true(consensus_score((a, a), (a, b),
                                        method, correction) <= 0)
            assert_true(consensus_score((b, b), (a, b),
                                        method, correction) <= 0)
            assert_true(consensus_score((b, b), (b, a),
                                        method, correction) <= 0)


def test_consensus_score_mismatched():
    """Test consensus score with differing number of biclusters."""
    a = [[True, True, False, False],
         [False, False, True, True]]
    b = [[True, True, False, False],
         [True, True, True, False],
         [True, False, True, False],
         [False, False, True, True]]
    assert_equal(consensus_score((a, a), (b, b)), 0.5)


def test_consensus_score_missing():
    """Test consensus score with empty biclusters"""
    a = [[True, True, False, False],
         [False, False, True, True],
         [False, False, False, False]]
    b = [[True, True, False, False],
         [True, True, True, False],
         [True, False, True, False],
         [False, False, True, True]]
    assert_equal(consensus_score((a, a), (b, b)), 0.5)
