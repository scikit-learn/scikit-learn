"""Testing for bicluster metrics module"""

import numpy as np

from sklearn.utils.testing import assert_equal

from ..bicluster_metrics import _jaccard
from ..bicluster_metrics import score_biclusters


def test_jaccard():
    a1 = np.array([True, True, False, False])
    a2 = np.array([True, True, True, True])
    a3 = np.array([False, True, True, False])
    a4 = np.array([False, False, True, True])

    assert_equal(_jaccard(a1, a1, a1, a1), 1)
    assert_equal(_jaccard(a1, a1, a2, a2), 0.25)
    assert_equal(_jaccard(a1, a1, a3, a3), 1.0 / 7)
    assert_equal(_jaccard(a1, a1, a4, a4), 0)


def test_score_biclusters():
    a = [[True, True, False, False],
         [False, False, True, True]]
    b = a[::-1]

    assert_equal(score_biclusters((a, a), (a, a)), 1)
    assert_equal(score_biclusters((a, a), (b, b)), 1)
    assert_equal(score_biclusters((a, b), (a, b)), 1)
    assert_equal(score_biclusters((a, b), (b, a)), 1)

    assert_equal(score_biclusters((a, a), (b, a)), 0)
    assert_equal(score_biclusters((a, a), (a, b)), 0)
    assert_equal(score_biclusters((b, b), (a, b)), 0)
    assert_equal(score_biclusters((b, b), (b, a)), 0)
