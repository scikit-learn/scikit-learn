"""
Tests for ig
"""

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from sklearn.feature_selection import SelectKBest, ig

from nose.tools import assert_raises
from numpy.testing import assert_equal

# Feature 0 is highly informative for class 1;
# feature 1 is the same everywhere;
# feature 2 is a bit informative for class 2.
X = [[2, 1, 2],
     [9, 1, 1],
     [6, 1, 2],
     [0, 1, 2]]
y = [0, 1, 2, 2]


def mk_ig(k):
    """Make k-best IG selector"""
    return SelectKBest(ig, k=k)


def test_ig():
    # Test IG feature extraction

    Xsp = csr_matrix(X, dtype=np.float)
    scores = mk_ig(k=2).fit(Xsp, y)
    assert_equal(sorted(scores.get_support(indices=True)), [0, 2])
    Xtrans = scores.transform(Xsp)
    assert_equal(Xtrans.shape, [Xsp.shape[0], 2])

    # == doesn't work on scipy.sparse matrices
    Xtrans = Xtrans.toarray()
    Xtrans2 = mk_ig(k=2).fit_transform(Xsp, y).toarray()
    assert_equal(Xtrans, Xtrans2)


def test_ig_coo():
    # Check that ig works with a COO matrix
    # (as returned by CountVectorizer, DictVectorizer)
    Xcoo = coo_matrix(X)
    mk_ig(k=2).fit_transform(Xcoo, y)
    # if we got here without an exception, we're safe


def test_ig_negative():
    # Check for proper error on negative numbers in the input X.
    X, y = [[0, 1], [-1e-20, 1]], [0, 1]
    assert_raises(ValueError, ig, csr_matrix(X), y)


def test_ig_nonsparse_matrix():
    assert_raises(ValueError, mk_ig(k=2).fit, X, y)
