"""
Tests for info_gain
"""

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from sklearn.feature_selection import SelectKBest, info_gain, info_gain_ratio

from sklearn.utils.testing import assert_raises_regex
from numpy.testing import assert_equal, assert_almost_equal

# Feature 0 is highly informative for class 1;
# feature 1 is the same everywhere;
# feature 2 is a bit informative for class 2.
X = [[2, 1, 2],
     [9, 1, 1],
     [6, 1, 2],
     [0, 1, 2]]
y = [0, 1, 1, 0]


def mk_info_gain(k):
    """Make k-best IG selector"""
    return SelectKBest(info_gain, k=k)


def test_info_gain_csr():
    # Test IG feature extraction

    Xsp = csr_matrix(X, dtype=np.float)
    scores = mk_info_gain(k=2).fit(Xsp, y)
    assert_equal(sorted(scores.get_support(indices=True)), [0, 2])
    Xtrans = scores.transform(Xsp)
    assert_equal(Xtrans.shape, [Xsp.shape[0], 2])

    # == doesn't work on scipy.sparse matrices
    Xtrans = Xtrans.toarray()
    Xtrans2 = mk_info_gain(k=2).fit_transform(Xsp, y).toarray()
    assert_equal(Xtrans, Xtrans2)


def test_info_gain_coo():
    # Check that ig works with a COO matrix
    # (as returned by CountVectorizer, DictVectorizer)
    Xcoo = coo_matrix(X)
    mk_info_gain(k=2).fit_transform(Xcoo, y)
    # if we got here without an exception, we're safe


def test_info_gain_dense():
    # Check IG works with a dense matrix

    Xden = np.array(X)
    scores = mk_info_gain(k=2).fit(Xden, y)
    assert_equal(sorted(scores.get_support(indices=True)), [0, 2])

    Xtrans = scores.transform(Xden)
    assert_equal(Xtrans.shape, [Xden.shape[0], 2])

    Xtrans2 = mk_info_gain(k=2).fit_transform(Xden, y)
    assert_equal(Xtrans, Xtrans2)


def test_info_gain_negative():
    # Check for proper error on negative numbers in the input X.
    X, y = [[0, 1], [-1e-20, 1]], [0, 1]
    assert_raises_regex(ValueError, "Input X must be non-negative.", info_gain,
                        csr_matrix(X), y)


def test_expected_value_info_gain():
    # Check calculation of an expected value of IG

    # two instances, three features, two classes
    X, y = [[1, 5, 9], [9, 5, 1]], [0, 1]

    # Counts:
    # f1: 10
    # f2: 10
    # f3: 10
    # c1: 15
    # c2: 15
    # total: 30

    # Probabilities:
    # f1: 0.33333
    # f2: 0.33333
    # c1: 0.5
    # c2: 0.5
    # f1, c1: 0.03333
    # nf1, c1: 0.46666
    # f1, c2: 0.3
    # nf1, c2: 0.2

    # Class-specific IG scores for f1:

    # f1, c1:
    # 0.03333 * log (0.03333 / (0.5 * 0.33333)) = -0.07739
    # f1, n_c1:
    # 0.3 * log (0.3 / (0.5 * 0.33333)) = 0.2544
    # n_f1, c1:
    # 0.46666 * log (0.46666 / (0.5 * 0.66666)) = 0.22652
    # n_f1, n_c1:
    # 0.2 * log (0.2 / (0.5 * 0.66666)) = -0.14739
    # sum:
    # -0.07739 + 0.2544 + 0.22652 + -0.14739 = 0.25614

    # f1, c2:
    # 0.3 * log (0.3 / (0.5 * 0.33333)) = 0.2544
    # f1, n_c2:
    # 0.03333 * log (0.03333 / (0.5 * 0.33333)) = -0.07739
    # n_f1, c2:
    # 0.2 * log (0.2 / (0.5 * 0.66666)) = -0.14739
    # n_f1, n_c2:
    # 0.46666 * log (0.46666 / (0.5 * 0.66666)) = 0.22652
    # sum:
    # -0.07739 + 0.2544 + 0.22652 + -0.14739 = 0.25614

    # Expected global max score for f1: max (0.25614, 0.25614)

    Xsp = csr_matrix(X, dtype=np.float)
    scores, probs = info_gain(Xsp, y)
    assert_almost_equal(scores[0], 0.25614, decimal=5)


def test_expected_value_info_gain_ratio():
    # Check calculation of an expected value of IGR

    # two instances, three features, two classes
    X, y = [[1, 5, 9], [9, 5, 1]], [0, 1]

    # Entropy of c1: -0.5 * log(0.5) = 0.5
    # Entropy of c2: -0.5 * log(0.5) = 0.5

    # Class-specific IG scores, calculated as in
    # `test_expected_value_info_gain`:
    # IG (f1, c1),
    # -0.07739 + 0.2544 + 0.22652 + -0.14739 = 0.25614
    # Normalize:
    # 0.25614 / (0.5 + 0.5) = 0.25614

    # IG (f1, c2):
    # -0.07739 + 0.2544 + 0.22652 + -0.14739 = 0.25614
    # Normalize:
    # 0.25614 / (0.5 + 0.5) = 0.25614

    # Expected global max score for f1: max (0.25614, 0.25614)

    Xsp = csr_matrix(X, dtype=np.float)
    scores, probs = info_gain_ratio(Xsp, y)
    assert_almost_equal(scores[0], 0.25614, decimal=5)
