import numpy as np
import scipy.sparse as sp

from sklearn.feature_extraction import FeatureHasher

from nose.tools import assert_raises, assert_true
from numpy.testing import assert_array_equal, assert_equal
from sklearn.utils.testing import assert_in


def test_feature_hasher():
    raw_X = [[u"foo", "bar", "baz", "foo"],    # note: duplicate
             [u"bar", "baz", "quux"]]

    for lg_n_features in (7, 9, 11, 16, 22):
        n_features = 2 ** lg_n_features

        it = (x for x in raw_X)                 # iterable

        h = FeatureHasher(n_features, non_negative=True)
        X = h.transform(it)

        assert_equal(X.shape[0], len(raw_X))
        assert_equal(X.shape[1], n_features)

        assert_true(np.all(X.data > 0))
        assert_equal(X[0].sum(), 4)
        assert_equal(X[1].sum(), 3)

        # .nnz is unreliable on coo_matrix
        assert_equal(X.tocsr().nnz, sum(len(set(x)) for x in raw_X))


def test_hash_empty_input():
    n_features = 16
    raw_X = [[], (), xrange(0)]

    h = FeatureHasher(n_features=n_features)
    X = h.transform(raw_X)

    assert_array_equal(X.A, np.zeros((len(raw_X), n_features)))


def test_hasher_invalid_input():
    assert_raises(ValueError, FeatureHasher, n_features=-1)
    assert_raises(ValueError, FeatureHasher, n_features=0)
    assert_raises(TypeError, FeatureHasher, n_features='ham')

    h = FeatureHasher(n_features=np.uint16(2**6))
    assert_raises(ValueError, h.transform, [])
    assert_raises(TypeError, h.transform, [[5.5]])
    assert_raises(TypeError, h.transform, [[None]])
