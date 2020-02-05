
import numpy as np
from numpy.testing import assert_array_equal
import pytest

from sklearn.feature_extraction import FeatureLightweightRandomIndexing
from sklearn.utils._testing import fails_if_pypy

pytestmark = fails_if_pypy


def test_feature_lri_dicts():
    h = FeatureLightweightRandomIndexing(n_features=16)
    assert "dict" == h.input_type

    raw_X = [{"foo": "bar", "dada": 42, "tzara": 37},
             {"foo": "baz", "gaga": "string1"}]
    X1 = FeatureLightweightRandomIndexing(n_features=16).transform(raw_X)
    gen = (iter(d.items()) for d in raw_X)
    X2 = FeatureLightweightRandomIndexing(n_features=16,
                                          input_type="pair").transform(gen)
    assert_array_equal(X1.toarray(), X2.toarray())


def test_feature_lri_strings():
    # mix byte and Unicode strings; note that "foo" is a duplicate in row 0
    raw_X = [["foo", "bar", "baz", "foo".encode("ascii")],
             ["bar".encode("ascii"), "baz", "quux"]]

    for lg_n_features in (7, 9, 11, 16, 22):
        n_features = 2 ** lg_n_features

        it = (x for x in raw_X)                 # iterable

        h = FeatureLightweightRandomIndexing(n_features, input_type="string")
        X = h.transform(it)

        assert X.shape[0] == len(raw_X)
        assert X.shape[1] == n_features


def test_lri_transform_seed():
    # check the influence of the seed when computing the hashes
    # import is here to avoid importing on pypy
    from sklearn.feature_extraction._lri_fast import (
            transform as _lri_transform)
    raw_X = [["foo", "bar", "baz", "foo".encode("ascii")],
             ["bar".encode("ascii"), "baz", "quux"]]

    raw_X_ = (((f, 1) for f in x) for x in raw_X)
    indices, indptr, _ = _lri_transform(raw_X_, 2 ** 7, str)

    raw_X_ = (((f, 1) for f in x) for x in raw_X)
    indices_0, indptr_0, _ = _lri_transform(raw_X_, 2 ** 7, str, seed=0)
    assert_array_equal(indices, indices_0)
    assert_array_equal(indptr, indptr_0)

    raw_X_ = (((f, 1) for f in x) for x in raw_X)
    indices_1, _, _ = _lri_transform(raw_X_, 2 ** 7, str, seed=1)
    with pytest.raises(AssertionError):
        assert_array_equal(indices, indices_1)


def test_feature_lri_pairs():
    raw_X = (iter(d.items()) for d in [{"foo": 1, "bar": 2},
                                       {"baz": 3, "quux": 4, "foo": -1}])
    h = FeatureLightweightRandomIndexing(n_features=16, input_type="pair")
    x1, x2 = h.transform(raw_X).toarray()
    x1_nz = set(np.abs(x1[x1 != 0]))
    x2_nz = set(np.abs(x2[x2 != 0]))
    assert {1, 2} == x1_nz
    assert {1, 3, 4} == x2_nz


def test_feature_lri_pairs_with_string_values():
    raw_X = (iter(d.items()) for d in [{"foo": 2, "bar": "a"},
                                       {"baz": "abc", "quux": 4, "foo": -1}])
    h = FeatureLightweightRandomIndexing(n_features=16, input_type="pair")
    x1, x2 = h.transform(raw_X).toarray()
    x1_nz = set(np.abs(x1[x1 != 0]))
    x2_nz = set(np.abs(x2[x2 != 0]))
    assert {1, 2} == x1_nz
    assert {1, 4} == x2_nz

    raw_X = (iter(d.items()) for d in [{"bax": "abc"},
                                       {"bax": "abc"}])
    x1, x2 = h.transform(raw_X).toarray()
    x1_nz = np.abs(x1[x1 != 0])
    assert_array_equal([1,1], x1_nz)
    assert_array_equal(x1, x2)


def test_lri_empty_input():
    n_features = 16
    raw_X = [[], (), iter(range(0))]

    h = FeatureLightweightRandomIndexing(n_features=n_features, input_type="string")
    X = h.transform(raw_X)

    assert_array_equal(X.A, np.zeros((len(raw_X), n_features)))


def test_lri_invalid_input():
    with pytest.raises(ValueError):
        FeatureLightweightRandomIndexing(input_type="gobbledygook")
    with pytest.raises(ValueError):
        FeatureLightweightRandomIndexing(n_features=-1)
    with pytest.raises(ValueError):
        FeatureLightweightRandomIndexing(n_features=0)
    with pytest.raises(TypeError):
        FeatureLightweightRandomIndexing(n_features='ham')

    h = FeatureLightweightRandomIndexing(n_features=np.uint16(2 ** 6))
    with pytest.raises(ValueError):
        h.transform([])
    with pytest.raises(Exception):
        h.transform([[5.5]])
    with pytest.raises(Exception):
        h.transform([[None]])


def test_lri_set_params():
    # Test delayed input validation in fit (useful for grid search).
    hasher = FeatureLightweightRandomIndexing()
    hasher.set_params(n_features=np.inf)
    with pytest.raises(TypeError):
        hasher.fit()


def test_lri_zeros():
    # Assert that no zeros are materialized in the output.
    X = FeatureLightweightRandomIndexing().transform([{'foo': 0}])
    assert X.data.shape == (0,)
