# Authors: Lars Buitinck
#          Dan Blanchard <dblanchard@ets.org>
# License: BSD 3 clause

from random import Random
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal

import pytest

from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_selection import SelectKBest, chi2


@pytest.mark.parametrize('sparse', (True, False))
@pytest.mark.parametrize('dtype', (int, np.float32, np.int16))
@pytest.mark.parametrize('sort', (True, False))
@pytest.mark.parametrize('iterable', (True, False))
def test_dictvectorizer(sparse, dtype, sort, iterable):
    D = [{"foo": 1, "bar": 3},
         {"bar": 4, "baz": 2},
         {"bar": 1, "quux": 1, "quuux": 2}]

    v = DictVectorizer(sparse=sparse, dtype=dtype, sort=sort)
    X = v.fit_transform(iter(D) if iterable else D)

    assert sp.issparse(X) == sparse
    assert X.shape == (3, 5)
    assert X.sum() == 14
    assert v.inverse_transform(X) == D

    if sparse:
        # CSR matrices can't be compared for equality
        assert_array_equal(X.A, v.transform(iter(D) if iterable
                                            else D).A)
    else:
        assert_array_equal(X, v.transform(iter(D) if iterable
                                          else D))

    if sort:
        assert (v.feature_names_ ==
                     sorted(v.feature_names_))


def test_feature_selection():
    # make two feature dicts with two useful features and a bunch of useless
    # ones, in terms of chi2
    d1 = dict([("useless%d" % i, 10) for i in range(20)],
              useful1=1, useful2=20)
    d2 = dict([("useless%d" % i, 10) for i in range(20)],
              useful1=20, useful2=1)

    for indices in (True, False):
        v = DictVectorizer().fit([d1, d2])
        X = v.transform([d1, d2])
        sel = SelectKBest(chi2, k=2).fit(X, [0, 1])

        v.restrict(sel.get_support(indices=indices), indices=indices)
        assert v.get_feature_names() == ["useful1", "useful2"]


def test_one_of_k():
    D_in = [{"version": "1", "ham": 2},
            {"version": "2", "spam": .3},
            {"version=3": True, "spam": -1}]
    v = DictVectorizer()
    X = v.fit_transform(D_in)
    assert X.shape == (3, 5)

    D_out = v.inverse_transform(X)
    assert D_out[0] == {"version=1": 1, "ham": 2}

    names = v.get_feature_names()
    assert "version=2" in names
    assert "version" not in names


def test_unseen_or_no_features():
    D = [{"camelot": 0, "spamalot": 1}]
    for sparse in [True, False]:
        v = DictVectorizer(sparse=sparse).fit(D)

        X = v.transform({"push the pram a lot": 2})
        if sparse:
            X = X.toarray()
        assert_array_equal(X, np.zeros((1, 2)))

        X = v.transform({})
        if sparse:
            X = X.toarray()
        assert_array_equal(X, np.zeros((1, 2)))

        try:
            v.transform([])
        except ValueError as e:
            assert "empty" in str(e)


def test_deterministic_vocabulary():
    # Generate equal dictionaries with different memory layouts
    items = [("%03d" % i, i) for i in range(1000)]
    rng = Random(42)
    d_sorted = dict(items)
    rng.shuffle(items)
    d_shuffled = dict(items)

    # check that the memory layout does not impact the resulting vocabulary
    v_1 = DictVectorizer().fit([d_sorted])
    v_2 = DictVectorizer().fit([d_shuffled])

    assert v_1.vocabulary_ == v_2.vocabulary_


def test_n_features_in():
    # For vectorizers, n_features_in_ does not make sense and it is always
    # None
    dv = DictVectorizer()
    assert dv.n_features_in_ is None
    d = [{'foo': 1, 'bar': 2}, {'foo': 3, 'baz': 1}]
    dv.fit(d)
    assert dv.n_features_in_ is None
