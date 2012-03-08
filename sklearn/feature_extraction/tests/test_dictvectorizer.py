# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: BSD-style.

import numpy as np
import scipy.sparse as sp

from nose.tools import assert_equal
from numpy.testing import assert_array_equal

from sklearn.feature_extraction import DictVectorizer


def test_dictvectorizer():
    D = [{"foo": 1, "bar": 3},
         {"bar": 4, "baz": 2},
         {"bar": 1, "quux": 1, "quuux": 2}]

    for sparse in (True, False):
        for dtype in (int, np.float32, np.int16):
            v = DictVectorizer(sparse=sparse, dtype=dtype)
            X = v.fit_transform(D)

            assert_equal(sp.issparse(X), sparse)
            assert_equal(X.shape, (3, 5))
            assert_equal(X.sum(), 14)
            assert_equal(v.inverse_transform(X), D)


def test_unseen_features():
    D = [{"camelot": 0, "spamalot": 1}]
    v = DictVectorizer(sparse=False).fit(D)
    X = v.transform({"push the pram a lot": 2})

    assert_array_equal(X, np.zeros((1, 2)))
