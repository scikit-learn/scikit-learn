# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import pickle
import numpy as np
import math

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal

from sklearn.utils.fixes import divide
from sklearn.utils.fixes import MaskedArray
from sklearn.utils.fixes import norm


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_masked_array_obj_dtype_pickleable():
    marr = MaskedArray([1, None, 'a'], dtype=object)

    for mask in (True, False, [0, 1, 0]):
        marr.mask = mask
        marr_pickled = pickle.loads(pickle.dumps(marr))
        assert_array_equal(marr.data, marr_pickled.data)
        assert_array_equal(marr.mask, marr_pickled.mask)


def test_norm():
    X = np.array([[-2, 4, 5],
                  [1, 3, -4],
                  [0, 0, 8],
                  [0, 0, 0]]).astype(float)

    # Test various axis and order
    assert_equal(math.sqrt(135), norm(X))
    assert_array_equal(
        np.array([math.sqrt(5), math.sqrt(25), math.sqrt(105)]),
        norm(X, axis=0)
    )
    assert_array_equal(np.array([3, 7, 17]), norm(X, axis=0, ord=1))
    assert_array_equal(np.array([2, 4, 8]), norm(X, axis=0, ord=np.inf))
    assert_array_equal(np.array([0, 0, 0]), norm(X, axis=0, ord=-np.inf))
    assert_array_equal(np.array([11, 8, 8, 0]), norm(X, axis=1, ord=1))

    # Test shapes
    assert_equal((), norm(X).shape)
    assert_equal((3,), norm(X, axis=0).shape)
    assert_equal((4,), norm(X, axis=1).shape)
