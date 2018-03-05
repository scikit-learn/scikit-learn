# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import pickle

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal

from sklearn.utils.fixes import divide
from sklearn.utils.fixes import MaskedArray
from sklearn.utils.fixes import unique


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_masked_array_obj_dtype_pickleable():
    marr = MaskedArray([1, None, 'a'], dtype=object)

    for mask in (True, False, [0, 1, 0]):
        marr.mask = mask
        marr_pickled = pickle.loads(pickle.dumps(marr))
        assert_array_equal(marr.data, marr_pickled.data)
        assert_array_equal(marr.mask, marr_pickled.mask)


def test_unique():
    ar = []

    # 0-length array, no optional_returns
    u_values = unique([])

    # 0-length array, all optional_returns
    u_values, ind, inv, counts = unique(
        ar, return_index=True, return_inverse=True, return_counts=True)

    ar = [4, 2, 5, 5, 3, 1, 4]

    # Normal array, no optional_returns
    u_values = unique(ar)
    assert_array_equal(u_values, [1, 2, 3, 4, 5])

    # Normal array, all optional_returns
    u_values, ind, inv, counts = unique(
        ar, return_index=True, return_inverse=True, return_counts=True)

    assert_array_equal(u_values, [1, 2, 3, 4, 5])
    assert_array_equal(ind, [5, 1, 4, 0, 2])
    assert_array_equal(inv, [3, 1, 4, 4, 2, 0, 3])
    assert_array_equal(counts, [1, 1, 1, 2, 2])
