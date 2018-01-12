# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import pickle

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose

from sklearn.utils.fixes import divide
from sklearn.utils.fixes import MaskedArray
from sklearn.utils.fixes import nanpercentile


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_masked_array_obj_dtype_pickleable():
    marr = MaskedArray([1, None, 'a'], dtype=object)

    for mask in (True, False, [0, 1, 0]):
        marr.mask = mask
        marr_pickled = pickle.loads(pickle.dumps(marr))
        assert_array_equal(marr.data, marr_pickled.data)
        assert_array_equal(marr.mask, marr_pickled.mask)


def test_nanpercentile():
    X = np.array([0, 1, 2, np.nan])
    percentile = nanpercentile(X, [0, 50, 100])
    assert_allclose(percentile, X[~np.isnan(X)])
