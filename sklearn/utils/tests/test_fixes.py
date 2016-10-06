# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import pickle
import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal

from sklearn.utils.fixes import divide, expit
from sklearn.utils.fixes import astype
from sklearn.utils.fixes import MaskedArray


def test_expit():
    # Check numerical stability of expit (logistic function).

    # Simulate our previous Cython implementation, based on
    #http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression
    assert_almost_equal(expit(1000.), 1. / (1. + np.exp(-1000.)), decimal=16)
    assert_almost_equal(expit(-1000.), np.exp(-1000.) / (1. + np.exp(-1000.)),
                        decimal=16)

    x = np.arange(10)
    out = np.zeros_like(x, dtype=np.float32)
    assert_array_almost_equal(expit(x), expit(x, out=out))


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_astype_copy_memory():
    a_int32 = np.ones(3, np.int32)

    # Check that dtype conversion works
    b_float32 = astype(a_int32, dtype=np.float32, copy=False)
    assert_equal(b_float32.dtype, np.float32)

    # Changing dtype forces a copy even if copy=False
    assert_false(np.may_share_memory(b_float32, a_int32))

    # Check that copy can be skipped if requested dtype match
    c_int32 = astype(a_int32, dtype=np.int32, copy=False)
    assert_true(c_int32 is a_int32)

    # Check that copy can be forced, and is the case by default:
    d_int32 = astype(a_int32, dtype=np.int32, copy=True)
    assert_false(np.may_share_memory(d_int32, a_int32))

    e_int32 = astype(a_int32, dtype=np.int32)
    assert_false(np.may_share_memory(e_int32, a_int32))


def test_masked_array_obj_dtype_pickleable():
    marr = MaskedArray([1, None, 'a'], dtype=object)

    for mask in (True, False, [0, 1, 0]):
        marr.mask = mask
        marr_pickled = pickle.loads(pickle.dumps(marr))
        assert_array_equal(marr.data, marr_pickled.data)
        assert_array_equal(marr.mask, marr_pickled.mask)
