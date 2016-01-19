# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import numpy as np
import math
from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_true
from numpy.testing import (assert_almost_equal,
                           assert_array_almost_equal)

from sklearn.utils.fixes import divide, expit
from sklearn.utils.fixes import astype
from sklearn.utils.fixes import norm_2d
from sklearn.utils.testing import assert_array_equal


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


def test_norm_2d():
    X = np.array([[-2, 4, 5],
                  [1, 3, -4],
                  [0, 0, 8]]).astype(float)

    # Test various axis and order
    assert_equal(math.sqrt(135), norm_2d(X))
    assert_array_equal(
        np.array([math.sqrt(5), math.sqrt(25), math.sqrt(105)]),
        norm_2d(X, axis=0)
    )
    assert_array_equal(np.array([3, 7, 17]), norm_2d(X, axis=0, ord=1))
    assert_array_equal(np.array([2, 4, 8]), norm_2d(X, axis=0, ord=math.inf))
    assert_array_equal(np.array([0, 0, 4]), norm_2d(X, axis=0, ord=-math.inf))
    assert_array_equal(np.array([11, 8, 8]), norm_2d(X, axis=1, ord=1))

    # Test keepdims
    assert_array_equal(np.array([[11], [8], [8]]), norm_2d(X, axis=1, ord=1, keepdims=True))
    assert_array_equal(np.array([[3, 7, 17]]), norm_2d(X, axis=0, ord=1, keepdims=True))
    assert_equal((1, 1), norm_2d(X, keepdims=True).shape)
