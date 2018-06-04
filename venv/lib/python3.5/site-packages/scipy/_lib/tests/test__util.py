from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal, assert_
from pytest import raises as assert_raises

from scipy._lib._util import _aligned_zeros, check_random_state


def test__aligned_zeros():
    niter = 10

    def check(shape, dtype, order, align):
        err_msg = repr((shape, dtype, order, align))
        x = _aligned_zeros(shape, dtype, order, align=align)
        if align is None:
            align = np.dtype(dtype).alignment
        assert_equal(x.__array_interface__['data'][0] % align, 0)
        if hasattr(shape, '__len__'):
            assert_equal(x.shape, shape, err_msg)
        else:
            assert_equal(x.shape, (shape,), err_msg)
        assert_equal(x.dtype, dtype)
        if order == "C":
            assert_(x.flags.c_contiguous, err_msg)
        elif order == "F":
            if x.size > 0:
                # Size-0 arrays get invalid flags on Numpy 1.5
                assert_(x.flags.f_contiguous, err_msg)
        elif order is None:
            assert_(x.flags.c_contiguous, err_msg)
        else:
            raise ValueError()

    # try various alignments
    for align in [1, 2, 3, 4, 8, 16, 32, 64, None]:
        for n in [0, 1, 3, 11]:
            for order in ["C", "F", None]:
                for dtype in [np.uint8, np.float64]:
                    for shape in [n, (1, 2, 3, n)]:
                        for j in range(niter):
                            check(shape, dtype, order, align)


def test_check_random_state():
    # If seed is None, return the RandomState singleton used by np.random.
    # If seed is an int, return a new RandomState instance seeded with seed.
    # If seed is already a RandomState instance, return it.
    # Otherwise raise ValueError.
    rsi = check_random_state(1)
    assert_equal(type(rsi), np.random.RandomState)
    rsi = check_random_state(rsi)
    assert_equal(type(rsi), np.random.RandomState)
    rsi = check_random_state(None)
    assert_equal(type(rsi), np.random.RandomState)
    assert_raises(ValueError, check_random_state, 'a')
