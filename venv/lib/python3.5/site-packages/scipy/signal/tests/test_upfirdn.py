# Code adapted from "upfirdn" python library with permission:
#
# Copyright (c) 2009, Motorola, Inc
#
# All Rights Reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# * Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# * Neither the name of Motorola nor the names of its contributors may be
# used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import numpy as np
from itertools import product

from numpy.testing import assert_equal, assert_allclose
from pytest import raises as assert_raises

from scipy.signal import upfirdn, firwin, lfilter
from scipy.signal._upfirdn import _output_len


def upfirdn_naive(x, h, up=1, down=1):
    """Naive upfirdn processing in Python

    Note: arg order (x, h) differs to facilitate apply_along_axis use.
    """
    h = np.asarray(h)
    out = np.zeros(len(x) * up, x.dtype)
    out[::up] = x
    out = np.convolve(h, out)[::down][:_output_len(len(h), len(x), up, down)]
    return out


class UpFIRDnCase(object):
    """Test _UpFIRDn object"""
    def __init__(self, up, down, h, x_dtype):
        self.up = up
        self.down = down
        self.h = np.atleast_1d(h)
        self.x_dtype = x_dtype
        self.rng = np.random.RandomState(17)

    def __call__(self):
        # tiny signal
        self.scrub(np.ones(1, self.x_dtype))
        # ones
        self.scrub(np.ones(10, self.x_dtype))  # ones
        # randn
        x = self.rng.randn(10).astype(self.x_dtype)
        if self.x_dtype in (np.complex64, np.complex128):
            x += 1j * self.rng.randn(10)
        self.scrub(x)
        # ramp
        self.scrub(np.arange(10).astype(self.x_dtype))
        # 3D, random
        size = (2, 3, 5)
        x = self.rng.randn(*size).astype(self.x_dtype)
        if self.x_dtype in (np.complex64, np.complex128):
            x += 1j * self.rng.randn(*size)
        for axis in range(len(size)):
            self.scrub(x, axis=axis)
        x = x[:, ::2, 1::3].T
        for axis in range(len(size)):
            self.scrub(x, axis=axis)

    def scrub(self, x, axis=-1):
        yr = np.apply_along_axis(upfirdn_naive, axis, x,
                                 self.h, self.up, self.down)
        y = upfirdn(self.h, x, self.up, self.down, axis=axis)
        dtypes = (self.h.dtype, x.dtype)
        if all(d == np.complex64 for d in dtypes):
            assert_equal(y.dtype, np.complex64)
        elif np.complex64 in dtypes and np.float32 in dtypes:
            assert_equal(y.dtype, np.complex64)
        elif all(d == np.float32 for d in dtypes):
            assert_equal(y.dtype, np.float32)
        elif np.complex128 in dtypes or np.complex64 in dtypes:
            assert_equal(y.dtype, np.complex128)
        else:
            assert_equal(y.dtype, np.float64)
        assert_allclose(yr, y)


class TestUpfirdn(object):

    def test_valid_input(self):
        assert_raises(ValueError, upfirdn, [1], [1], 1, 0)  # up or down < 1
        assert_raises(ValueError, upfirdn, [], [1], 1, 1)  # h.ndim != 1
        assert_raises(ValueError, upfirdn, [[1]], [1], 1, 1)

    def test_vs_lfilter(self):
        # Check that up=1.0 gives same answer as lfilter + slicing
        random_state = np.random.RandomState(17)
        try_types = (int, np.float32, np.complex64, float, complex)
        size = 10000
        down_factors = [2, 11, 79]

        for dtype in try_types:
            x = random_state.randn(size).astype(dtype)
            if dtype in (np.complex64, np.complex128):
                x += 1j * random_state.randn(size)

            for down in down_factors:
                h = firwin(31, 1. / down, window='hamming')
                yl = lfilter(h, 1.0, x)[::down]
                y = upfirdn(h, x, up=1, down=down)
                assert_allclose(yl, y[:yl.size], atol=1e-7, rtol=1e-7)

    def test_vs_naive(self):
        tests = []
        try_types = (int, np.float32, np.complex64, float, complex)

        # Simple combinations of factors
        for x_dtype, h in product(try_types, (1., 1j)):
                tests.append(UpFIRDnCase(1, 1, h, x_dtype))
                tests.append(UpFIRDnCase(2, 2, h, x_dtype))
                tests.append(UpFIRDnCase(3, 2, h, x_dtype))
                tests.append(UpFIRDnCase(2, 3, h, x_dtype))

        # mixture of big, small, and both directions (net up and net down)
        # use all combinations of data and filter dtypes
        factors = (100, 10)  # up/down factors
        cases = product(factors, factors, try_types, try_types)
        for case in cases:
            tests += self._random_factors(*case)

        for test in tests:
            test()

    def _random_factors(self, p_max, q_max, h_dtype, x_dtype):
        n_rep = 3
        longest_h = 25
        random_state = np.random.RandomState(17)
        tests = []

        for _ in range(n_rep):
            # Randomize the up/down factors somewhat
            p_add = q_max if p_max > q_max else 1
            q_add = p_max if q_max > p_max else 1
            p = random_state.randint(p_max) + p_add
            q = random_state.randint(q_max) + q_add

            # Generate random FIR coefficients
            len_h = random_state.randint(longest_h) + 1
            h = np.atleast_1d(random_state.randint(len_h))
            h = h.astype(h_dtype)
            if h_dtype == complex:
                h += 1j * random_state.randint(len_h)

            tests.append(UpFIRDnCase(p, q, h, x_dtype))

        return tests
