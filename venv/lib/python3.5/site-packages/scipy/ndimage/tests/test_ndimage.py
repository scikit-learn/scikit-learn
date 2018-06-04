# Copyright (C) 2003-2005 Peter J. Verveer
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import division, print_function, absolute_import

import math
import sys

import numpy
from numpy import fft
from numpy.testing import (assert_, assert_equal, assert_array_equal,
                           assert_array_almost_equal, assert_almost_equal)
import pytest
from pytest import raises as assert_raises
from scipy._lib._numpy_compat import suppress_warnings
import scipy.ndimage as ndimage


eps = 1e-12


def sumsq(a, b):
    return math.sqrt(((a - b)**2).sum())


class TestNdimage:
    def setup_method(self):
        # list of numarray data types
        self.integer_types = [
            numpy.int8, numpy.uint8, numpy.int16, numpy.uint16,
            numpy.int32, numpy.uint32, numpy.int64, numpy.uint64]

        self.float_types = [numpy.float32, numpy.float64]

        self.types = self.integer_types + self.float_types

        # list of boundary modes:
        self.modes = ['nearest', 'wrap', 'reflect', 'mirror', 'constant']

    def test_correlate01(self):
        array = numpy.array([1, 2])
        weights = numpy.array([2])
        expected = [2, 4]

        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, expected)

    def test_correlate02(self):
        array = numpy.array([1, 2, 3])
        kernel = numpy.array([1])

        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal(array, output)

        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal(array, output)

        output = ndimage.correlate1d(array, kernel)
        assert_array_almost_equal(array, output)

        output = ndimage.convolve1d(array, kernel)
        assert_array_almost_equal(array, output)

    def test_correlate03(self):
        array = numpy.array([1])
        weights = numpy.array([1, 1])
        expected = [2]

        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, expected)

        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, expected)

    def test_correlate04(self):
        array = numpy.array([1, 2])
        tcor = [2, 3]
        tcov = [3, 4]
        weights = numpy.array([1, 1])
        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, tcov)
        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, tcov)

    def test_correlate05(self):
        array = numpy.array([1, 2, 3])
        tcor = [2, 3, 5]
        tcov = [3, 5, 6]
        kernel = numpy.array([1, 1])
        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal(tcor, output)
        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal(tcov, output)
        output = ndimage.correlate1d(array, kernel)
        assert_array_almost_equal(tcor, output)
        output = ndimage.convolve1d(array, kernel)
        assert_array_almost_equal(tcov, output)

    def test_correlate06(self):
        array = numpy.array([1, 2, 3])
        tcor = [9, 14, 17]
        tcov = [7, 10, 15]
        weights = numpy.array([1, 2, 3])
        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, tcov)
        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, tcov)

    def test_correlate07(self):
        array = numpy.array([1, 2, 3])
        expected = [5, 8, 11]
        weights = numpy.array([1, 2, 1])
        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, expected)
        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, expected)
        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, expected)
        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, expected)

    def test_correlate08(self):
        array = numpy.array([1, 2, 3])
        tcor = [1, 2, 5]
        tcov = [3, 6, 7]
        weights = numpy.array([1, 2, -1])
        output = ndimage.correlate(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve(array, weights)
        assert_array_almost_equal(output, tcov)
        output = ndimage.correlate1d(array, weights)
        assert_array_almost_equal(output, tcor)
        output = ndimage.convolve1d(array, weights)
        assert_array_almost_equal(output, tcov)

    def test_correlate09(self):
        array = []
        kernel = numpy.array([1, 1])
        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal(array, output)
        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal(array, output)
        output = ndimage.correlate1d(array, kernel)
        assert_array_almost_equal(array, output)
        output = ndimage.convolve1d(array, kernel)
        assert_array_almost_equal(array, output)

    def test_correlate10(self):
        array = [[]]
        kernel = numpy.array([[1, 1]])
        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal(array, output)
        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal(array, output)

    def test_correlate11(self):
        array = numpy.array([[1, 2, 3],
                             [4, 5, 6]])
        kernel = numpy.array([[1, 1],
                              [1, 1]])
        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal([[4, 6, 10], [10, 12, 16]], output)
        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal([[12, 16, 18], [18, 22, 24]], output)

    def test_correlate12(self):
        array = numpy.array([[1, 2, 3],
                             [4, 5, 6]])
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        output = ndimage.correlate(array, kernel)
        assert_array_almost_equal([[2, 3, 5], [5, 6, 8]], output)
        output = ndimage.convolve(array, kernel)
        assert_array_almost_equal([[6, 8, 9], [9, 11, 12]], output)

    def test_correlate13(self):
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [4, 5, 6]], type1)
            for type2 in self.types:
                output = ndimage.correlate(array, kernel, output=type2)
                assert_array_almost_equal([[2, 3, 5], [5, 6, 8]], output)
                assert_equal(output.dtype.type, type2)

                output = ndimage.convolve(array, kernel,
                                          output=type2)
                assert_array_almost_equal([[6, 8, 9], [9, 11, 12]], output)
                assert_equal(output.dtype.type, type2)

    def test_correlate14(self):
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [4, 5, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros(array.shape, type2)
                ndimage.correlate(array, kernel,
                                  output=output)
                assert_array_almost_equal([[2, 3, 5], [5, 6, 8]], output)
                assert_equal(output.dtype.type, type2)

                ndimage.convolve(array, kernel, output=output)
                assert_array_almost_equal([[6, 8, 9], [9, 11, 12]], output)
                assert_equal(output.dtype.type, type2)

    def test_correlate15(self):
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [4, 5, 6]], type1)
            output = ndimage.correlate(array, kernel,
                                       output=numpy.float32)
            assert_array_almost_equal([[2, 3, 5], [5, 6, 8]], output)
            assert_equal(output.dtype.type, numpy.float32)

            output = ndimage.convolve(array, kernel,
                                      output=numpy.float32)
            assert_array_almost_equal([[6, 8, 9], [9, 11, 12]], output)
            assert_equal(output.dtype.type, numpy.float32)

    def test_correlate16(self):
        kernel = numpy.array([[0.5, 0],
                              [0, 0.5]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3], [4, 5, 6]], type1)
            output = ndimage.correlate(array, kernel, output=numpy.float32)
            assert_array_almost_equal([[1, 1.5, 2.5], [2.5, 3, 4]], output)
            assert_equal(output.dtype.type, numpy.float32)

            output = ndimage.convolve(array, kernel, output=numpy.float32)
            assert_array_almost_equal([[3, 4, 4.5], [4.5, 5.5, 6]], output)
            assert_equal(output.dtype.type, numpy.float32)

    def test_correlate17(self):
        array = numpy.array([1, 2, 3])
        tcor = [3, 5, 6]
        tcov = [2, 3, 5]
        kernel = numpy.array([1, 1])
        output = ndimage.correlate(array, kernel, origin=-1)
        assert_array_almost_equal(tcor, output)
        output = ndimage.convolve(array, kernel, origin=-1)
        assert_array_almost_equal(tcov, output)
        output = ndimage.correlate1d(array, kernel, origin=-1)
        assert_array_almost_equal(tcor, output)
        output = ndimage.convolve1d(array, kernel, origin=-1)
        assert_array_almost_equal(tcov, output)

    def test_correlate18(self):
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [4, 5, 6]], type1)
            output = ndimage.correlate(array, kernel,
                                       output=numpy.float32,
                                       mode='nearest', origin=-1)
            assert_array_almost_equal([[6, 8, 9], [9, 11, 12]], output)
            assert_equal(output.dtype.type, numpy.float32)

            output = ndimage.convolve(array, kernel,
                                      output=numpy.float32,
                                      mode='nearest', origin=-1)
            assert_array_almost_equal([[2, 3, 5], [5, 6, 8]], output)
            assert_equal(output.dtype.type, numpy.float32)

    def test_correlate19(self):
        kernel = numpy.array([[1, 0],
                              [0, 1]])
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [4, 5, 6]], type1)
            output = ndimage.correlate(array, kernel,
                                       output=numpy.float32,
                                       mode='nearest', origin=[-1, 0])
            assert_array_almost_equal([[5, 6, 8], [8, 9, 11]], output)
            assert_equal(output.dtype.type, numpy.float32)

            output = ndimage.convolve(array, kernel,
                                      output=numpy.float32,
                                      mode='nearest', origin=[-1, 0])
            assert_array_almost_equal([[3, 5, 6], [6, 8, 9]], output)
            assert_equal(output.dtype.type, numpy.float32)

    def test_correlate20(self):
        weights = numpy.array([1, 2, 1])
        expected = [[5, 10, 15], [7, 14, 21]]
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [2, 4, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros((2, 3), type2)
                ndimage.correlate1d(array, weights, axis=0,
                                    output=output)
                assert_array_almost_equal(output, expected)
                ndimage.convolve1d(array, weights, axis=0,
                                   output=output)
                assert_array_almost_equal(output, expected)

    def test_correlate21(self):
        array = numpy.array([[1, 2, 3],
                             [2, 4, 6]])
        expected = [[5, 10, 15], [7, 14, 21]]
        weights = numpy.array([1, 2, 1])
        output = ndimage.correlate1d(array, weights, axis=0)
        assert_array_almost_equal(output, expected)
        output = ndimage.convolve1d(array, weights, axis=0)
        assert_array_almost_equal(output, expected)

    def test_correlate22(self):
        weights = numpy.array([1, 2, 1])
        expected = [[6, 12, 18], [6, 12, 18]]
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [2, 4, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros((2, 3), type2)
                ndimage.correlate1d(array, weights, axis=0,
                                    mode='wrap', output=output)
                assert_array_almost_equal(output, expected)
                ndimage.convolve1d(array, weights, axis=0,
                                   mode='wrap', output=output)
                assert_array_almost_equal(output, expected)

    def test_correlate23(self):
        weights = numpy.array([1, 2, 1])
        expected = [[5, 10, 15], [7, 14, 21]]
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [2, 4, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros((2, 3), type2)
                ndimage.correlate1d(array, weights, axis=0,
                                    mode='nearest', output=output)
                assert_array_almost_equal(output, expected)
                ndimage.convolve1d(array, weights, axis=0,
                                   mode='nearest', output=output)
                assert_array_almost_equal(output, expected)

    def test_correlate24(self):
        weights = numpy.array([1, 2, 1])
        tcor = [[7, 14, 21], [8, 16, 24]]
        tcov = [[4, 8, 12], [5, 10, 15]]
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [2, 4, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros((2, 3), type2)
                ndimage.correlate1d(array, weights, axis=0,
                                    mode='nearest', output=output, origin=-1)
                assert_array_almost_equal(output, tcor)
                ndimage.convolve1d(array, weights, axis=0,
                                   mode='nearest', output=output, origin=-1)
                assert_array_almost_equal(output, tcov)

    def test_correlate25(self):
        weights = numpy.array([1, 2, 1])
        tcor = [[4, 8, 12], [5, 10, 15]]
        tcov = [[7, 14, 21], [8, 16, 24]]
        for type1 in self.types:
            array = numpy.array([[1, 2, 3],
                                 [2, 4, 6]], type1)
            for type2 in self.types:
                output = numpy.zeros((2, 3), type2)
                ndimage.correlate1d(array, weights, axis=0,
                                    mode='nearest', output=output, origin=1)
                assert_array_almost_equal(output, tcor)
                ndimage.convolve1d(array, weights, axis=0,
                                   mode='nearest', output=output, origin=1)
                assert_array_almost_equal(output, tcov)

    def test_gauss01(self):
        input = numpy.array([[1, 2, 3],
                             [2, 4, 6]], numpy.float32)
        output = ndimage.gaussian_filter(input, 0)
        assert_array_almost_equal(output, input)

    def test_gauss02(self):
        input = numpy.array([[1, 2, 3],
                             [2, 4, 6]], numpy.float32)
        output = ndimage.gaussian_filter(input, 1.0)
        assert_equal(input.dtype, output.dtype)
        assert_equal(input.shape, output.shape)

    def test_gauss03(self):
        # single precision data"
        input = numpy.arange(100 * 100).astype(numpy.float32)
        input.shape = (100, 100)
        output = ndimage.gaussian_filter(input, [1.0, 1.0])

        assert_equal(input.dtype, output.dtype)
        assert_equal(input.shape, output.shape)

        # input.sum() is 49995000.0.  With single precision floats, we can't
        # expect more than 8 digits of accuracy, so use decimal=0 in this test.
        assert_almost_equal(output.sum(dtype='d'), input.sum(dtype='d'),
                            decimal=0)
        assert_(sumsq(input, output) > 1.0)

    def test_gauss04(self):
        input = numpy.arange(100 * 100).astype(numpy.float32)
        input.shape = (100, 100)
        otype = numpy.float64
        output = ndimage.gaussian_filter(input, [1.0, 1.0], output=otype)
        assert_equal(output.dtype.type, numpy.float64)
        assert_equal(input.shape, output.shape)
        assert_(sumsq(input, output) > 1.0)

    def test_gauss05(self):
        input = numpy.arange(100 * 100).astype(numpy.float32)
        input.shape = (100, 100)
        otype = numpy.float64
        output = ndimage.gaussian_filter(input, [1.0, 1.0],
                                         order=1, output=otype)
        assert_equal(output.dtype.type, numpy.float64)
        assert_equal(input.shape, output.shape)
        assert_(sumsq(input, output) > 1.0)

    def test_gauss06(self):
        input = numpy.arange(100 * 100).astype(numpy.float32)
        input.shape = (100, 100)
        otype = numpy.float64
        output1 = ndimage.gaussian_filter(input, [1.0, 1.0], output=otype)
        output2 = ndimage.gaussian_filter(input, 1.0, output=otype)
        assert_array_almost_equal(output1, output2)

    def test_prewitt01(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 0)
            t = ndimage.correlate1d(t, [1.0, 1.0, 1.0], 1)
            output = ndimage.prewitt(array, 0)
            assert_array_almost_equal(t, output)

    def test_prewitt02(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 0)
            t = ndimage.correlate1d(t, [1.0, 1.0, 1.0], 1)
            output = numpy.zeros(array.shape, type_)
            ndimage.prewitt(array, 0, output)
            assert_array_almost_equal(t, output)

    def test_prewitt03(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 1)
            t = ndimage.correlate1d(t, [1.0, 1.0, 1.0], 0)
            output = ndimage.prewitt(array, 1)
            assert_array_almost_equal(t, output)

    def test_prewitt04(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.prewitt(array, -1)
            output = ndimage.prewitt(array, 1)
            assert_array_almost_equal(t, output)

    def test_sobel01(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 0)
            t = ndimage.correlate1d(t, [1.0, 2.0, 1.0], 1)
            output = ndimage.sobel(array, 0)
            assert_array_almost_equal(t, output)

    def test_sobel02(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 0)
            t = ndimage.correlate1d(t, [1.0, 2.0, 1.0], 1)
            output = numpy.zeros(array.shape, type_)
            ndimage.sobel(array, 0, output)
            assert_array_almost_equal(t, output)

    def test_sobel03(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.correlate1d(array, [-1.0, 0.0, 1.0], 1)
            t = ndimage.correlate1d(t, [1.0, 2.0, 1.0], 0)
            output = numpy.zeros(array.shape, type_)
            output = ndimage.sobel(array, 1)
            assert_array_almost_equal(t, output)

    def test_sobel04(self):
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            t = ndimage.sobel(array, -1)
            output = ndimage.sobel(array, 1)
            assert_array_almost_equal(t, output)

    def test_laplace01(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.correlate1d(array, [1, -2, 1], 0)
            tmp2 = ndimage.correlate1d(array, [1, -2, 1], 1)
            output = ndimage.laplace(array)
            assert_array_almost_equal(tmp1 + tmp2, output)

    def test_laplace02(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.correlate1d(array, [1, -2, 1], 0)
            tmp2 = ndimage.correlate1d(array, [1, -2, 1], 1)
            output = numpy.zeros(array.shape, type_)
            ndimage.laplace(array, output=output)
            assert_array_almost_equal(tmp1 + tmp2, output)

    def test_gaussian_laplace01(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.gaussian_filter(array, 1.0, [2, 0])
            tmp2 = ndimage.gaussian_filter(array, 1.0, [0, 2])
            output = ndimage.gaussian_laplace(array, 1.0)
            assert_array_almost_equal(tmp1 + tmp2, output)

    def test_gaussian_laplace02(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.gaussian_filter(array, 1.0, [2, 0])
            tmp2 = ndimage.gaussian_filter(array, 1.0, [0, 2])
            output = numpy.zeros(array.shape, type_)
            ndimage.gaussian_laplace(array, 1.0, output)
            assert_array_almost_equal(tmp1 + tmp2, output)

    def test_generic_laplace01(self):
        def derivative2(input, axis, output, mode, cval, a, b):
            sigma = [a, b / 2.0]
            input = numpy.asarray(input)
            order = [0] * input.ndim
            order[axis] = 2
            return ndimage.gaussian_filter(input, sigma, order,
                                           output, mode, cval)
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = numpy.zeros(array.shape, type_)
            tmp = ndimage.generic_laplace(array, derivative2,
                                          extra_arguments=(1.0,),
                                          extra_keywords={'b': 2.0})
            ndimage.gaussian_laplace(array, 1.0, output)
            assert_array_almost_equal(tmp, output)

    def test_gaussian_gradient_magnitude01(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.gaussian_filter(array, 1.0, [1, 0])
            tmp2 = ndimage.gaussian_filter(array, 1.0, [0, 1])
            output = ndimage.gaussian_gradient_magnitude(array, 1.0)
            expected = tmp1 * tmp1 + tmp2 * tmp2
            expected = numpy.sqrt(expected).astype(type_)
            assert_array_almost_equal(expected, output)

    def test_gaussian_gradient_magnitude02(self):
        for type_ in [numpy.int32, numpy.float32, numpy.float64]:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_) * 100
            tmp1 = ndimage.gaussian_filter(array, 1.0, [1, 0])
            tmp2 = ndimage.gaussian_filter(array, 1.0, [0, 1])
            output = numpy.zeros(array.shape, type_)
            ndimage.gaussian_gradient_magnitude(array, 1.0, output)
            expected = tmp1 * tmp1 + tmp2 * tmp2
            expected = numpy.sqrt(expected).astype(type_)
            assert_array_almost_equal(expected, output)

    def test_generic_gradient_magnitude01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [5, 8, 3, 7, 1],
                             [5, 6, 9, 3, 5]], numpy.float64)

        def derivative(input, axis, output, mode, cval, a, b):
            sigma = [a, b / 2.0]
            input = numpy.asarray(input)
            order = [0] * input.ndim
            order[axis] = 1
            return ndimage.gaussian_filter(input, sigma, order,
                                           output, mode, cval)
        tmp1 = ndimage.gaussian_gradient_magnitude(array, 1.0)
        tmp2 = ndimage.generic_gradient_magnitude(
            array, derivative, extra_arguments=(1.0,),
            extra_keywords={'b': 2.0})
        assert_array_almost_equal(tmp1, tmp2)

    def test_uniform01(self):
        array = numpy.array([2, 4, 6])
        size = 2
        output = ndimage.uniform_filter1d(array, size, origin=-1)
        assert_array_almost_equal([3, 5, 6], output)

    def test_uniform02(self):
        array = numpy.array([1, 2, 3])
        filter_shape = [0]
        output = ndimage.uniform_filter(array, filter_shape)
        assert_array_almost_equal(array, output)

    def test_uniform03(self):
        array = numpy.array([1, 2, 3])
        filter_shape = [1]
        output = ndimage.uniform_filter(array, filter_shape)
        assert_array_almost_equal(array, output)

    def test_uniform04(self):
        array = numpy.array([2, 4, 6])
        filter_shape = [2]
        output = ndimage.uniform_filter(array, filter_shape)
        assert_array_almost_equal([2, 3, 5], output)

    def test_uniform05(self):
        array = []
        filter_shape = [1]
        output = ndimage.uniform_filter(array, filter_shape)
        assert_array_almost_equal([], output)

    def test_uniform06(self):
        filter_shape = [2, 2]
        for type1 in self.types:
            array = numpy.array([[4, 8, 12],
                                 [16, 20, 24]], type1)
            for type2 in self.types:
                output = ndimage.uniform_filter(
                    array, filter_shape, output=type2)
                assert_array_almost_equal([[4, 6, 10], [10, 12, 16]], output)
                assert_equal(output.dtype.type, type2)

    def test_minimum_filter01(self):
        array = numpy.array([1, 2, 3, 4, 5])
        filter_shape = numpy.array([2])
        output = ndimage.minimum_filter(array, filter_shape)
        assert_array_almost_equal([1, 1, 2, 3, 4], output)

    def test_minimum_filter02(self):
        array = numpy.array([1, 2, 3, 4, 5])
        filter_shape = numpy.array([3])
        output = ndimage.minimum_filter(array, filter_shape)
        assert_array_almost_equal([1, 1, 2, 3, 4], output)

    def test_minimum_filter03(self):
        array = numpy.array([3, 2, 5, 1, 4])
        filter_shape = numpy.array([2])
        output = ndimage.minimum_filter(array, filter_shape)
        assert_array_almost_equal([3, 2, 2, 1, 1], output)

    def test_minimum_filter04(self):
        array = numpy.array([3, 2, 5, 1, 4])
        filter_shape = numpy.array([3])
        output = ndimage.minimum_filter(array, filter_shape)
        assert_array_almost_equal([2, 2, 1, 1, 1], output)

    def test_minimum_filter05(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        filter_shape = numpy.array([2, 3])
        output = ndimage.minimum_filter(array, filter_shape)
        assert_array_almost_equal([[2, 2, 1, 1, 1],
                                   [2, 2, 1, 1, 1],
                                   [5, 3, 3, 1, 1]], output)

    def test_minimum_filter06(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 1, 1], [1, 1, 1]]
        output = ndimage.minimum_filter(array, footprint=footprint)
        assert_array_almost_equal([[2, 2, 1, 1, 1],
                                   [2, 2, 1, 1, 1],
                                   [5, 3, 3, 1, 1]], output)

    def test_minimum_filter07(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.minimum_filter(array, footprint=footprint)
        assert_array_almost_equal([[2, 2, 1, 1, 1],
                                   [2, 3, 1, 3, 1],
                                   [5, 5, 3, 3, 1]], output)

    def test_minimum_filter08(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.minimum_filter(array, footprint=footprint, origin=-1)
        assert_array_almost_equal([[3, 1, 3, 1, 1],
                                   [5, 3, 3, 1, 1],
                                   [3, 3, 1, 1, 1]], output)

    def test_minimum_filter09(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.minimum_filter(array, footprint=footprint,
                                        origin=[-1, 0])
        assert_array_almost_equal([[2, 3, 1, 3, 1],
                                   [5, 5, 3, 3, 1],
                                   [5, 3, 3, 1, 1]], output)

    def test_maximum_filter01(self):
        array = numpy.array([1, 2, 3, 4, 5])
        filter_shape = numpy.array([2])
        output = ndimage.maximum_filter(array, filter_shape)
        assert_array_almost_equal([1, 2, 3, 4, 5], output)

    def test_maximum_filter02(self):
        array = numpy.array([1, 2, 3, 4, 5])
        filter_shape = numpy.array([3])
        output = ndimage.maximum_filter(array, filter_shape)
        assert_array_almost_equal([2, 3, 4, 5, 5], output)

    def test_maximum_filter03(self):
        array = numpy.array([3, 2, 5, 1, 4])
        filter_shape = numpy.array([2])
        output = ndimage.maximum_filter(array, filter_shape)
        assert_array_almost_equal([3, 3, 5, 5, 4], output)

    def test_maximum_filter04(self):
        array = numpy.array([3, 2, 5, 1, 4])
        filter_shape = numpy.array([3])
        output = ndimage.maximum_filter(array, filter_shape)
        assert_array_almost_equal([3, 5, 5, 5, 4], output)

    def test_maximum_filter05(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        filter_shape = numpy.array([2, 3])
        output = ndimage.maximum_filter(array, filter_shape)
        assert_array_almost_equal([[3, 5, 5, 5, 4],
                                   [7, 9, 9, 9, 5],
                                   [8, 9, 9, 9, 7]], output)

    def test_maximum_filter06(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 1, 1], [1, 1, 1]]
        output = ndimage.maximum_filter(array, footprint=footprint)
        assert_array_almost_equal([[3, 5, 5, 5, 4],
                                   [7, 9, 9, 9, 5],
                                   [8, 9, 9, 9, 7]], output)

    def test_maximum_filter07(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.maximum_filter(array, footprint=footprint)
        assert_array_almost_equal([[3, 5, 5, 5, 4],
                                   [7, 7, 9, 9, 5],
                                   [7, 9, 8, 9, 7]], output)

    def test_maximum_filter08(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.maximum_filter(array, footprint=footprint, origin=-1)
        assert_array_almost_equal([[7, 9, 9, 5, 5],
                                   [9, 8, 9, 7, 5],
                                   [8, 8, 7, 7, 7]], output)

    def test_maximum_filter09(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                                [7, 6, 9, 3, 5],
                                [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.maximum_filter(array, footprint=footprint,
                                        origin=[-1, 0])
        assert_array_almost_equal([[7, 7, 9, 9, 5],
                                   [7, 9, 8, 9, 7],
                                   [8, 8, 8, 7, 7]], output)

    def test_rank01(self):
        array = numpy.array([1, 2, 3, 4, 5])
        output = ndimage.rank_filter(array, 1, size=2)
        assert_array_almost_equal(array, output)
        output = ndimage.percentile_filter(array, 100, size=2)
        assert_array_almost_equal(array, output)
        output = ndimage.median_filter(array, 2)
        assert_array_almost_equal(array, output)

    def test_rank02(self):
        array = numpy.array([1, 2, 3, 4, 5])
        output = ndimage.rank_filter(array, 1, size=[3])
        assert_array_almost_equal(array, output)
        output = ndimage.percentile_filter(array, 50, size=3)
        assert_array_almost_equal(array, output)
        output = ndimage.median_filter(array, (3,))
        assert_array_almost_equal(array, output)

    def test_rank03(self):
        array = numpy.array([3, 2, 5, 1, 4])
        output = ndimage.rank_filter(array, 1, size=[2])
        assert_array_almost_equal([3, 3, 5, 5, 4], output)
        output = ndimage.percentile_filter(array, 100, size=2)
        assert_array_almost_equal([3, 3, 5, 5, 4], output)

    def test_rank04(self):
        array = numpy.array([3, 2, 5, 1, 4])
        expected = [3, 3, 2, 4, 4]
        output = ndimage.rank_filter(array, 1, size=3)
        assert_array_almost_equal(expected, output)
        output = ndimage.percentile_filter(array, 50, size=3)
        assert_array_almost_equal(expected, output)
        output = ndimage.median_filter(array, size=3)
        assert_array_almost_equal(expected, output)

    def test_rank05(self):
        array = numpy.array([3, 2, 5, 1, 4])
        expected = [3, 3, 2, 4, 4]
        output = ndimage.rank_filter(array, -2, size=3)
        assert_array_almost_equal(expected, output)

    def test_rank06(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [5, 8, 3, 7, 1],
                             [5, 6, 9, 3, 5]])
        expected = [[2, 2, 1, 1, 1],
                    [3, 3, 2, 1, 1],
                    [5, 5, 3, 3, 1]]
        output = ndimage.rank_filter(array, 1, size=[2, 3])
        assert_array_almost_equal(expected, output)
        output = ndimage.percentile_filter(array, 17, size=(2, 3))
        assert_array_almost_equal(expected, output)

    def test_rank07(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [5, 8, 3, 7, 1],
                             [5, 6, 9, 3, 5]])
        expected = [[3, 5, 5, 5, 4],
                    [5, 5, 7, 5, 4],
                    [6, 8, 8, 7, 5]]
        output = ndimage.rank_filter(array, -2, size=[2, 3])
        assert_array_almost_equal(expected, output)

    def test_rank08(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [5, 8, 3, 7, 1],
                             [5, 6, 9, 3, 5]])
        expected = [[3, 3, 2, 4, 4],
                    [5, 5, 5, 4, 4],
                    [5, 6, 7, 5, 5]]
        output = ndimage.percentile_filter(array, 50.0, size=(2, 3))
        assert_array_almost_equal(expected, output)
        output = ndimage.rank_filter(array, 3, size=(2, 3))
        assert_array_almost_equal(expected, output)
        output = ndimage.median_filter(array, size=(2, 3))
        assert_array_almost_equal(expected, output)

    def test_rank09(self):
        expected = [[3, 3, 2, 4, 4],
                    [3, 5, 2, 5, 1],
                    [5, 5, 8, 3, 5]]
        footprint = [[1, 0, 1], [0, 1, 0]]
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = ndimage.rank_filter(array, 1, footprint=footprint)
            assert_array_almost_equal(expected, output)
            output = ndimage.percentile_filter(array, 35, footprint=footprint)
            assert_array_almost_equal(expected, output)

    def test_rank10(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        expected = [[2, 2, 1, 1, 1],
                    [2, 3, 1, 3, 1],
                    [5, 5, 3, 3, 1]]
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.rank_filter(array, 0, footprint=footprint)
        assert_array_almost_equal(expected, output)
        output = ndimage.percentile_filter(array, 0.0, footprint=footprint)
        assert_array_almost_equal(expected, output)

    def test_rank11(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        expected = [[3, 5, 5, 5, 4],
                    [7, 7, 9, 9, 5],
                    [7, 9, 8, 9, 7]]
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.rank_filter(array, -1, footprint=footprint)
        assert_array_almost_equal(expected, output)
        output = ndimage.percentile_filter(array, 100.0, footprint=footprint)
        assert_array_almost_equal(expected, output)

    def test_rank12(self):
        expected = [[3, 3, 2, 4, 4],
                    [3, 5, 2, 5, 1],
                    [5, 5, 8, 3, 5]]
        footprint = [[1, 0, 1], [0, 1, 0]]
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = ndimage.rank_filter(array, 1, footprint=footprint)
            assert_array_almost_equal(expected, output)
            output = ndimage.percentile_filter(array, 50.0,
                                               footprint=footprint)
            assert_array_almost_equal(expected, output)
            output = ndimage.median_filter(array, footprint=footprint)
            assert_array_almost_equal(expected, output)

    def test_rank13(self):
        expected = [[5, 2, 5, 1, 1],
                    [5, 8, 3, 5, 5],
                    [6, 6, 5, 5, 5]]
        footprint = [[1, 0, 1], [0, 1, 0]]
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = ndimage.rank_filter(array, 1, footprint=footprint,
                                         origin=-1)
            assert_array_almost_equal(expected, output)

    def test_rank14(self):
        expected = [[3, 5, 2, 5, 1],
                    [5, 5, 8, 3, 5],
                    [5, 6, 6, 5, 5]]
        footprint = [[1, 0, 1], [0, 1, 0]]
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = ndimage.rank_filter(array, 1, footprint=footprint,
                                         origin=[-1, 0])
            assert_array_almost_equal(expected, output)

    def test_rank15(self):
        "rank filter 15"
        expected = [[2, 3, 1, 4, 1],
                    [5, 3, 7, 1, 1],
                    [5, 5, 3, 3, 3]]
        footprint = [[1, 0, 1], [0, 1, 0]]
        for type_ in self.types:
            array = numpy.array([[3, 2, 5, 1, 4],
                                 [5, 8, 3, 7, 1],
                                 [5, 6, 9, 3, 5]], type_)
            output = ndimage.rank_filter(array, 0, footprint=footprint,
                                         origin=[-1, 0])
            assert_array_almost_equal(expected, output)

    def test_generic_filter1d01(self):
        weights = numpy.array([1.1, 2.2, 3.3])

        def _filter_func(input, output, fltr, total):
            fltr = fltr / total
            for ii in range(input.shape[0] - 2):
                output[ii] = input[ii] * fltr[0]
                output[ii] += input[ii + 1] * fltr[1]
                output[ii] += input[ii + 2] * fltr[2]
        for type_ in self.types:
            a = numpy.arange(12, dtype=type_)
            a.shape = (3, 4)
            r1 = ndimage.correlate1d(a, weights / weights.sum(), 0, origin=-1)
            r2 = ndimage.generic_filter1d(
                a, _filter_func, 3, axis=0, origin=-1,
                extra_arguments=(weights,),
                extra_keywords={'total': weights.sum()})
            assert_array_almost_equal(r1, r2)

    def test_generic_filter01(self):
        filter_ = numpy.array([[1.0, 2.0], [3.0, 4.0]])
        footprint = numpy.array([[1, 0], [0, 1]])
        cf = numpy.array([1., 4.])

        def _filter_func(buffer, weights, total=1.0):
            weights = cf / total
            return (buffer * weights).sum()
        for type_ in self.types:
            a = numpy.arange(12, dtype=type_)
            a.shape = (3, 4)
            r1 = ndimage.correlate(a, filter_ * footprint)
            if type_ in self.float_types:
                r1 /= 5
            else:
                r1 //= 5
            r2 = ndimage.generic_filter(
                a, _filter_func, footprint=footprint, extra_arguments=(cf,),
                extra_keywords={'total': cf.sum()})
            assert_array_almost_equal(r1, r2)

    def test_extend01(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([1, 0])
        expected_values = [[1, 1, 2],
                           [3, 1, 2],
                           [1, 1, 2],
                           [2, 1, 2],
                           [0, 1, 2]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate1d(array, weights, 0,
                                         mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend02(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([1, 0, 0, 0, 0, 0, 0, 0])
        expected_values = [[1, 1, 1],
                           [3, 1, 2],
                           [3, 3, 2],
                           [1, 2, 3],
                           [0, 0, 0]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate1d(array, weights, 0,
                                         mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend03(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([0, 0, 1])
        expected_values = [[2, 3, 3],
                           [2, 3, 1],
                           [2, 3, 3],
                           [2, 3, 2],
                           [2, 3, 0]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate1d(array, weights, 0,
                                         mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend04(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 1])
        expected_values = [[3, 3, 3],
                           [2, 3, 1],
                           [2, 1, 1],
                           [1, 2, 3],
                           [0, 0, 0]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate1d(array, weights, 0,
                                         mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend05(self):
        array = numpy.array([[1, 2, 3],
                             [4, 5, 6],
                             [7, 8, 9]])
        weights = numpy.array([[1, 0], [0, 0]])
        expected_values = [[[1, 1, 2], [1, 1, 2], [4, 4, 5]],
                           [[9, 7, 8], [3, 1, 2], [6, 4, 5]],
                           [[1, 1, 2], [1, 1, 2], [4, 4, 5]],
                           [[5, 4, 5], [2, 1, 2], [5, 4, 5]],
                           [[0, 0, 0], [0, 1, 2], [0, 4, 5]]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights,
                                       mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend06(self):
        array = numpy.array([[1, 2, 3],
                             [4, 5, 6],
                             [7, 8, 9]])
        weights = numpy.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
        expected_values = [[[5, 6, 6], [8, 9, 9], [8, 9, 9]],
                           [[5, 6, 4], [8, 9, 7], [2, 3, 1]],
                           [[5, 6, 6], [8, 9, 9], [8, 9, 9]],
                           [[5, 6, 5], [8, 9, 8], [5, 6, 5]],
                           [[5, 6, 0], [8, 9, 0], [0, 0, 0]]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights,
                                       mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend07(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 1])
        expected_values = [[3, 3, 3],
                           [2, 3, 1],
                           [2, 1, 1],
                           [1, 2, 3],
                           [0, 0, 0]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights, mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend08(self):
        array = numpy.array([[1], [2], [3]])
        weights = numpy.array([[0], [0], [0], [0], [0], [0], [0], [0], [1]])
        expected_values = [[[3], [3], [3]],
                           [[2], [3], [1]],
                           [[2], [1], [1]],
                           [[1], [2], [3]],
                           [[0], [0], [0]]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights, mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend09(self):
        array = numpy.array([1, 2, 3])
        weights = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 1])
        expected_values = [[3, 3, 3],
                           [2, 3, 1],
                           [2, 1, 1],
                           [1, 2, 3],
                           [0, 0, 0]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights,
                                       mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_extend10(self):
        array = numpy.array([[1], [2], [3]])
        weights = numpy.array([[0], [0], [0], [0], [0], [0], [0], [0], [1]])
        expected_values = [[[3], [3], [3]],
                           [[2], [3], [1]],
                           [[2], [1], [1]],
                           [[1], [2], [3]],
                           [[0], [0], [0]]]
        for mode, expected_value in zip(self.modes, expected_values):
            output = ndimage.correlate(array, weights,
                                       mode=mode, cval=0)
            assert_array_equal(output, expected_value)

    def test_boundaries(self):
        def shift(x):
            return (x[0] + 0.5,)

        data = numpy.array([1, 2, 3, 4.])
        expected = {'constant': [1.5, 2.5, 3.5, -1, -1, -1, -1],
                    'wrap': [1.5, 2.5, 3.5, 1.5, 2.5, 3.5, 1.5],
                    'mirror': [1.5, 2.5, 3.5, 3.5, 2.5, 1.5, 1.5],
                    'nearest': [1.5, 2.5, 3.5, 4, 4, 4, 4]}

        for mode in expected:
            assert_array_equal(
                expected[mode],
                ndimage.geometric_transform(data, shift, cval=-1, mode=mode,
                                            output_shape=(7,), order=1))

    def test_boundaries2(self):
        def shift(x):
            return (x[0] - 0.9,)

        data = numpy.array([1, 2, 3, 4])
        expected = {'constant': [-1, 1, 2, 3],
                    'wrap': [3, 1, 2, 3],
                    'mirror': [2, 1, 2, 3],
                    'nearest': [1, 1, 2, 3]}

        for mode in expected:
            assert_array_equal(
                expected[mode],
                ndimage.geometric_transform(data, shift, cval=-1, mode=mode,
                                            output_shape=(4,)))

    def test_fourier_gaussian_real01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.float32, numpy.float64], [6, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.rfft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_gaussian(a, [5.0, 2.5], shape[0], 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.irfft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a), 1, decimal=dec)

    def test_fourier_gaussian_complex01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.complex64, numpy.complex128], [6, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.fft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_gaussian(a, [5.0, 2.5], -1, 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.ifft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    def test_fourier_uniform_real01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.float32, numpy.float64], [6, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.rfft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_uniform(a, [5.0, 2.5], shape[0], 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.irfft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a), 1.0, decimal=dec)

    def test_fourier_uniform_complex01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.complex64, numpy.complex128], [6, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.fft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_uniform(a, [5.0, 2.5], -1, 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.ifft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    def test_fourier_shift_real01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.float32, numpy.float64], [4, 11]):
                expected = numpy.arange(shape[0] * shape[1], dtype=type_)
                expected.shape = shape
                a = fft.rfft(expected, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_shift(a, [1, 1], shape[0], 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.irfft(a, shape[0], 0)
                assert_array_almost_equal(a[1:, 1:], expected[:-1, :-1],
                                          decimal=dec)
                assert_array_almost_equal(a.imag, numpy.zeros(shape),
                                          decimal=dec)

    def test_fourier_shift_complex01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.complex64, numpy.complex128], [4, 11]):
                expected = numpy.arange(shape[0] * shape[1], dtype=type_)
                expected.shape = shape
                a = fft.fft(expected, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_shift(a, [1, 1], -1, 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.ifft(a, shape[0], 0)
                assert_array_almost_equal(a.real[1:, 1:], expected[:-1, :-1],
                                          decimal=dec)
                assert_array_almost_equal(a.imag, numpy.zeros(shape),
                                          decimal=dec)

    def test_fourier_ellipsoid_real01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.float32, numpy.float64], [5, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.rfft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_ellipsoid(a, [5.0, 2.5],
                                              shape[0], 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.irfft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a), 1.0, decimal=dec)

    def test_fourier_ellipsoid_complex01(self):
        for shape in [(32, 16), (31, 15)]:
            for type_, dec in zip([numpy.complex64, numpy.complex128],
                                  [5, 14]):
                a = numpy.zeros(shape, type_)
                a[0, 0] = 1.0
                a = fft.fft(a, shape[0], 0)
                a = fft.fft(a, shape[1], 1)
                a = ndimage.fourier_ellipsoid(a, [5.0, 2.5], -1, 0)
                a = fft.ifft(a, shape[1], 1)
                a = fft.ifft(a, shape[0], 0)
                assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    def test_spline01(self):
        for type_ in self.types:
            data = numpy.ones([], type_)
            for order in range(2, 6):
                out = ndimage.spline_filter(data, order=order)
                assert_array_almost_equal(out, 1)

    def test_spline02(self):
        for type_ in self.types:
            data = numpy.array([1], type_)
            for order in range(2, 6):
                out = ndimage.spline_filter(data, order=order)
                assert_array_almost_equal(out, [1])

    def test_spline03(self):
        for type_ in self.types:
            data = numpy.ones([], type_)
            for order in range(2, 6):
                out = ndimage.spline_filter(data, order,
                                            output=type_)
                assert_array_almost_equal(out, 1)

    def test_spline04(self):
        for type_ in self.types:
            data = numpy.ones([4], type_)
            for order in range(2, 6):
                out = ndimage.spline_filter(data, order)
                assert_array_almost_equal(out, [1, 1, 1, 1])

    def test_spline05(self):
        for type_ in self.types:
            data = numpy.ones([4, 4], type_)
            for order in range(2, 6):
                out = ndimage.spline_filter(data, order=order)
                assert_array_almost_equal(out, [[1, 1, 1, 1],
                                                [1, 1, 1, 1],
                                                [1, 1, 1, 1],
                                                [1, 1, 1, 1]])

    def test_geometric_transform01(self):
        data = numpy.array([1])

        def mapping(x):
            return x
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [1])

    def test_geometric_transform02(self):
        data = numpy.ones([4])

        def mapping(x):
            return x
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [1, 1, 1, 1])

    def test_geometric_transform03(self):
        data = numpy.ones([4])

        def mapping(x):
            return (x[0] - 1,)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [0, 1, 1, 1])

    def test_geometric_transform04(self):
        data = numpy.array([4, 1, 3, 2])

        def mapping(x):
            return (x[0] - 1,)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [0, 4, 1, 3])

    def test_geometric_transform05(self):
        data = numpy.array([[1, 1, 1, 1],
                            [1, 1, 1, 1],
                            [1, 1, 1, 1]])

        def mapping(x):
            return (x[0], x[1] - 1)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [[0, 1, 1, 1],
                                            [0, 1, 1, 1],
                                            [0, 1, 1, 1]])

    def test_geometric_transform06(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])

        def mapping(x):
            return (x[0], x[1] - 1)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [[0, 4, 1, 3],
                                            [0, 7, 6, 8],
                                            [0, 3, 5, 3]])

    def test_geometric_transform07(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])

        def mapping(x):
            return (x[0] - 1, x[1])
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [4, 1, 3, 2],
                                            [7, 6, 8, 5]])

    def test_geometric_transform08(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])

        def mapping(x):
            return (x[0] - 1, x[1] - 1)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, data.shape,
                                              order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_geometric_transform10(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])

        def mapping(x):
            return (x[0] - 1, x[1] - 1)
        for order in range(0, 6):
            if (order > 1):
                filtered = ndimage.spline_filter(data, order=order)
            else:
                filtered = data
            out = ndimage.geometric_transform(filtered, mapping, data.shape,
                                              order=order, prefilter=False)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_geometric_transform13(self):
        data = numpy.ones([2], numpy.float64)

        def mapping(x):
            return (x[0] // 2,)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, [4], order=order)
            assert_array_almost_equal(out, [1, 1, 1, 1])

    def test_geometric_transform14(self):
        data = [1, 5, 2, 6, 3, 7, 4, 4]

        def mapping(x):
            return (2 * x[0],)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, [4], order=order)
            assert_array_almost_equal(out, [1, 2, 3, 4])

    def test_geometric_transform15(self):
        data = [1, 2, 3, 4]

        def mapping(x):
            return (x[0] / 2,)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, [8], order=order)
            assert_array_almost_equal(out[::2], [1, 2, 3, 4])

    def test_geometric_transform16(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9.0, 10, 11, 12]]

        def mapping(x):
            return (x[0], x[1] * 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (3, 2),
                                              order=order)
            assert_array_almost_equal(out, [[1, 3], [5, 7], [9, 11]])

    def test_geometric_transform17(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (x[0] * 2, x[1])
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (1, 4),
                                              order=order)
            assert_array_almost_equal(out, [[1, 2, 3, 4]])

    def test_geometric_transform18(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (x[0] * 2, x[1] * 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (1, 2),
                                              order=order)
            assert_array_almost_equal(out, [[1, 3]])

    def test_geometric_transform19(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (x[0], x[1] / 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (3, 8),
                                              order=order)
            assert_array_almost_equal(out[..., ::2], data)

    def test_geometric_transform20(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (x[0] / 2, x[1])
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (6, 4),
                                              order=order)
            assert_array_almost_equal(out[::2, ...], data)

    def test_geometric_transform21(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (x[0] / 2, x[1] / 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (6, 8),
                                              order=order)
            assert_array_almost_equal(out[::2, ::2], data)

    def test_geometric_transform22(self):
        data = numpy.array([[1, 2, 3, 4],
                            [5, 6, 7, 8],
                            [9, 10, 11, 12]], numpy.float64)

        def mapping1(x):
            return (x[0] / 2, x[1] / 2)

        def mapping2(x):
            return (x[0] * 2, x[1] * 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping1,
                                              (6, 8), order=order)
            out = ndimage.geometric_transform(out, mapping2,
                                              (3, 4), order=order)
            assert_array_almost_equal(out, data)

    def test_geometric_transform23(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x):
            return (1, x[0] * 2)
        for order in range(0, 6):
            out = ndimage.geometric_transform(data, mapping, (2,), order=order)
            out = out.astype(numpy.int32)
            assert_array_almost_equal(out, [5, 7])

    def test_geometric_transform24(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]

        def mapping(x, a, b):
            return (a, x[0] * b)
        for order in range(0, 6):
            out = ndimage.geometric_transform(
                data, mapping, (2,), order=order, extra_arguments=(1,),
                extra_keywords={'b': 2})
            assert_array_almost_equal(out, [5, 7])

    def test_geometric_transform_endianness_with_output_parameter(self):
        # geometric transform given output ndarray or dtype with
        # non-native endianness. see issue #4127
        data = numpy.array([1])

        def mapping(x):
            return x

        for out in [data.dtype, data.dtype.newbyteorder(),
                    numpy.empty_like(data),
                    numpy.empty_like(data).astype(data.dtype.newbyteorder())]:
            returned = ndimage.geometric_transform(data, mapping, data.shape,
                                                   output=out)
            result = out if returned is None else returned
            assert_array_almost_equal(result, [1])

    def test_map_coordinates01(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        idx = numpy.indices(data.shape)
        idx -= 1
        for order in range(0, 6):
            out = ndimage.map_coordinates(data, idx, order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_map_coordinates02(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        idx = numpy.indices(data.shape, numpy.float64)
        idx -= 0.5
        for order in range(0, 6):
            out1 = ndimage.shift(data, 0.5, order=order)
            out2 = ndimage.map_coordinates(data, idx, order=order)
            assert_array_almost_equal(out1, out2)

    def test_map_coordinates03(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]], order='F')
        idx = numpy.indices(data.shape) - 1
        out = ndimage.map_coordinates(data, idx)
        assert_array_almost_equal(out, [[0, 0, 0, 0],
                                        [0, 4, 1, 3],
                                        [0, 7, 6, 8]])
        assert_array_almost_equal(out, ndimage.shift(data, (1, 1)))
        idx = numpy.indices(data[::2].shape) - 1
        out = ndimage.map_coordinates(data[::2], idx)
        assert_array_almost_equal(out, [[0, 0, 0, 0],
                                        [0, 4, 1, 3]])
        assert_array_almost_equal(out, ndimage.shift(data[::2], (1, 1)))
        idx = numpy.indices(data[:, ::2].shape) - 1
        out = ndimage.map_coordinates(data[:, ::2], idx)
        assert_array_almost_equal(out, [[0, 0], [0, 4], [0, 7]])
        assert_array_almost_equal(out, ndimage.shift(data[:, ::2], (1, 1)))

    def test_map_coordinates_endianness_with_output_parameter(self):
        # output parameter given as array or dtype with either endianness
        # see issue #4127
        data = numpy.array([[1, 2], [7, 6]])
        expected = numpy.array([[0, 0], [0, 1]])
        idx = numpy.indices(data.shape)
        idx -= 1
        for out in [data.dtype, data.dtype.newbyteorder(), numpy.empty_like(expected),
                    numpy.empty_like(expected).astype(expected.dtype.newbyteorder())]:
            returned = ndimage.map_coordinates(data, idx, output=out)
            result = out if returned is None else returned
            assert_array_almost_equal(result, expected)

    @pytest.mark.skipif('win32' in sys.platform or numpy.intp(0).itemsize < 8,
                        reason="do not run on 32 bit or windows (no sparse memory)")
    def test_map_coordinates_large_data(self):
        # check crash on large data
        try:
            n = 30000
            a = numpy.empty(n**2, dtype=numpy.float32).reshape(n, n)
            # fill the part we might read
            a[n-3:, n-3:] = 0
            ndimage.map_coordinates(a, [[n - 1.5], [n - 1.5]], order=1)
        except MemoryError:
            raise pytest.skip("Not enough memory available")

    def test_affine_transform01(self):
        data = numpy.array([1])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1]], order=order)
            assert_array_almost_equal(out, [1])

    def test_affine_transform02(self):
        data = numpy.ones([4])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1]], order=order)
            assert_array_almost_equal(out, [1, 1, 1, 1])

    def test_affine_transform03(self):
        data = numpy.ones([4])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1]], -1, order=order)
            assert_array_almost_equal(out, [0, 1, 1, 1])

    def test_affine_transform04(self):
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1]], -1, order=order)
            assert_array_almost_equal(out, [0, 4, 1, 3])

    def test_affine_transform05(self):
        data = numpy.array([[1, 1, 1, 1],
                            [1, 1, 1, 1],
                            [1, 1, 1, 1]])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0], [0, 1]],
                                           [0, -1], order=order)
            assert_array_almost_equal(out, [[0, 1, 1, 1],
                                            [0, 1, 1, 1],
                                            [0, 1, 1, 1]])

    def test_affine_transform06(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0], [0, 1]],
                                           [0, -1], order=order)
            assert_array_almost_equal(out, [[0, 4, 1, 3],
                                            [0, 7, 6, 8],
                                            [0, 3, 5, 3]])

    def test_affine_transform07(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0], [0, 1]],
                                           [-1, 0], order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [4, 1, 3, 2],
                                            [7, 6, 8, 5]])

    def test_affine_transform08(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0], [0, 1]],
                                           [-1, -1], order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_affine_transform09(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            if (order > 1):
                filtered = ndimage.spline_filter(data, order=order)
            else:
                filtered = data
            out = ndimage.affine_transform(filtered, [[1, 0], [0, 1]],
                                           [-1, -1], order=order,
                                           prefilter=False)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_affine_transform10(self):
        data = numpy.ones([2], numpy.float64)
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5]], output_shape=(4,),
                                           order=order)
            assert_array_almost_equal(out, [1, 1, 1, 0])

    def test_affine_transform11(self):
        data = [1, 5, 2, 6, 3, 7, 4, 4]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[2]], 0, (4,), order=order)
            assert_array_almost_equal(out, [1, 2, 3, 4])

    def test_affine_transform12(self):
        data = [1, 2, 3, 4]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5]], 0, (8,), order=order)
            assert_array_almost_equal(out[::2], [1, 2, 3, 4])

    def test_affine_transform13(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9.0, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0], [0, 2]], 0, (3, 2),
                                           order=order)
            assert_array_almost_equal(out, [[1, 3], [5, 7], [9, 11]])

    def test_affine_transform14(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[2, 0], [0, 1]], 0, (1, 4),
                                           order=order)
            assert_array_almost_equal(out, [[1, 2, 3, 4]])

    def test_affine_transform15(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[2, 0], [0, 2]], 0, (1, 2),
                                           order=order)
            assert_array_almost_equal(out, [[1, 3]])

    def test_affine_transform16(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[1, 0.0], [0, 0.5]], 0,
                                           (3, 8), order=order)
            assert_array_almost_equal(out[..., ::2], data)

    def test_affine_transform17(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5, 0], [0, 1]], 0,
                                           (6, 4), order=order)
            assert_array_almost_equal(out[::2, ...], data)

    def test_affine_transform18(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5, 0], [0, 0.5]], 0,
                                           (6, 8), order=order)
            assert_array_almost_equal(out[::2, ::2], data)

    def test_affine_transform19(self):
        data = numpy.array([[1, 2, 3, 4],
                            [5, 6, 7, 8],
                            [9, 10, 11, 12]], numpy.float64)
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5, 0], [0, 0.5]], 0,
                                           (6, 8), order=order)
            out = ndimage.affine_transform(out, [[2.0, 0], [0, 2.0]], 0,
                                           (3, 4), order=order)
            assert_array_almost_equal(out, data)

    def test_affine_transform20(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0], [2]], 0, (2,),
                                           order=order)
            assert_array_almost_equal(out, [1, 3])

    def test_affine_transform21(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[2], [0]], 0, (2,),
                                           order=order)
            assert_array_almost_equal(out, [1, 9])

    def test_affine_transform22(self):
        # shift and offset interaction; see issue #1547
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[2]], [-1], (3,),
                                           order=order)
            assert_array_almost_equal(out, [0, 1, 2])

    def test_affine_transform23(self):
        # shift and offset interaction; see issue #1547
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            out = ndimage.affine_transform(data, [[0.5]], [-1], (8,),
                                           order=order)
            assert_array_almost_equal(out[::2], [0, 4, 1, 3])

    def test_affine_transform24(self):
        # consistency between diagonal and non-diagonal case; see issue #1547
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            with suppress_warnings() as sup:
                sup.filter(UserWarning,
                           "The behaviour of affine_transform with a one-dimensional array .* has changed")
                out1 = ndimage.affine_transform(data, [2], -1, order=order)
            out2 = ndimage.affine_transform(data, [[2]], -1, order=order)
            assert_array_almost_equal(out1, out2)

    def test_affine_transform25(self):
        # consistency between diagonal and non-diagonal case; see issue #1547
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            with suppress_warnings() as sup:
                sup.filter(UserWarning,
                           "The behaviour of affine_transform with a one-dimensional array .* has changed")
                out1 = ndimage.affine_transform(data, [0.5], -1, order=order)
            out2 = ndimage.affine_transform(data, [[0.5]], -1, order=order)
            assert_array_almost_equal(out1, out2)

    def test_affine_transform26(self):
        # test homogeneous coordinates
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            if (order > 1):
                filtered = ndimage.spline_filter(data, order=order)
            else:
                filtered = data
            tform_original = numpy.eye(2)
            offset_original = -numpy.ones((2, 1))
            tform_h1 = numpy.hstack((tform_original, offset_original))
            tform_h2 = numpy.vstack((tform_h1, [[0, 0, 1]]))
            out1 = ndimage.affine_transform(filtered, tform_original,
                                            offset_original.ravel(),
                                            order=order, prefilter=False)
            out2 = ndimage.affine_transform(filtered, tform_h1, order=order,
                                            prefilter=False)
            out3 = ndimage.affine_transform(filtered, tform_h2, order=order,
                                            prefilter=False)
            for out in [out1, out2, out3]:
                assert_array_almost_equal(out, [[0, 0, 0, 0],
                                                [0, 4, 1, 3],
                                                [0, 7, 6, 8]])

    def test_affine_transform27(self):
        # test valid homogeneous transformation matrix
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        tform_h1 = numpy.hstack((numpy.eye(2), -numpy.ones((2, 1))))
        tform_h2 = numpy.vstack((tform_h1, [[5, 2, 1]]))
        assert_raises(ValueError, ndimage.affine_transform, data, tform_h2)

    def test_affine_transform_1d_endianness_with_output_parameter(self):
        # 1d affine transform given output ndarray or dtype with
        # either endianness. see issue #7388
        data = numpy.ones((2, 2))
        for out in [numpy.empty_like(data),
                    numpy.empty_like(data).astype(data.dtype.newbyteorder()),
                    data.dtype, data.dtype.newbyteorder()]:
            with suppress_warnings() as sup:
                sup.filter(UserWarning,
                           "The behaviour of affine_transform with a one-dimensional array .* has changed")
                returned = ndimage.affine_transform(data, [1, 1], output=out)
            result = out if returned is None else returned
            assert_array_almost_equal(result, [[1, 1], [1, 1]])

    def test_affine_transform_multi_d_endianness_with_output_parameter(self):
        # affine transform given output ndarray or dtype with either endianness
        # see issue #4127
        data = numpy.array([1])
        for out in [data.dtype, data.dtype.newbyteorder(),
                    numpy.empty_like(data),
                    numpy.empty_like(data).astype(data.dtype.newbyteorder())]:
            returned = ndimage.affine_transform(data, [[1]], output=out)
            result = out if returned is None else returned
            assert_array_almost_equal(result, [1])

    def test_shift01(self):
        data = numpy.array([1])
        for order in range(0, 6):
            out = ndimage.shift(data, [1], order=order)
            assert_array_almost_equal(out, [0])

    def test_shift02(self):
        data = numpy.ones([4])
        for order in range(0, 6):
            out = ndimage.shift(data, [1], order=order)
            assert_array_almost_equal(out, [0, 1, 1, 1])

    def test_shift03(self):
        data = numpy.ones([4])
        for order in range(0, 6):
            out = ndimage.shift(data, -1, order=order)
            assert_array_almost_equal(out, [1, 1, 1, 0])

    def test_shift04(self):
        data = numpy.array([4, 1, 3, 2])
        for order in range(0, 6):
            out = ndimage.shift(data, 1, order=order)
            assert_array_almost_equal(out, [0, 4, 1, 3])

    def test_shift05(self):
        data = numpy.array([[1, 1, 1, 1],
                            [1, 1, 1, 1],
                            [1, 1, 1, 1]])
        for order in range(0, 6):
            out = ndimage.shift(data, [0, 1], order=order)
            assert_array_almost_equal(out, [[0, 1, 1, 1],
                                            [0, 1, 1, 1],
                                            [0, 1, 1, 1]])

    def test_shift06(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.shift(data, [0, 1], order=order)
            assert_array_almost_equal(out, [[0, 4, 1, 3],
                                            [0, 7, 6, 8],
                                            [0, 3, 5, 3]])

    def test_shift07(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.shift(data, [1, 0], order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [4, 1, 3, 2],
                                            [7, 6, 8, 5]])

    def test_shift08(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            out = ndimage.shift(data, [1, 1], order=order)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_shift09(self):
        data = numpy.array([[4, 1, 3, 2],
                            [7, 6, 8, 5],
                            [3, 5, 3, 6]])
        for order in range(0, 6):
            if (order > 1):
                filtered = ndimage.spline_filter(data, order=order)
            else:
                filtered = data
            out = ndimage.shift(filtered, [1, 1], order=order, prefilter=False)
            assert_array_almost_equal(out, [[0, 0, 0, 0],
                                            [0, 4, 1, 3],
                                            [0, 7, 6, 8]])

    def test_zoom1(self):
        for order in range(0, 6):
            for z in [2, [2, 2]]:
                arr = numpy.array(list(range(25))).reshape((5, 5)).astype(float)
                arr = ndimage.zoom(arr, z, order=order)
                assert_equal(arr.shape, (10, 10))
                assert_(numpy.all(arr[-1, :] != 0))
                assert_(numpy.all(arr[-1, :] >= (20 - eps)))
                assert_(numpy.all(arr[0, :] <= (5 + eps)))
                assert_(numpy.all(arr >= (0 - eps)))
                assert_(numpy.all(arr <= (24 + eps)))

    def test_zoom2(self):
        arr = numpy.arange(12).reshape((3, 4))
        out = ndimage.zoom(ndimage.zoom(arr, 2), 0.5)
        assert_array_equal(out, arr)

    def test_zoom3(self):
        arr = numpy.array([[1, 2]])
        out1 = ndimage.zoom(arr, (2, 1))
        out2 = ndimage.zoom(arr, (1, 2))

        assert_array_almost_equal(out1, numpy.array([[1, 2], [1, 2]]))
        assert_array_almost_equal(out2, numpy.array([[1, 1, 2, 2]]))

    def test_zoom_affine01(self):
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12]]
        for order in range(0, 6):
            with suppress_warnings() as sup:
                sup.filter(UserWarning,
                           "The behaviour of affine_transform with a one-dimensional array .* has changed")
                out = ndimage.affine_transform(data, [0.5, 0.5], 0,
                                               (6, 8), order=order)
            assert_array_almost_equal(out[::2, ::2], data)

    def test_zoom_infinity(self):
        # Ticket #1419 regression test
        dim = 8
        ndimage.zoom(numpy.zeros((dim, dim)), 1./dim, mode='nearest')

    def test_zoom_zoomfactor_one(self):
        # Ticket #1122 regression test
        arr = numpy.zeros((1, 5, 5))
        zoom = (1.0, 2.0, 2.0)

        out = ndimage.zoom(arr, zoom, cval=7)
        ref = numpy.zeros((1, 10, 10))
        assert_array_almost_equal(out, ref)

    def test_zoom_output_shape_roundoff(self):
        arr = numpy.zeros((3, 11, 25))
        zoom = (4.0 / 3, 15.0 / 11, 29.0 / 25)
        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "From scipy 0.13.0, the output shape of zoom.. is calculated with round.. instead of int")
            out = ndimage.zoom(arr, zoom)
        assert_array_equal(out.shape, (4, 15, 29))

    def test_rotate01(self):
        data = numpy.array([[0, 0, 0, 0],
                            [0, 1, 1, 0],
                            [0, 0, 0, 0]], dtype=numpy.float64)
        for order in range(0, 6):
            out = ndimage.rotate(data, 0)
            assert_array_almost_equal(out, data)

    def test_rotate02(self):
        data = numpy.array([[0, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 0, 0]], dtype=numpy.float64)
        expected = numpy.array([[0, 0, 0],
                               [0, 0, 0],
                               [0, 1, 0],
                               [0, 0, 0]], dtype=numpy.float64)
        for order in range(0, 6):
            out = ndimage.rotate(data, 90)
            assert_array_almost_equal(out, expected)

    def test_rotate03(self):
        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0]], dtype=numpy.float64)
        expected = numpy.array([[0, 0, 0],
                               [0, 0, 0],
                               [0, 1, 0],
                               [0, 1, 0],
                               [0, 0, 0]], dtype=numpy.float64)
        for order in range(0, 6):
            out = ndimage.rotate(data, 90)
            assert_array_almost_equal(out, expected)

    def test_rotate04(self):
        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0]], dtype=numpy.float64)
        expected = numpy.array([[0, 0, 0, 0, 0],
                                [0, 0, 1, 0, 0],
                                [0, 0, 1, 0, 0]], dtype=numpy.float64)
        for order in range(0, 6):
            out = ndimage.rotate(data, 90, reshape=False)
            assert_array_almost_equal(out, expected)

    def test_rotate05(self):
        data = numpy.empty((4, 3, 3))
        for i in range(3):
            data[:, :, i] = numpy.array([[0, 0, 0],
                                         [0, 1, 0],
                                         [0, 1, 0],
                                         [0, 0, 0]], dtype=numpy.float64)

        expected = numpy.array([[0, 0, 0, 0],
                                [0, 1, 1, 0],
                                [0, 0, 0, 0]], dtype=numpy.float64)

        for order in range(0, 6):
            out = ndimage.rotate(data, 90)
            for i in range(3):
                assert_array_almost_equal(out[:, :, i], expected)

    def test_rotate06(self):
        data = numpy.empty((3, 4, 3))
        for i in range(3):
            data[:, :, i] = numpy.array([[0, 0, 0, 0],
                                         [0, 1, 1, 0],
                                         [0, 0, 0, 0]], dtype=numpy.float64)

        expected = numpy.array([[0, 0, 0],
                                [0, 1, 0],
                                [0, 1, 0],
                                [0, 0, 0]], dtype=numpy.float64)

        for order in range(0, 6):
            out = ndimage.rotate(data, 90)
            for i in range(3):
                assert_array_almost_equal(out[:, :, i], expected)

    def test_rotate07(self):
        data = numpy.array([[[0, 0, 0, 0, 0],
                             [0, 1, 1, 0, 0],
                             [0, 0, 0, 0, 0]]] * 2, dtype=numpy.float64)
        data = data.transpose()
        expected = numpy.array([[[0, 0, 0],
                                 [0, 1, 0],
                                 [0, 1, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]]] * 2, dtype=numpy.float64)
        expected = expected.transpose([2, 1, 0])

        for order in range(0, 6):
            out = ndimage.rotate(data, 90, axes=(0, 1))
            assert_array_almost_equal(out, expected)

    def test_rotate08(self):
        data = numpy.array([[[0, 0, 0, 0, 0],
                             [0, 1, 1, 0, 0],
                             [0, 0, 0, 0, 0]]] * 2, dtype=numpy.float64)
        data = data.transpose()
        expected = numpy.array([[[0, 0, 1, 0, 0],
                                 [0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 0]]] * 2, dtype=numpy.float64)
        expected = expected.transpose()
        for order in range(0, 6):
            out = ndimage.rotate(data, 90, axes=(0, 1), reshape=False)
            assert_array_almost_equal(out, expected)

    def test_watershed_ift01(self):
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[-1, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 1, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]], numpy.int8)
        out = ndimage.watershed_ift(data, markers, structure=[[1, 1, 1],
                                                              [1, 1, 1],
                                                              [1, 1, 1]])
        expected = [[-1, -1, -1, -1, -1, -1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift02(self):
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[-1, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 1, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]], numpy.int8)
        out = ndimage.watershed_ift(data, markers)
        expected = [[-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, 1, 1, 1, -1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, -1, 1, 1, 1, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift03(self):
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 2, 0, 3, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, -1]], numpy.int8)
        out = ndimage.watershed_ift(data, markers)
        expected = [[-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, 2, -1, 3, -1, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, -1, 2, -1, 3, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift04(self):
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 2, 0, 3, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, -1]],
                              numpy.int8)
        out = ndimage.watershed_ift(data, markers,
                                    structure=[[1, 1, 1],
                                               [1, 1, 1],
                                               [1, 1, 1]])
        expected = [[-1, -1, -1, -1, -1, -1, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, 2, 2, 3, 3, 3, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift05(self):
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 0, 1, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 3, 0, 2, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, -1]],
                              numpy.int8)
        out = ndimage.watershed_ift(data, markers,
                                    structure=[[1, 1, 1],
                                               [1, 1, 1],
                                               [1, 1, 1]])
        expected = [[-1, -1, -1, -1, -1, -1, -1],
                    [-1, 3, 3, 2, 2, 2, -1],
                    [-1, 3, 3, 2, 2, 2, -1],
                    [-1, 3, 3, 2, 2, 2, -1],
                    [-1, 3, 3, 2, 2, 2, -1],
                    [-1, 3, 3, 2, 2, 2, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift06(self):
        data = numpy.array([[0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 0, 0, 0, 1, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[-1, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 1, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]], numpy.int8)
        out = ndimage.watershed_ift(data, markers,
                                    structure=[[1, 1, 1],
                                               [1, 1, 1],
                                               [1, 1, 1]])
        expected = [[-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_watershed_ift07(self):
        shape = (7, 6)
        data = numpy.zeros(shape, dtype=numpy.uint8)
        data = data.transpose()
        data[...] = numpy.array([[0, 1, 0, 0, 0, 1, 0],
                                 [0, 1, 0, 0, 0, 1, 0],
                                 [0, 1, 0, 0, 0, 1, 0],
                                 [0, 1, 1, 1, 1, 1, 0],
                                 [0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0]], numpy.uint8)
        markers = numpy.array([[-1, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 1, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]], numpy.int8)
        out = numpy.zeros(shape, dtype=numpy.int16)
        out = out.transpose()
        ndimage.watershed_ift(data, markers,
                              structure=[[1, 1, 1],
                                         [1, 1, 1],
                                         [1, 1, 1]],
                              output=out)
        expected = [[-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, 1, 1, 1, 1, 1, -1],
                    [-1, -1, -1, -1, -1, -1, -1],
                    [-1, -1, -1, -1, -1, -1, -1]]
        assert_array_almost_equal(out, expected)

    def test_distance_transform_bf01(self):
        # brute force (bf) distance transform
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_bf(data, 'euclidean',
                                                return_indices=True)
        expected = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 1, 2, 4, 2, 1, 0, 0],
                    [0, 0, 1, 4, 8, 4, 1, 0, 0],
                    [0, 0, 1, 2, 4, 2, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        assert_array_almost_equal(out * out, expected)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2, 1, 2, 2, 2, 2],
                     [3, 3, 3, 2, 1, 2, 3, 3, 3],
                     [4, 4, 4, 4, 6, 4, 4, 4, 4],
                     [5, 5, 6, 6, 7, 6, 6, 5, 5],
                     [6, 6, 6, 7, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 1, 2, 4, 6, 7, 7, 8],
                     [0, 1, 1, 1, 6, 7, 7, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_bf02(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_bf(data, 'cityblock',
                                                return_indices=True)

        expected = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 1, 2, 2, 2, 1, 0, 0],
                    [0, 0, 1, 2, 3, 2, 1, 0, 0],
                    [0, 0, 1, 2, 2, 2, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        assert_array_almost_equal(out, expected)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2, 1, 2, 2, 2, 2],
                     [3, 3, 3, 3, 1, 3, 3, 3, 3],
                     [4, 4, 4, 4, 7, 4, 4, 4, 4],
                     [5, 5, 6, 7, 7, 7, 6, 5, 5],
                     [6, 6, 6, 7, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 1, 1, 4, 7, 7, 7, 8],
                     [0, 1, 1, 1, 4, 7, 7, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(expected, ft)

    def test_distance_transform_bf03(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_bf(data, 'chessboard',
                                                return_indices=True)

        expected = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 1, 1, 2, 1, 1, 0, 0],
                    [0, 0, 1, 2, 2, 2, 1, 0, 0],
                    [0, 0, 1, 1, 2, 1, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        assert_array_almost_equal(out, expected)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2, 1, 2, 2, 2, 2],
                     [3, 3, 4, 2, 2, 2, 4, 3, 3],
                     [4, 4, 5, 6, 6, 6, 5, 4, 4],
                     [5, 5, 6, 6, 7, 6, 6, 5, 5],
                     [6, 6, 6, 7, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 5, 6, 6, 7, 8],
                     [0, 1, 1, 2, 6, 6, 7, 7, 8],
                     [0, 1, 1, 2, 6, 7, 7, 7, 8],
                     [0, 1, 2, 2, 6, 6, 7, 7, 8],
                     [0, 1, 2, 4, 5, 6, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_bf04(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        tdt, tft = ndimage.distance_transform_bf(data, return_indices=1)
        dts = []
        fts = []
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ndimage.distance_transform_bf(data, distances=dt)
        dts.append(dt)
        ft = ndimage.distance_transform_bf(
            data, return_distances=False, return_indices=1)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_bf(
            data, return_distances=False, return_indices=True, indices=ft)
        fts.append(ft)
        dt, ft = ndimage.distance_transform_bf(
            data, return_indices=1)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ft = ndimage.distance_transform_bf(
            data, distances=dt, return_indices=True)
        dts.append(dt)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        dt = ndimage.distance_transform_bf(
            data, return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_bf(
            data, distances=dt, return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        for dt in dts:
            assert_array_almost_equal(tdt, dt)
        for ft in fts:
            assert_array_almost_equal(tft, ft)

    def test_distance_transform_bf05(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_bf(
            data, 'euclidean', return_indices=True, sampling=[2, 2])
        expected = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 4, 4, 4, 0, 0, 0],
                    [0, 0, 4, 8, 16, 8, 4, 0, 0],
                    [0, 0, 4, 16, 32, 16, 4, 0, 0],
                    [0, 0, 4, 8, 16, 8, 4, 0, 0],
                    [0, 0, 0, 4, 4, 4, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        assert_array_almost_equal(out * out, expected)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2, 1, 2, 2, 2, 2],
                     [3, 3, 3, 2, 1, 2, 3, 3, 3],
                     [4, 4, 4, 4, 6, 4, 4, 4, 4],
                     [5, 5, 6, 6, 7, 6, 6, 5, 5],
                     [6, 6, 6, 7, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 1, 2, 4, 6, 7, 7, 8],
                     [0, 1, 1, 1, 6, 7, 7, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_bf06(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_bf(
            data, 'euclidean', return_indices=True, sampling=[2, 1])
        expected = [[0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 4, 1, 0, 0, 0],
                    [0, 0, 1, 4, 8, 4, 1, 0, 0],
                    [0, 0, 1, 4, 9, 4, 1, 0, 0],
                    [0, 0, 1, 4, 8, 4, 1, 0, 0],
                    [0, 0, 0, 1, 4, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0]]
        assert_array_almost_equal(out * out, expected)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 2, 2, 2, 2, 2, 2],
                     [3, 3, 3, 3, 2, 3, 3, 3, 3],
                     [4, 4, 4, 4, 4, 4, 4, 4, 4],
                     [5, 5, 5, 5, 6, 5, 5, 5, 5],
                     [6, 6, 6, 6, 7, 6, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 6, 6, 6, 7, 8],
                     [0, 1, 1, 1, 6, 7, 7, 7, 8],
                     [0, 1, 1, 1, 7, 7, 7, 7, 8],
                     [0, 1, 1, 1, 6, 7, 7, 7, 8],
                     [0, 1, 2, 2, 4, 6, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_cdt01(self):
        # chamfer type distance (cdt) transform
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_cdt(
            data, 'cityblock', return_indices=True)
        bf = ndimage.distance_transform_bf(data, 'cityblock')
        assert_array_almost_equal(bf, out)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 1, 1, 1, 2, 2, 2],
                     [3, 3, 2, 1, 1, 1, 2, 3, 3],
                     [4, 4, 4, 4, 1, 4, 4, 4, 4],
                     [5, 5, 5, 5, 7, 7, 6, 5, 5],
                     [6, 6, 6, 6, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 1, 1, 4, 7, 7, 7, 8],
                     [0, 1, 1, 1, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_cdt02(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_cdt(data, 'chessboard',
                                                 return_indices=True)
        bf = ndimage.distance_transform_bf(data, 'chessboard')
        assert_array_almost_equal(bf, out)

        expected = [[[0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [2, 2, 2, 1, 1, 1, 2, 2, 2],
                     [3, 3, 2, 2, 1, 2, 2, 3, 3],
                     [4, 4, 3, 2, 2, 2, 3, 4, 4],
                     [5, 5, 4, 6, 7, 6, 4, 5, 5],
                     [6, 6, 6, 6, 7, 7, 6, 6, 6],
                     [7, 7, 7, 7, 7, 7, 7, 7, 7],
                     [8, 8, 8, 8, 8, 8, 8, 8, 8]],
                    [[0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 2, 3, 4, 6, 7, 8],
                     [0, 1, 1, 2, 2, 6, 6, 7, 8],
                     [0, 1, 1, 1, 2, 6, 7, 7, 8],
                     [0, 1, 1, 2, 6, 6, 7, 7, 8],
                     [0, 1, 2, 2, 5, 6, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8],
                     [0, 1, 2, 3, 4, 5, 6, 7, 8]]]
        assert_array_almost_equal(ft, expected)

    def test_distance_transform_cdt03(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        tdt, tft = ndimage.distance_transform_cdt(data, return_indices=True)
        dts = []
        fts = []
        dt = numpy.zeros(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_cdt(data, distances=dt)
        dts.append(dt)
        ft = ndimage.distance_transform_cdt(
            data, return_distances=False, return_indices=True)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_cdt(
            data, return_distances=False, return_indices=True, indices=ft)
        fts.append(ft)
        dt, ft = ndimage.distance_transform_cdt(
            data, return_indices=True)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.int32)
        ft = ndimage.distance_transform_cdt(
            data, distances=dt, return_indices=True)
        dts.append(dt)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        dt = ndimage.distance_transform_cdt(
            data, return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.int32)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_cdt(data, distances=dt,
                                       return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        for dt in dts:
            assert_array_almost_equal(tdt, dt)
        for ft in fts:
            assert_array_almost_equal(tft, ft)

    def test_distance_transform_edt01(self):
        # euclidean distance transform (edt)
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        out, ft = ndimage.distance_transform_edt(data, return_indices=True)
        bf = ndimage.distance_transform_bf(data, 'euclidean')
        assert_array_almost_equal(bf, out)

        dt = ft - numpy.indices(ft.shape[1:], dtype=ft.dtype)
        dt = dt.astype(numpy.float64)
        numpy.multiply(dt, dt, dt)
        dt = numpy.add.reduce(dt, axis=0)
        numpy.sqrt(dt, dt)

        assert_array_almost_equal(bf, dt)

    def test_distance_transform_edt02(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        tdt, tft = ndimage.distance_transform_edt(data, return_indices=True)
        dts = []
        fts = []
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ndimage.distance_transform_edt(data, distances=dt)
        dts.append(dt)
        ft = ndimage.distance_transform_edt(
            data, return_distances=0, return_indices=True)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_edt(
            data, return_distances=False, return_indices=True, indices=ft)
        fts.append(ft)
        dt, ft = ndimage.distance_transform_edt(
            data, return_indices=True)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ft = ndimage.distance_transform_edt(
            data, distances=dt, return_indices=True)
        dts.append(dt)
        fts.append(ft)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        dt = ndimage.distance_transform_edt(
            data, return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        dt = numpy.zeros(data.shape, dtype=numpy.float64)
        ft = numpy.indices(data.shape, dtype=numpy.int32)
        ndimage.distance_transform_edt(
            data, distances=dt, return_indices=True, indices=ft)
        dts.append(dt)
        fts.append(ft)
        for dt in dts:
            assert_array_almost_equal(tdt, dt)
        for ft in fts:
            assert_array_almost_equal(tft, ft)

    def test_distance_transform_edt03(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        ref = ndimage.distance_transform_bf(data, 'euclidean', sampling=[2, 2])
        out = ndimage.distance_transform_edt(data, sampling=[2, 2])
        assert_array_almost_equal(ref, out)

    def test_distance_transform_edt4(self):
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0]], type_)
        ref = ndimage.distance_transform_bf(data, 'euclidean', sampling=[2, 1])
        out = ndimage.distance_transform_edt(data, sampling=[2, 1])
        assert_array_almost_equal(ref, out)

    def test_distance_transform_edt5(self):
        # Ticket #954 regression test
        out = ndimage.distance_transform_edt(False)
        assert_array_almost_equal(out, [0.])

    def test_generate_structure01(self):
        struct = ndimage.generate_binary_structure(0, 1)
        assert_array_almost_equal(struct, 1)

    def test_generate_structure02(self):
        struct = ndimage.generate_binary_structure(1, 1)
        assert_array_almost_equal(struct, [1, 1, 1])

    def test_generate_structure03(self):
        struct = ndimage.generate_binary_structure(2, 1)
        assert_array_almost_equal(struct, [[0, 1, 0],
                                           [1, 1, 1],
                                           [0, 1, 0]])

    def test_generate_structure04(self):
        struct = ndimage.generate_binary_structure(2, 2)
        assert_array_almost_equal(struct, [[1, 1, 1],
                                           [1, 1, 1],
                                           [1, 1, 1]])

    def test_iterate_structure01(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        out = ndimage.iterate_structure(struct, 2)
        assert_array_almost_equal(out, [[0, 0, 1, 0, 0],
                                        [0, 1, 1, 1, 0],
                                        [1, 1, 1, 1, 1],
                                        [0, 1, 1, 1, 0],
                                        [0, 0, 1, 0, 0]])

    def test_iterate_structure02(self):
        struct = [[0, 1],
                  [1, 1],
                  [0, 1]]
        out = ndimage.iterate_structure(struct, 2)
        assert_array_almost_equal(out, [[0, 0, 1],
                                        [0, 1, 1],
                                        [1, 1, 1],
                                        [0, 1, 1],
                                        [0, 0, 1]])

    def test_iterate_structure03(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        out = ndimage.iterate_structure(struct, 2, 1)
        expected = [[0, 0, 1, 0, 0],
                    [0, 1, 1, 1, 0],
                    [1, 1, 1, 1, 1],
                    [0, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0]]
        assert_array_almost_equal(out[0], expected)
        assert_equal(out[1], [2, 2])

    def test_binary_erosion01(self):
        for type_ in self.types:
            data = numpy.ones([], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, 1)

    def test_binary_erosion02(self):
        for type_ in self.types:
            data = numpy.ones([], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, 1)

    def test_binary_erosion03(self):
        for type_ in self.types:
            data = numpy.ones([1], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [0])

    def test_binary_erosion04(self):
        for type_ in self.types:
            data = numpy.ones([1], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [1])

    def test_binary_erosion05(self):
        for type_ in self.types:
            data = numpy.ones([3], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [0, 1, 0])

    def test_binary_erosion06(self):
        for type_ in self.types:
            data = numpy.ones([3], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [1, 1, 1])

    def test_binary_erosion07(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [0, 1, 1, 1, 0])

    def test_binary_erosion08(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [1, 1, 1, 1, 1])

    def test_binary_erosion09(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [0, 0, 0, 0, 0])

    def test_binary_erosion10(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [1, 0, 0, 0, 1])

    def test_binary_erosion11(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            struct = [1, 0, 1]
            out = ndimage.binary_erosion(data, struct, border_value=1)
            assert_array_almost_equal(out, [1, 0, 1, 0, 1])

    def test_binary_erosion12(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            struct = [1, 0, 1]
            out = ndimage.binary_erosion(data, struct, border_value=1,
                                         origin=-1)
            assert_array_almost_equal(out, [0, 1, 0, 1, 1])

    def test_binary_erosion13(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            struct = [1, 0, 1]
            out = ndimage.binary_erosion(data, struct, border_value=1,
                                         origin=1)
            assert_array_almost_equal(out, [1, 1, 0, 1, 0])

    def test_binary_erosion14(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            struct = [1, 1]
            out = ndimage.binary_erosion(data, struct, border_value=1)
            assert_array_almost_equal(out, [1, 1, 0, 0, 1])

    def test_binary_erosion15(self):
        for type_ in self.types:
            data = numpy.ones([5], type_)
            data[2] = 0
            struct = [1, 1]
            out = ndimage.binary_erosion(data, struct, border_value=1,
                                         origin=-1)
            assert_array_almost_equal(out, [1, 0, 0, 1, 1])

    def test_binary_erosion16(self):
        for type_ in self.types:
            data = numpy.ones([1, 1], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [[1]])

    def test_binary_erosion17(self):
        for type_ in self.types:
            data = numpy.ones([1, 1], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [[0]])

    def test_binary_erosion18(self):
        for type_ in self.types:
            data = numpy.ones([1, 3], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [[0, 0, 0]])

    def test_binary_erosion19(self):
        for type_ in self.types:
            data = numpy.ones([1, 3], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [[1, 1, 1]])

    def test_binary_erosion20(self):
        for type_ in self.types:
            data = numpy.ones([3, 3], type_)
            out = ndimage.binary_erosion(data)
            assert_array_almost_equal(out, [[0, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 0]])

    def test_binary_erosion21(self):
        for type_ in self.types:
            data = numpy.ones([3, 3], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, [[1, 1, 1],
                                            [1, 1, 1],
                                            [1, 1, 1]])

    def test_binary_erosion22(self):
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_erosion(data, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_erosion23(self):
        struct = ndimage.generate_binary_structure(2, 2)
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_erosion(data, struct, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_erosion24(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 1, 1],
                    [0, 0, 0, 1, 1, 1, 0, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 1, 1],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_erosion(data, struct, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_erosion25(self):
        struct = [[0, 1, 0],
                  [1, 0, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 1, 0, 1, 1],
                                [0, 0, 1, 0, 1, 1, 0, 0],
                                [0, 1, 0, 1, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_erosion(data, struct, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_erosion26(self):
        struct = [[0, 1, 0],
                  [1, 0, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0, 0, 1],
                    [0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1]]
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 1, 1],
                                [0, 0, 1, 1, 1, 0, 1, 1],
                                [0, 0, 1, 0, 1, 1, 0, 0],
                                [0, 1, 0, 1, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_erosion(data, struct, border_value=1,
                                         origin=(-1, -1))
            assert_array_almost_equal(out, expected)

    def test_binary_erosion27(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_erosion(data, struct, border_value=1,
                                     iterations=2)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion28(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], bool)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_erosion(data, struct, border_value=1,
                               iterations=2, output=out)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion29(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [1, 1, 1, 1, 1, 1, 1],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0]], bool)
        out = ndimage.binary_erosion(data, struct,
                                     border_value=1, iterations=3)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion30(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [1, 1, 1, 1, 1, 1, 1],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0]], bool)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_erosion(data, struct, border_value=1,
                               iterations=3, output=out)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion31(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 1, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 0, 1],
                    [0, 1, 1, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 1]]
        data = numpy.array([[0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [1, 1, 1, 1, 1, 1, 1],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0]], bool)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_erosion(data, struct, border_value=1,
                               iterations=1, output=out, origin=(-1, -1))
        assert_array_almost_equal(out, expected)

    def test_binary_erosion32(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_erosion(data, struct,
                                     border_value=1, iterations=2)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion33(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 1, 1],
                    [0, 0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        mask = [[1, 1, 1, 1, 1, 0, 0],
                [1, 1, 1, 1, 1, 1, 0],
                [1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1, 1, 1, 1]]
        data = numpy.array([[0, 0, 0, 0, 0, 1, 1],
                            [0, 0, 0, 1, 0, 0, 1],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_erosion(data, struct,
                                     border_value=1, mask=mask, iterations=-1)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion34(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 1, 1, 1, 1, 1, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0]]
        mask = [[0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 1, 0, 1, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_erosion(data, struct,
                                     border_value=1, mask=mask)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion35(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        mask = [[0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 1, 0, 1, 0, 0],
                [0, 0, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0]]
        data = numpy.array([[0, 0, 0, 1, 0, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0],
                            [1, 1, 1, 1, 1, 1, 1],
                            [0, 1, 1, 1, 1, 1, 0],
                            [0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 1, 0, 0, 0]], bool)
        tmp = [[0, 0, 1, 0, 0, 0, 0],
               [0, 1, 1, 1, 0, 0, 0],
               [1, 1, 1, 1, 1, 0, 1],
               [0, 1, 1, 1, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 1]]
        expected = numpy.logical_and(tmp, mask)
        tmp = numpy.logical_and(data, numpy.logical_not(mask))
        expected = numpy.logical_or(expected, tmp)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_erosion(data, struct, border_value=1,
                               iterations=1, output=out,
                               origin=(-1, -1), mask=mask)
        assert_array_almost_equal(out, expected)

    def test_binary_erosion36(self):
        struct = [[0, 1, 0],
                  [1, 0, 1],
                  [0, 1, 0]]
        mask = [[0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0, 0],
                [0, 0, 1, 0, 1, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0, 0],
                [0, 0, 1, 1, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0]]
        tmp = [[0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 0, 1, 0, 0, 1],
               [0, 0, 1, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 1]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 1, 1],
                            [0, 0, 1, 1, 1, 0, 1, 1],
                            [0, 0, 1, 0, 1, 1, 0, 0],
                            [0, 1, 0, 1, 1, 1, 1, 0],
                            [0, 1, 1, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]])
        expected = numpy.logical_and(tmp, mask)
        tmp = numpy.logical_and(data, numpy.logical_not(mask))
        expected = numpy.logical_or(expected, tmp)
        out = ndimage.binary_erosion(data, struct, mask=mask,
                                     border_value=1, origin=(-1, -1))
        assert_array_almost_equal(out, expected)

    def test_binary_erosion37(self):
        a = numpy.array([[1, 0, 1],
                         [0, 1, 0],
                         [1, 0, 1]], dtype=bool)
        b = numpy.zeros_like(a)
        out = ndimage.binary_erosion(a, structure=a, output=b, iterations=0,
                                     border_value=True, brute_force=True)
        assert_(out is b)
        assert_array_equal(
            ndimage.binary_erosion(a, structure=a, iterations=0,
                                   border_value=True),
            b)

    def test_binary_dilation01(self):
        for type_ in self.types:
            data = numpy.ones([], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, 1)

    def test_binary_dilation02(self):
        for type_ in self.types:
            data = numpy.zeros([], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, 0)

    def test_binary_dilation03(self):
        for type_ in self.types:
            data = numpy.ones([1], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [1])

    def test_binary_dilation04(self):
        for type_ in self.types:
            data = numpy.zeros([1], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [0])

    def test_binary_dilation05(self):
        for type_ in self.types:
            data = numpy.ones([3], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [1, 1, 1])

    def test_binary_dilation06(self):
        for type_ in self.types:
            data = numpy.zeros([3], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [0, 0, 0])

    def test_binary_dilation07(self):
        for type_ in self.types:
            data = numpy.zeros([3], type_)
            data[1] = 1
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [1, 1, 1])

    def test_binary_dilation08(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            data[3] = 1
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [1, 1, 1, 1, 1])

    def test_binary_dilation09(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [1, 1, 1, 0, 0])

    def test_binary_dilation10(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            out = ndimage.binary_dilation(data, origin=-1)
            assert_array_almost_equal(out, [0, 1, 1, 1, 0])

    def test_binary_dilation11(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            out = ndimage.binary_dilation(data, origin=1)
            assert_array_almost_equal(out, [1, 1, 0, 0, 0])

    def test_binary_dilation12(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            struct = [1, 0, 1]
            out = ndimage.binary_dilation(data, struct)
            assert_array_almost_equal(out, [1, 0, 1, 0, 0])

    def test_binary_dilation13(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            struct = [1, 0, 1]
            out = ndimage.binary_dilation(data, struct, border_value=1)
            assert_array_almost_equal(out, [1, 0, 1, 0, 1])

    def test_binary_dilation14(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            struct = [1, 0, 1]
            out = ndimage.binary_dilation(data, struct, origin=-1)
            assert_array_almost_equal(out, [0, 1, 0, 1, 0])

    def test_binary_dilation15(self):
        for type_ in self.types:
            data = numpy.zeros([5], type_)
            data[1] = 1
            struct = [1, 0, 1]
            out = ndimage.binary_dilation(data, struct,
                                          origin=-1, border_value=1)
            assert_array_almost_equal(out, [1, 1, 0, 1, 0])

    def test_binary_dilation16(self):
        for type_ in self.types:
            data = numpy.ones([1, 1], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [[1]])

    def test_binary_dilation17(self):
        for type_ in self.types:
            data = numpy.zeros([1, 1], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [[0]])

    def test_binary_dilation18(self):
        for type_ in self.types:
            data = numpy.ones([1, 3], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [[1, 1, 1]])

    def test_binary_dilation19(self):
        for type_ in self.types:
            data = numpy.ones([3, 3], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [[1, 1, 1],
                                            [1, 1, 1],
                                            [1, 1, 1]])

    def test_binary_dilation20(self):
        for type_ in self.types:
            data = numpy.zeros([3, 3], type_)
            data[1, 1] = 1
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, [[0, 1, 0],
                                            [1, 1, 1],
                                            [0, 1, 0]])

    def test_binary_dilation21(self):
        struct = ndimage.generate_binary_structure(2, 2)
        for type_ in self.types:
            data = numpy.zeros([3, 3], type_)
            data[1, 1] = 1
            out = ndimage.binary_dilation(data, struct)
            assert_array_almost_equal(out, [[1, 1, 1],
                                            [1, 1, 1],
                                            [1, 1, 1]])

    def test_binary_dilation22(self):
        expected = [[0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 1, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation23(self):
        expected = [[1, 1, 1, 1, 1, 1, 1, 1],
                    [1, 1, 1, 0, 0, 0, 0, 1],
                    [1, 1, 0, 0, 0, 1, 0, 1],
                    [1, 0, 0, 1, 1, 1, 1, 1],
                    [1, 0, 1, 1, 1, 1, 0, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1],
                    [1, 0, 1, 0, 0, 1, 0, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation24(self):
        expected = [[1, 1, 0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 1, 0, 0, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [0, 1, 1, 1, 1, 0, 0, 0],
                    [1, 1, 1, 1, 1, 1, 0, 0],
                    [0, 1, 0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, origin=(1, 1))
            assert_array_almost_equal(out, expected)

    def test_binary_dilation25(self):
        expected = [[1, 1, 0, 0, 0, 0, 1, 1],
                    [1, 0, 0, 0, 1, 0, 1, 1],
                    [0, 0, 1, 1, 1, 1, 1, 1],
                    [0, 1, 1, 1, 1, 0, 1, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1],
                    [0, 1, 0, 0, 1, 0, 1, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, origin=(1, 1), border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation26(self):
        struct = ndimage.generate_binary_structure(2, 2)
        expected = [[1, 1, 1, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 1, 1, 1, 0],
                    [0, 0, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, struct)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation27(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 1, 1, 1, 0, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [0, 1, 1, 0, 1, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, struct)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation28(self):
        expected = [[1, 1, 1, 1],
                    [1, 0, 0, 1],
                    [1, 0, 0, 1],
                    [1, 1, 1, 1]]

        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_dilation29(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0]]

        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_dilation(data, struct, iterations=2)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation30(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0]]

        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0]], bool)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_dilation(data, struct, iterations=2, output=out)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation31(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 0, 0, 1, 0],
                    [0, 0, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [1, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0]]

        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_dilation(data, struct, iterations=3)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation32(self):
        struct = [[0, 1],
                  [1, 1]]
        expected = [[0, 0, 0, 1, 0],
                    [0, 0, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [1, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0]]

        data = numpy.array([[0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 0]], bool)
        out = numpy.zeros(data.shape, bool)
        ndimage.binary_dilation(data, struct, iterations=3, output=out)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation33(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 0, 0, 0],
                                [0, 1, 1, 0, 1, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        mask = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0, 0],
                            [0, 1, 1, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)

        out = ndimage.binary_dilation(data, struct, iterations=-1,
                                      mask=mask, border_value=0)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation34(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        mask = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.zeros(mask.shape, bool)
        out = ndimage.binary_dilation(data, struct, iterations=-1,
                                      mask=mask, border_value=1)
        assert_array_almost_equal(out, expected)

    def test_binary_dilation35(self):
        tmp = [[1, 1, 0, 0, 0, 0, 1, 1],
               [1, 0, 0, 0, 1, 0, 1, 1],
               [0, 0, 1, 1, 1, 1, 1, 1],
               [0, 1, 1, 1, 1, 0, 1, 1],
               [1, 1, 1, 1, 1, 1, 1, 1],
               [0, 1, 0, 0, 1, 0, 1, 1],
               [1, 1, 1, 1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1, 1, 1, 1]]
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]])
        mask = [[0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 1, 0, 0],
                [0, 0, 1, 1, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0]]
        expected = numpy.logical_and(tmp, mask)
        tmp = numpy.logical_and(data, numpy.logical_not(mask))
        expected = numpy.logical_or(expected, tmp)
        for type_ in self.types:
            data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_dilation(data, mask=mask,
                                          origin=(1, 1), border_value=1)
            assert_array_almost_equal(out, expected)

    def test_binary_propagation01(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 1, 1, 0, 0],
                               [0, 0, 1, 1, 1, 0, 0, 0],
                               [0, 1, 1, 0, 1, 1, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        mask = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 1, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 1, 1, 1, 0, 0, 0],
                            [0, 1, 1, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)

        out = ndimage.binary_propagation(data, struct,
                                         mask=mask, border_value=0)
        assert_array_almost_equal(out, expected)

    def test_binary_propagation02(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        mask = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.zeros(mask.shape, bool)
        out = ndimage.binary_propagation(data, struct,
                                         mask=mask, border_value=1)
        assert_array_almost_equal(out, expected)

    def test_binary_opening01(self):
        expected = [[0, 1, 0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 1, 0],
                                [0, 0, 1, 1, 0, 1, 0, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_opening(data)
            assert_array_almost_equal(out, expected)

    def test_binary_opening02(self):
        struct = ndimage.generate_binary_structure(2, 2)
        expected = [[1, 1, 1, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0],
                    [1, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[1, 1, 1, 0, 0, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0, 0, 0],
                                [1, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 1, 0, 1, 1, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_opening(data, struct)
            assert_array_almost_equal(out, expected)

    def test_binary_closing01(self):
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 1, 1, 1, 0, 1, 0, 0],
                    [0, 0, 1, 1, 1, 1, 1, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 1, 0, 0, 0, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 1, 1, 1, 0],
                                [0, 0, 1, 1, 0, 1, 0, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_closing(data)
            assert_array_almost_equal(out, expected)

    def test_binary_closing02(self):
        struct = ndimage.generate_binary_structure(2, 2)
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 1, 1, 1, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[1, 1, 1, 0, 0, 0, 0, 0],
                                [1, 1, 1, 0, 0, 0, 0, 0],
                                [1, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 1, 0, 1, 1, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_closing(data, struct)
            assert_array_almost_equal(out, expected)

    def test_binary_fill_holes01(self):
        expected = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 1, 1, 1, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 1, 1, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_fill_holes(data)
        assert_array_almost_equal(out, expected)

    def test_binary_fill_holes02(self):
        expected = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 1, 1, 1, 1, 0, 0],
                                [0, 0, 0, 1, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 1, 0, 0, 1, 0, 0],
                            [0, 0, 0, 1, 1, 0, 0, 0],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_fill_holes(data)
        assert_array_almost_equal(out, expected)

    def test_binary_fill_holes03(self):
        expected = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 1, 0, 0, 0, 0, 0],
                                [0, 1, 1, 1, 0, 1, 1, 1],
                                [0, 1, 1, 1, 0, 1, 1, 1],
                                [0, 1, 1, 1, 0, 1, 1, 1],
                                [0, 0, 1, 0, 0, 1, 1, 1],
                                [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        data = numpy.array([[0, 0, 0, 0, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0, 0, 0],
                            [0, 1, 0, 1, 0, 1, 1, 1],
                            [0, 1, 0, 1, 0, 1, 0, 1],
                            [0, 1, 0, 1, 0, 1, 0, 1],
                            [0, 0, 1, 0, 0, 1, 1, 1],
                            [0, 0, 0, 0, 0, 0, 0, 0]], bool)
        out = ndimage.binary_fill_holes(data)
        assert_array_almost_equal(out, expected)

    def test_grey_erosion01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        output = ndimage.grey_erosion(array, footprint=footprint)
        assert_array_almost_equal([[2, 2, 1, 1, 1],
                                   [2, 3, 1, 3, 1],
                                   [5, 5, 3, 3, 1]], output)

    def test_grey_erosion02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        output = ndimage.grey_erosion(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal([[2, 2, 1, 1, 1],
                                   [2, 3, 1, 3, 1],
                                   [5, 5, 3, 3, 1]], output)

    def test_grey_erosion03(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[1, 1, 1], [1, 1, 1]]
        output = ndimage.grey_erosion(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal([[1, 1, 0, 0, 0],
                                   [1, 2, 0, 2, 0],
                                   [4, 4, 2, 2, 0]], output)

    def test_grey_dilation01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[0, 1, 1], [1, 0, 1]]
        output = ndimage.grey_dilation(array, footprint=footprint)
        assert_array_almost_equal([[7, 7, 9, 9, 5],
                                   [7, 9, 8, 9, 7],
                                   [8, 8, 8, 7, 7]], output)

    def test_grey_dilation02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[0, 1, 1], [1, 0, 1]]
        structure = [[0, 0, 0], [0, 0, 0]]
        output = ndimage.grey_dilation(array, footprint=footprint,
                                       structure=structure)
        assert_array_almost_equal([[7, 7, 9, 9, 5],
                                   [7, 9, 8, 9, 7],
                                   [8, 8, 8, 7, 7]], output)

    def test_grey_dilation03(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[0, 1, 1], [1, 0, 1]]
        structure = [[1, 1, 1], [1, 1, 1]]
        output = ndimage.grey_dilation(array, footprint=footprint,
                                       structure=structure)
        assert_array_almost_equal([[8, 8, 10, 10, 6],
                                   [8, 10, 9, 10, 8],
                                   [9, 9, 9, 8, 8]], output)

    def test_grey_opening01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        tmp = ndimage.grey_erosion(array, footprint=footprint)
        expected = ndimage.grey_dilation(tmp, footprint=footprint)
        output = ndimage.grey_opening(array, footprint=footprint)
        assert_array_almost_equal(expected, output)

    def test_grey_opening02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_erosion(array, footprint=footprint,
                                   structure=structure)
        expected = ndimage.grey_dilation(tmp, footprint=footprint,
                                         structure=structure)
        output = ndimage.grey_opening(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal(expected, output)

    def test_grey_closing01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        tmp = ndimage.grey_dilation(array, footprint=footprint)
        expected = ndimage.grey_erosion(tmp, footprint=footprint)
        output = ndimage.grey_closing(array, footprint=footprint)
        assert_array_almost_equal(expected, output)

    def test_grey_closing02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_dilation(array, footprint=footprint,
                                    structure=structure)
        expected = ndimage.grey_erosion(tmp, footprint=footprint,
                                        structure=structure)
        output = ndimage.grey_closing(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal(expected, output)

    def test_morphological_gradient01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp1 = ndimage.grey_dilation(array, footprint=footprint,
                                     structure=structure)
        tmp2 = ndimage.grey_erosion(array, footprint=footprint,
                                    structure=structure)
        expected = tmp1 - tmp2
        output = numpy.zeros(array.shape, array.dtype)
        ndimage.morphological_gradient(array, footprint=footprint,
                                       structure=structure, output=output)
        assert_array_almost_equal(expected, output)

    def test_morphological_gradient02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp1 = ndimage.grey_dilation(array, footprint=footprint,
                                     structure=structure)
        tmp2 = ndimage.grey_erosion(array, footprint=footprint,
                                    structure=structure)
        expected = tmp1 - tmp2
        output = ndimage.morphological_gradient(array, footprint=footprint,
                                                structure=structure)
        assert_array_almost_equal(expected, output)

    def test_morphological_laplace01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp1 = ndimage.grey_dilation(array, footprint=footprint,
                                     structure=structure)
        tmp2 = ndimage.grey_erosion(array, footprint=footprint,
                                    structure=structure)
        expected = tmp1 + tmp2 - 2 * array
        output = numpy.zeros(array.shape, array.dtype)
        ndimage.morphological_laplace(array, footprint=footprint,
                                      structure=structure, output=output)
        assert_array_almost_equal(expected, output)

    def test_morphological_laplace02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp1 = ndimage.grey_dilation(array, footprint=footprint,
                                     structure=structure)
        tmp2 = ndimage.grey_erosion(array, footprint=footprint,
                                    structure=structure)
        expected = tmp1 + tmp2 - 2 * array
        output = ndimage.morphological_laplace(array, footprint=footprint,
                                               structure=structure)
        assert_array_almost_equal(expected, output)

    def test_white_tophat01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_opening(array, footprint=footprint,
                                   structure=structure)
        expected = array - tmp
        output = numpy.zeros(array.shape, array.dtype)
        ndimage.white_tophat(array, footprint=footprint,
                             structure=structure, output=output)
        assert_array_almost_equal(expected, output)

    def test_white_tophat02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_opening(array, footprint=footprint,
                                   structure=structure)
        expected = array - tmp
        output = ndimage.white_tophat(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal(expected, output)

    def test_white_tophat03(self):
        array = numpy.array([[1, 0, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 0, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 0, 1]], dtype=numpy.bool_)
        structure = numpy.ones((3, 3), dtype=numpy.bool_)
        expected = numpy.array([[0, 1, 1, 0, 0, 0, 0],
                                [1, 0, 0, 1, 1, 1, 0],
                                [1, 0, 0, 1, 1, 1, 0],
                                [0, 1, 1, 0, 0, 0, 1],
                                [0, 1, 1, 0, 1, 0, 1],
                                [0, 1, 1, 0, 0, 0, 1],
                                [0, 0, 0, 1, 1, 1, 1]], dtype=numpy.bool_)

        output = ndimage.white_tophat(array, structure=structure)
        assert_array_equal(expected, output)

    def test_white_tophat04(self):
        array = numpy.eye(5, dtype=numpy.bool_)
        structure = numpy.ones((3, 3), dtype=numpy.bool_)

        # Check that type missmatch is properly handled
        output = numpy.empty_like(array, dtype=numpy.float)
        ndimage.white_tophat(array, structure=structure, output=output)

    def test_black_tophat01(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_closing(array, footprint=footprint,
                                   structure=structure)
        expected = tmp - array
        output = numpy.zeros(array.shape, array.dtype)
        ndimage.black_tophat(array, footprint=footprint,
                             structure=structure, output=output)
        assert_array_almost_equal(expected, output)

    def test_black_tophat02(self):
        array = numpy.array([[3, 2, 5, 1, 4],
                             [7, 6, 9, 3, 5],
                             [5, 8, 3, 7, 1]])
        footprint = [[1, 0, 1], [1, 1, 0]]
        structure = [[0, 0, 0], [0, 0, 0]]
        tmp = ndimage.grey_closing(array, footprint=footprint,
                                   structure=structure)
        expected = tmp - array
        output = ndimage.black_tophat(array, footprint=footprint,
                                      structure=structure)
        assert_array_almost_equal(expected, output)

    def test_black_tophat03(self):
        array = numpy.array([[1, 0, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 1, 1, 1, 0, 1, 0],
                             [0, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 0, 1]], dtype=numpy.bool_)
        structure = numpy.ones((3, 3), dtype=numpy.bool_)
        expected = numpy.array([[0, 1, 1, 1, 1, 1, 1],
                                [1, 0, 0, 0, 0, 0, 1],
                                [1, 0, 0, 0, 0, 0, 1],
                                [1, 0, 0, 0, 0, 0, 1],
                                [1, 0, 0, 0, 1, 0, 1],
                                [1, 0, 0, 0, 0, 0, 1],
                                [1, 1, 1, 1, 1, 1, 0]], dtype=numpy.bool_)

        output = ndimage.black_tophat(array, structure=structure)
        assert_array_equal(expected, output)

    def test_black_tophat04(self):
        array = numpy.eye(5, dtype=numpy.bool_)
        structure = numpy.ones((3, 3), dtype=numpy.bool_)

        # Check that type missmatch is properly handled
        output = numpy.empty_like(array, dtype=numpy.float)
        ndimage.black_tophat(array, structure=structure, output=output)

    def test_hit_or_miss01(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 1, 0, 0, 0],
                                [1, 1, 1, 0, 0],
                                [0, 1, 0, 1, 1],
                                [0, 0, 1, 1, 1],
                                [0, 1, 1, 1, 0],
                                [0, 1, 1, 1, 1],
                                [0, 1, 1, 1, 1],
                                [0, 0, 0, 0, 0]], type_)
            out = numpy.zeros(data.shape, bool)
            ndimage.binary_hit_or_miss(data, struct, output=out)
            assert_array_almost_equal(expected, out)

    def test_hit_or_miss02(self):
        struct = [[0, 1, 0],
                  [1, 1, 1],
                  [0, 1, 0]]
        expected = [[0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 1, 0, 0, 1, 1, 1, 0],
                                [1, 1, 1, 0, 0, 1, 0, 0],
                                [0, 1, 0, 1, 1, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_hit_or_miss(data, struct)
            assert_array_almost_equal(expected, out)

    def test_hit_or_miss03(self):
        struct1 = [[0, 0, 0],
                   [1, 1, 1],
                   [0, 0, 0]]
        struct2 = [[1, 1, 1],
                   [0, 0, 0],
                   [1, 1, 1]]
        expected = [[0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0]]
        for type_ in self.types:
            data = numpy.array([[0, 1, 0, 0, 1, 1, 1, 0],
                                [1, 1, 1, 0, 0, 0, 0, 0],
                                [0, 1, 0, 1, 1, 1, 1, 0],
                                [0, 0, 1, 1, 1, 1, 1, 0],
                                [0, 1, 1, 1, 0, 1, 1, 0],
                                [0, 0, 0, 0, 1, 1, 1, 0],
                                [0, 1, 1, 1, 1, 1, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0]], type_)
            out = ndimage.binary_hit_or_miss(data, struct1, struct2)
            assert_array_almost_equal(expected, out)


class TestDilateFix:

    def setup_method(self):
        # dilation related setup
        self.array = numpy.array([[0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 0],
                                  [0, 0, 0, 1, 0],
                                  [0, 0, 1, 1, 0],
                                  [0, 0, 0, 0, 0]], dtype=numpy.uint8)

        self.sq3x3 = numpy.ones((3, 3))
        dilated3x3 = ndimage.binary_dilation(self.array, structure=self.sq3x3)
        self.dilated3x3 = dilated3x3.view(numpy.uint8)

    def test_dilation_square_structure(self):
        result = ndimage.grey_dilation(self.array, structure=self.sq3x3)
        # +1 accounts for difference between grey and binary dilation
        assert_array_almost_equal(result, self.dilated3x3 + 1)

    def test_dilation_scalar_size(self):
        result = ndimage.grey_dilation(self.array, size=3)
        assert_array_almost_equal(result, self.dilated3x3)

class TestBinaryOpeningClosing:

    def setup_method(self):
        a = numpy.zeros((5,5), dtype=bool)
        a[1:4, 1:4] = True
        a[4,4] = True
        self.array = a
        self.sq3x3 = numpy.ones((3,3))
        self.opened_old = ndimage.binary_opening(self.array, self.sq3x3,
                                                 1, None, 0)
        self.closed_old = ndimage.binary_closing(self.array, self.sq3x3,
                                                 1, None, 0)

    def test_opening_new_arguments(self):
        opened_new = ndimage.binary_opening(self.array, self.sq3x3, 1, None,
                                            0, None, 0, False)
        assert_array_equal(opened_new, self.opened_old)

    def test_closing_new_arguments(self):
        closed_new = ndimage.binary_closing(self.array, self.sq3x3, 1, None,
                                            0, None, 0, False)
        assert_array_equal(closed_new, self.closed_old)
