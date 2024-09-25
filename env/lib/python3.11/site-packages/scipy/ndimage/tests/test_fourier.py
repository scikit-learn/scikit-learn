import numpy as np
from numpy import fft
from numpy.testing import (assert_almost_equal, assert_array_almost_equal,
                           assert_equal)

import pytest

from scipy import ndimage


class TestNdimageFourier:

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [(np.float32, 6), (np.float64, 14)])
    def test_fourier_gaussian_real01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_gaussian(a, [5.0, 2.5], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), 1, decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [(np.complex64, 6), (np.complex128, 14)])
    def test_fourier_gaussian_complex01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_gaussian(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [(np.float32, 6), (np.float64, 14)])
    def test_fourier_uniform_real01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_uniform(a, [5.0, 2.5], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), 1.0, decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [(np.complex64, 6), (np.complex128, 14)])
    def test_fourier_uniform_complex01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_uniform(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [(np.float32, 4), (np.float64, 11)])
    def test_fourier_shift_real01(self, shape, dtype, dec):
        expected = np.arange(shape[0] * shape[1], dtype=dtype)
        expected.shape = shape
        a = fft.rfft(expected, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_shift(a, [1, 1], shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_array_almost_equal(a[1:, 1:], expected[:-1, :-1], decimal=dec)
        assert_array_almost_equal(a.imag, np.zeros(shape), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [(np.complex64, 4), (np.complex128, 11)])
    def test_fourier_shift_complex01(self, shape, dtype, dec):
        expected = np.arange(shape[0] * shape[1], dtype=dtype)
        expected.shape = shape
        a = fft.fft(expected, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_shift(a, [1, 1], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_array_almost_equal(a.real[1:, 1:], expected[:-1, :-1], decimal=dec)
        assert_array_almost_equal(a.imag, np.zeros(shape), decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15), (1, 10)])
    @pytest.mark.parametrize('dtype, dec', [(np.float32, 5), (np.float64, 14)])
    def test_fourier_ellipsoid_real01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.rfft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_ellipsoid(a, [5.0, 2.5],
                                      shape[0], 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.irfft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a), 1.0, decimal=dec)

    @pytest.mark.parametrize('shape', [(32, 16), (31, 15)])
    @pytest.mark.parametrize('dtype, dec', [(np.complex64, 5), (np.complex128, 14)])
    def test_fourier_ellipsoid_complex01(self, shape, dtype, dec):
        a = np.zeros(shape, dtype)
        a[0, 0] = 1.0
        a = fft.fft(a, shape[0], 0)
        a = fft.fft(a, shape[1], 1)
        a = ndimage.fourier_ellipsoid(a, [5.0, 2.5], -1, 0)
        a = fft.ifft(a, shape[1], 1)
        a = fft.ifft(a, shape[0], 0)
        assert_almost_equal(ndimage.sum(a.real), 1.0, decimal=dec)

    def test_fourier_ellipsoid_unimplemented_ndim(self):
        # arrays with ndim > 3 raise NotImplementedError
        x = np.ones((4, 6, 8, 10), dtype=np.complex128)
        with pytest.raises(NotImplementedError):
            ndimage.fourier_ellipsoid(x, 3)

    def test_fourier_ellipsoid_1d_complex(self):
        # expected result of 1d ellipsoid is the same as for fourier_uniform
        for shape in [(32, ), (31, )]:
            for type_, dec in zip([np.complex64, np.complex128], [5, 14]):
                x = np.ones(shape, dtype=type_)
                a = ndimage.fourier_ellipsoid(x, 5, -1, 0)
                b = ndimage.fourier_uniform(x, 5, -1, 0)
                assert_array_almost_equal(a, b, decimal=dec)

    @pytest.mark.parametrize('shape', [(0, ), (0, 10), (10, 0)])
    @pytest.mark.parametrize('dtype', [np.float32, np.float64,
                                       np.complex64, np.complex128])
    @pytest.mark.parametrize('test_func',
                             [ndimage.fourier_ellipsoid,
                              ndimage.fourier_gaussian,
                              ndimage.fourier_uniform])
    def test_fourier_zero_length_dims(self, shape, dtype, test_func):
        a = np.ones(shape, dtype)
        b = test_func(a, 3)
        assert_equal(a, b)
