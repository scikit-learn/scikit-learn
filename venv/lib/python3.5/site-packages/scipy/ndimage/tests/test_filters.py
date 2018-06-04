''' Some tests for filters '''
from __future__ import division, print_function, absolute_import

import sys
import numpy as np

from numpy.testing import (assert_equal, assert_allclose,
                           assert_array_equal, assert_almost_equal)
from pytest import raises as assert_raises

import scipy.ndimage as sndi
from scipy.ndimage.filters import _gaussian_kernel1d


def test_ticket_701():
    # Test generic filter sizes
    arr = np.arange(4).reshape((2,2))
    func = lambda x: np.min(x)
    res = sndi.generic_filter(arr, func, size=(1,1))
    # The following raises an error unless ticket 701 is fixed
    res2 = sndi.generic_filter(arr, func, size=1)
    assert_equal(res, res2)


def test_gh_5430():
    # At least one of these raises an error unless gh-5430 is
    # fixed. In py2k an int is implemented using a C long, so
    # which one fails depends on your system. In py3k there is only
    # one arbitrary precision integer type, so both should fail.
    sigma = np.int32(1)
    out = sndi._ni_support._normalize_sequence(sigma, 1)
    assert_equal(out, [sigma])
    sigma = np.int64(1)
    out = sndi._ni_support._normalize_sequence(sigma, 1)
    assert_equal(out, [sigma])
    # This worked before; make sure it still works
    sigma = 1
    out = sndi._ni_support._normalize_sequence(sigma, 1)
    assert_equal(out, [sigma])
    # This worked before; make sure it still works
    sigma = [1, 1]
    out = sndi._ni_support._normalize_sequence(sigma, 2)
    assert_equal(out, sigma)
    # Also include the OPs original example to make sure we fixed the issue
    x = np.random.normal(size=(256, 256))
    perlin = np.zeros_like(x)
    for i in 2**np.arange(6):
        perlin += sndi.filters.gaussian_filter(x, i, mode="wrap") * i**2
    # This also fixes gh-4106, show that the OPs example now runs.
    x = np.int64(21)
    sndi._ni_support._normalize_sequence(x, 0)


def test_gaussian_kernel1d():
    radius = 10
    sigma = 2
    sigma2 = sigma * sigma
    x = np.arange(-radius, radius + 1, dtype=np.double)
    phi_x = np.exp(-0.5 * x * x / sigma2)
    phi_x /= phi_x.sum()
    assert_allclose(phi_x, _gaussian_kernel1d(sigma, 0, radius))
    assert_allclose(-phi_x * x / sigma2, _gaussian_kernel1d(sigma, 1, radius))
    assert_allclose(phi_x * (x * x / sigma2 - 1) / sigma2,
                    _gaussian_kernel1d(sigma, 2, radius))
    assert_allclose(phi_x * (3 - x * x / sigma2) * x / (sigma2 * sigma2),
                    _gaussian_kernel1d(sigma, 3, radius))


def test_orders_gauss():
    # Check order inputs to Gaussians
    arr = np.zeros((1,))
    assert_equal(0, sndi.gaussian_filter(arr, 1, order=0))
    assert_equal(0, sndi.gaussian_filter(arr, 1, order=3))
    assert_raises(ValueError, sndi.gaussian_filter, arr, 1, -1)
    assert_equal(0, sndi.gaussian_filter1d(arr, 1, axis=-1, order=0))
    assert_equal(0, sndi.gaussian_filter1d(arr, 1, axis=-1, order=3))
    assert_raises(ValueError, sndi.gaussian_filter1d, arr, 1, -1, -1)


def test_valid_origins():
    """Regression test for #1311."""
    func = lambda x: np.mean(x)
    data = np.array([1,2,3,4,5], dtype=np.float64)
    assert_raises(ValueError, sndi.generic_filter, data, func, size=3,
                  origin=2)
    func2 = lambda x, y: np.mean(x + y)
    assert_raises(ValueError, sndi.generic_filter1d, data, func,
                  filter_size=3, origin=2)
    assert_raises(ValueError, sndi.percentile_filter, data, 0.2, size=3,
                  origin=2)

    for filter in [sndi.uniform_filter, sndi.minimum_filter,
                   sndi.maximum_filter, sndi.maximum_filter1d,
                   sndi.median_filter, sndi.minimum_filter1d]:
        # This should work, since for size == 3, the valid range for origin is
        # -1 to 1.
        list(filter(data, 3, origin=-1))
        list(filter(data, 3, origin=1))
        # Just check this raises an error instead of silently accepting or
        # segfaulting.
        assert_raises(ValueError, filter, data, 3, origin=2)


def test_multiple_modes():
    # Test that the filters with multiple mode cababilities for different
    # dimensions give the same result as applying a single mode.
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    mode1 = 'reflect'
    mode2 = ['reflect', 'reflect']

    assert_equal(sndi.gaussian_filter(arr, 1, mode=mode1),
                 sndi.gaussian_filter(arr, 1, mode=mode2))
    assert_equal(sndi.prewitt(arr, mode=mode1),
                 sndi.prewitt(arr, mode=mode2))
    assert_equal(sndi.sobel(arr, mode=mode1),
                 sndi.sobel(arr, mode=mode2))
    assert_equal(sndi.laplace(arr, mode=mode1),
                 sndi.laplace(arr, mode=mode2))
    assert_equal(sndi.gaussian_laplace(arr, 1, mode=mode1),
                 sndi.gaussian_laplace(arr, 1, mode=mode2))
    assert_equal(sndi.maximum_filter(arr, size=5, mode=mode1),
                 sndi.maximum_filter(arr, size=5, mode=mode2))
    assert_equal(sndi.minimum_filter(arr, size=5, mode=mode1),
                 sndi.minimum_filter(arr, size=5, mode=mode2))
    assert_equal(sndi.gaussian_gradient_magnitude(arr, 1, mode=mode1),
                 sndi.gaussian_gradient_magnitude(arr, 1, mode=mode2))
    assert_equal(sndi.uniform_filter(arr, 5, mode=mode1),
                 sndi.uniform_filter(arr, 5, mode=mode2))


def test_multiple_modes_sequentially():
    # Test that the filters with multiple mode cababilities for different
    # dimensions give the same result as applying the filters with
    # different modes sequentially
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    modes = ['reflect', 'wrap']

    expected = sndi.gaussian_filter1d(arr, 1, axis=0, mode=modes[0])
    expected = sndi.gaussian_filter1d(expected, 1, axis=1, mode=modes[1])
    assert_equal(expected,
                 sndi.gaussian_filter(arr, 1, mode=modes))

    expected = sndi.uniform_filter1d(arr, 5, axis=0, mode=modes[0])
    expected = sndi.uniform_filter1d(expected, 5, axis=1, mode=modes[1])
    assert_equal(expected,
                 sndi.uniform_filter(arr, 5, mode=modes))

    expected = sndi.maximum_filter1d(arr, size=5, axis=0, mode=modes[0])
    expected = sndi.maximum_filter1d(expected, size=5, axis=1, mode=modes[1])
    assert_equal(expected,
                 sndi.maximum_filter(arr, size=5, mode=modes))

    expected = sndi.minimum_filter1d(arr, size=5, axis=0, mode=modes[0])
    expected = sndi.minimum_filter1d(expected, size=5, axis=1, mode=modes[1])
    assert_equal(expected,
                 sndi.minimum_filter(arr, size=5, mode=modes))


def test_multiple_modes_prewitt():
    # Test prewitt filter for multiple extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[1., -3., 2.],
                         [1., -2., 1.],
                         [1., -1., 0.]])

    modes = ['reflect', 'wrap']

    assert_equal(expected,
                 sndi.prewitt(arr, mode=modes))


def test_multiple_modes_sobel():
    # Test sobel filter for multiple extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[1., -4., 3.],
                         [2., -3., 1.],
                         [1., -1., 0.]])

    modes = ['reflect', 'wrap']

    assert_equal(expected,
                 sndi.sobel(arr, mode=modes))


def test_multiple_modes_laplace():
    # Test laplace filter for multiple extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[-2., 2., 1.],
                         [-2., -3., 2.],
                         [1., 1., 0.]])

    modes = ['reflect', 'wrap']

    assert_equal(expected,
                 sndi.laplace(arr, mode=modes))


def test_multiple_modes_gaussian_laplace():
    # Test gaussian_laplace filter for multiple extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[-0.28438687, 0.01559809, 0.19773499],
                         [-0.36630503, -0.20069774, 0.07483620],
                         [0.15849176, 0.18495566, 0.21934094]])

    modes = ['reflect', 'wrap']

    assert_almost_equal(expected,
                        sndi.gaussian_laplace(arr, 1, mode=modes))


def test_multiple_modes_gaussian_gradient_magnitude():
    # Test gaussian_gradient_magnitude filter for multiple
    # extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[0.04928965, 0.09745625, 0.06405368],
                         [0.23056905, 0.14025305, 0.04550846],
                         [0.19894369, 0.14950060, 0.06796850]])

    modes = ['reflect', 'wrap']

    calculated = sndi.gaussian_gradient_magnitude(arr, 1, mode=modes)

    assert_almost_equal(expected, calculated)


def test_multiple_modes_uniform():
    # Test uniform filter for multiple extrapolation modes
    arr = np.array([[1., 0., 0.],
                    [1., 1., 0.],
                    [0., 0., 0.]])

    expected = np.array([[0.32, 0.40, 0.48],
                         [0.20, 0.28, 0.32],
                         [0.28, 0.32, 0.40]])

    modes = ['reflect', 'wrap']

    assert_almost_equal(expected,
                        sndi.uniform_filter(arr, 5, mode=modes))


def test_gaussian_truncate():
    # Test that Gaussian filters can be truncated at different widths.
    # These tests only check that the result has the expected number
    # of nonzero elements.
    arr = np.zeros((100, 100), float)
    arr[50, 50] = 1
    num_nonzeros_2 = (sndi.gaussian_filter(arr, 5, truncate=2) > 0).sum()
    assert_equal(num_nonzeros_2, 21**2)
    num_nonzeros_5 = (sndi.gaussian_filter(arr, 5, truncate=5) > 0).sum()
    assert_equal(num_nonzeros_5, 51**2)

    # Test truncate when sigma is a sequence.
    f = sndi.gaussian_filter(arr, [0.5, 2.5], truncate=3.5)
    fpos = f > 0
    n0 = fpos.any(axis=0).sum()
    # n0 should be 2*int(2.5*3.5 + 0.5) + 1
    assert_equal(n0, 19)
    n1 = fpos.any(axis=1).sum()
    # n1 should be 2*int(0.5*3.5 + 0.5) + 1
    assert_equal(n1, 5)

    # Test gaussian_filter1d.
    x = np.zeros(51)
    x[25] = 1
    f = sndi.gaussian_filter1d(x, sigma=2, truncate=3.5)
    n = (f > 0).sum()
    assert_equal(n, 15)

    # Test gaussian_laplace
    y = sndi.gaussian_laplace(x, sigma=2, truncate=3.5)
    nonzero_indices = np.where(y != 0)[0]
    n = nonzero_indices.ptp() + 1
    assert_equal(n, 15)

    # Test gaussian_gradient_magnitude
    y = sndi.gaussian_gradient_magnitude(x, sigma=2, truncate=3.5)
    nonzero_indices = np.where(y != 0)[0]
    n = nonzero_indices.ptp() + 1
    assert_equal(n, 15)


class TestThreading(object):
    def check_func_thread(self, n, fun, args, out):
        from threading import Thread
        thrds = [Thread(target=fun, args=args, kwargs={'output': out[x]}) for x in range(n)]
        [t.start() for t in thrds]
        [t.join() for t in thrds]

    def check_func_serial(self, n, fun, args, out):
        for i in range(n):
            fun(*args, output=out[i])

    def test_correlate1d(self):
        d = np.random.randn(5000)
        os = np.empty((4, d.size))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.correlate1d, (d, np.arange(5)), os)
        self.check_func_thread(4, sndi.correlate1d, (d, np.arange(5)), ot)
        assert_array_equal(os, ot)

    def test_correlate(self):
        d = np.random.randn(500, 500)
        k = np.random.randn(10, 10)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.correlate, (d, k), os)
        self.check_func_thread(4, sndi.correlate, (d, k), ot)
        assert_array_equal(os, ot)

    def test_median_filter(self):
        d = np.random.randn(500, 500)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.median_filter, (d, 3), os)
        self.check_func_thread(4, sndi.median_filter, (d, 3), ot)
        assert_array_equal(os, ot)

    def test_uniform_filter1d(self):
        d = np.random.randn(5000)
        os = np.empty((4, d.size))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.uniform_filter1d, (d, 5), os)
        self.check_func_thread(4, sndi.uniform_filter1d, (d, 5), ot)
        assert_array_equal(os, ot)

    def test_minmax_filter(self):
        d = np.random.randn(500, 500)
        os = np.empty([4] + list(d.shape))
        ot = np.empty_like(os)
        self.check_func_serial(4, sndi.maximum_filter, (d, 3), os)
        self.check_func_thread(4, sndi.maximum_filter, (d, 3), ot)
        assert_array_equal(os, ot)
        self.check_func_serial(4, sndi.minimum_filter, (d, 3), os)
        self.check_func_thread(4, sndi.minimum_filter, (d, 3), ot)
        assert_array_equal(os, ot)


def test_minmaximum_filter1d():
    # Regression gh-3898
    in_ = np.arange(10)
    out = sndi.minimum_filter1d(in_, 1)
    assert_equal(in_, out)
    out = sndi.maximum_filter1d(in_, 1)
    assert_equal(in_, out)
    # Test reflect
    out = sndi.minimum_filter1d(in_, 5, mode='reflect')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 6, 7], out)
    out = sndi.maximum_filter1d(in_, 5, mode='reflect')
    assert_equal([2, 3, 4, 5, 6, 7, 8, 9, 9, 9], out)
    #Test constant
    out = sndi.minimum_filter1d(in_, 5, mode='constant', cval=-1)
    assert_equal([-1, -1, 0, 1, 2, 3, 4, 5, -1, -1], out)
    out = sndi.maximum_filter1d(in_, 5, mode='constant', cval=10)
    assert_equal([10, 10, 4, 5, 6, 7, 8, 9, 10, 10], out)
    # Test nearest
    out = sndi.minimum_filter1d(in_, 5, mode='nearest')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 6, 7], out)
    out = sndi.maximum_filter1d(in_, 5, mode='nearest')
    assert_equal([2, 3, 4, 5, 6, 7, 8, 9, 9, 9], out)
    # Test wrap
    out = sndi.minimum_filter1d(in_, 5, mode='wrap')
    assert_equal([0, 0, 0, 1, 2, 3, 4, 5, 0, 0], out)
    out = sndi.maximum_filter1d(in_, 5, mode='wrap')
    assert_equal([9, 9, 4, 5, 6, 7, 8, 9, 9, 9], out)


def test_uniform_filter1d_roundoff_errors():
    # gh-6930
    in_ = np.repeat([0, 1, 0], [9, 9, 9])
    for filter_size in range(3, 10):
        out = sndi.uniform_filter1d(in_, filter_size)
        assert_equal(out.sum(), 10 - filter_size)


def test_footprint_all_zeros():
    # regression test for gh-6876: footprint of all zeros segfaults
    arr = np.random.randint(0, 100, (100, 100))
    kernel = np.zeros((3, 3), bool)
    with assert_raises(ValueError):
        sndi.maximum_filter(arr, footprint=kernel)

def test_gaussian_filter():
    # Test gaussian filter with np.float16
    # gh-8207
    data = np.array([1],dtype = np.float16)
    sigma = 1.0
    with assert_raises(RuntimeError):
        sndi.gaussian_filter(data,sigma)
